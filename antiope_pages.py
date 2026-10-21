#!/usr/bin/env python3
"""Assemble Antiope pages from Didon packets (output of antiope_demod.py).

Didon packet (CCIR 653 system A): byte 3 = E7 framing, bytes 4-6 channel
address (Hamming 8/4), byte 7 continuity index (Hamming), byte 8 format
(Hamming), bytes 9.. data block with odd parity.
Antiope pages are data groups SOH RS <record header> ... ETX EOT inside the
byte stream of one data channel.  Pages are carouselled, so all receptions
of a page are combined by majority vote per byte.
"""
import argparse
import collections
import html
import os
import re

HAM = [0x15, 0x02, 0x49, 0x5E, 0x64, 0x73, 0x38, 0x2F,
       0xD0, 0xC7, 0x8C, 0x9B, 0xA1, 0xB6, 0xFD, 0xEA]


def _ham_table():
    """Hamming 8/4 decode with single-bit correction (None = uncorrectable)."""
    t = {}
    for v in range(256):
        best = [(bin(v ^ c).count("1"), i) for i, c in enumerate(HAM)]
        best.sort()
        d, i = best[0]
        if d == 0 or (d == 1 and best[1][0] > 1):
            t[v] = i
        else:
            t[v] = None
    return t


HDEC = _ham_table()
# Data block length per format nibble F (TDF Didon spec ch. 1 annex 2, list 1;
# the list in use is a magazine convention, list 1 is the one seen on air).
# A 40-byte line holds at most 32 data bytes (F = B).
DATA_LEN = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 8, 6: 12, 7: 16, 8: 20, 9: 24, 10: 28, 11: 32}
MISSING = None            # placeholder for bytes of lost packets


def odd(b):
    return bin(b).count("1") & 1 == 1


def read_packets(path):
    for ln in open(path):
        f = ln.split()
        yield int(f[0]), int(f[1]), float(f[2]), bytes(int(x, 16) for x in f[3:])


def channel_streams(packets, max_gap_frames=2):
    """Yield per channel lists of (frame, data bytes or None) in continuity order.

    Addresses with one unreadable/odd nibble are matched to a known channel,
    an unreadable continuity index is taken as previous + 1.
    """
    packets = list(packets)
    cnt = collections.Counter()
    for fi, li, q, b in packets:
        a = tuple(HDEC[x] for x in b[1:4])
        if None not in a:
            cnt[a] += 1
    known = [a for a, k in cnt.items() if k >= max(20, 0.002 * len(packets))]
    streams = collections.defaultdict(list)
    last = {}
    for fi, li, q, b in packets:
        a = tuple(HDEC[x] for x in b[1:4])
        if a not in known:
            m = [k for k in known if sum(x == y for x, y in zip(a, k)) == 2]
            if len(m) != 1:
                continue
            a = m[0]
        fmt = HDEC[b[5]]
        if fmt is None:
            fmt = 11
        if fmt not in DATA_LEN:
            continue
        if fmt == 0:
            # empty (filler) packet: no data, and per spec ch. 1 sec. 5.3 a
            # continuity break next to it is not a loss
            if a in last:
                last[a] = (fi, None)
            continue
        ci = HDEC[b[4]]
        if ci is None:
            ci = (last[a][1] + 1) % 16 if a in last and last[a][1] is not None else None
            if ci is None:
                continue
        data = list(b[6:6 + DATA_LEN[fmt]])
        st = streams[a]
        if a in last:
            pf, pci = last[a]
            miss = 0 if pci is None else (ci - pci - 1) % 16
            if fi - pf > max_gap_frames:
                st.append((fi, "BREAK"))
            elif miss == 15:          # duplicate / repeated packet
                continue
            else:
                for _ in range(miss):
                    st.append((fi, None))
        st.append((fi, data))
        last[a] = (fi, ci)
    return streams


def flatten(st):
    """Stream entries -> list of segments (bytes, frame of each byte); None = lost byte."""
    segs, cur, fr = [], [], []
    for fi, d in st:
        if d == "BREAK":
            segs.append((cur, fr))
            cur, fr = [], []
        elif d is None:
            cur.extend([MISSING] * 32)
            fr.extend([fi] * 32)
        else:
            cur.extend(d)
            fr.extend([fi] * len(d))
    segs.append((cur, fr))
    return segs


def split_pages(segment):
    """Find SOH RS ... ETX EOT groups; return list of (bytes from SOH, first frame)."""
    seg, fr = segment
    out = []
    def is_start(i):
        # SOH RS, or RS + record header (3 Hamming address bytes, then 8) when SOH was lost
        if seg[i] is None or seg[i] & 0x7F != 0x1E:
            return False
        if i > 0 and seg[i - 1] is not None and seg[i - 1] & 0x7F == 0x01:
            return True
        h = seg[i + 1:i + 5]
        return len(h) == 4 and None not in h and all(x in HAM for x in h) and HDEC[h[3]] == 8

    starts = [i - 1 if i > 0 and seg[i - 1] is not None and seg[i - 1] & 0x7F == 0x01 else i
              for i in range(len(seg)) if is_start(i)]
    for k, s in enumerate(starts):
        end = starts[k + 1] if k + 1 < len(starts) else len(seg)
        body = seg[s:end]
        for i in range(len(body) - 1):
            if body[i] is not None and body[i + 1] is not None and \
               body[i] & 0x7F == 0x03 and body[i + 1] & 0x7F == 0x04:
                body = body[:i + 2]
                break
        else:
            # end marker lost: keep only up to the first lost packet, the
            # rest may already belong to the next group
            if None in body:
                body = body[:body.index(None)]
            elif k + 1 == len(starts):
                continue          # incomplete trailing group
            body = body + [MISSING] * 2
        if body[0] & 0x7F == 0x1E:
            body = [0x01] + body              # SOH was lost
        out.append((body, fr[s]))
    return out


def similarity(a, b):
    n = same = 0
    for x, y in zip(a, b):
        if x is None or y is None:
            continue
        n += 1
        same += x == y
    return same / n if n >= 32 else 0.0


def cluster(recs, thr=0.85):
    """Greedy clustering of receptions (body, frame) of one page number."""
    cl = []
    for body, fi in sorted(recs, key=lambda r: -sum(v is not None for v in r[0])):
        for c in cl:
            if abs(len(c["ref"]) - len(body)) <= 64 and similarity(c["ref"], body) >= thr:
                c["recs"].append((body, fi))
                break
        else:
            cl.append({"ref": body, "recs": [(body, fi)]})
    for c in cl:
        c["recs"].sort(key=lambda r: r[1])
    return sorted(cl, key=lambda c: -len(c["recs"]))


def page_header(body):
    """Record header after SOH RS: Hamming nibbles until first non-Hamming byte."""
    nib = []
    k = 2 if body and body[0] is not None and body[0] & 0x7F == 0x01 else 1
    for x in body[k:k + 18]:
        if x is None:
            nib.append(None)
            continue
        v = HDEC.get(x)
        if x not in HAM and (v is None or odd(x) and x & 0x7F >= 0x20 and len(nib) >= 9):
            break
        nib.append(v)
        if x == 0x1F or x == 0x9F:
            break
    return nib



# Defaults of the interpretation bytes (spec ch. 4 annex 2, figure 1) as nibble
# values, bit 0 = b2 ... bit 3 = b8.
Y_DEFAULTS = {"11": 0, "12": 0, "13": 0, "14": 8, "15": 0, "16": 8,
              "21": 0, "22": 0, "23": 0, "24": 7}


def parse_header(nib):
    """Article header after SOH RS (spec ch. 4 sec. 2.2).

    C1 C2 C3 (classification, decimal page number for 0-9), L (link byte:
    b2 b4 b6 = part number, b8 = last part), then designation bytes Y0i.  In
    Y0i, b4/b6/b8 announce the pairs (Yi1,Yi2)/(Yi3,Yi4)/(Yi5,Yi6) and b2 a
    further designation byte.  Returns a dict with C, L, Y (all Yij with
    defaults filled in) and n (number of header nibbles).
    """
    h = {"C": nib[:3], "L": nib[3] if len(nib) > 3 else None, "Y": dict(Y_DEFAULTS), "n": len(nib)}
    i, k = 4, 1
    while i < len(nib):
        y0 = nib[i]
        if y0 is None:
            break
        i += 1
        for bit, pair in ((2, 1), (4, 3), (8, 5)):
            if y0 & bit:
                for j in (pair, pair + 1):
                    if i < len(nib):
                        if nib[i] is not None:
                            h["Y"]["%d%d" % (k, j)] = nib[i]
                        i += 1
        k += 1
        if not y0 & 1:
            break
    h["n"] = i
    return h


def header_defaults(h):
    """Article-level ('level 2') attribute values from Y21..Y24."""
    y = h["Y"]
    return {
        "bg": (y["21"] & 7) if y["21"] & 8 else None,   # b8 = 0: transparent
        "conceal": bool(y["22"] & 2),
        "dw": bool(y["23"] & 1), "dh": bool(y["23"] & 2),
        "lining": bool(y["23"] & 4), "flash": bool(y["23"] & 8),
        "fg": y["24"] & 7, "inv": bool(y["24"] & 8),
    }

def vote(bodies):
    n = max(len(b) for b in bodies)
    res, conf = [], []
    for i in range(n):
        c = collections.Counter()
        for b in bodies:
            if i < len(b) and b[i] is not None:
                w = 2 if odd(b[i]) else 1
                c[b[i]] += w
        if not c:
            res.append(None)
            conf.append(0)
            continue
        v, k = c.most_common(1)[0]
        res.append(v)
        conf.append(k)
    return res, conf


# ---------------------------------------------------------------- rendering
# G2 (spec ch. 3 fig. 5): column 4 = non-spacing accents, the rest symbols.
G2_ACCENT = {0x41: "\u0300", 0x42: "\u0301", 0x43: "\u0302", 0x44: "\u0303",
             0x45: "\u0304", 0x46: "\u0306", 0x47: "\u0307", 0x48: "\u0308",
             0x4A: "\u030A", 0x4B: "\u0327", 0x4D: "\u030B", 0x4E: "\u0328",
             0x4F: "\u030C"}
G2_CHARS = dict(zip(range(0x21, 0x40), "¡¢£$¥#§¤‘“«←↑→↓°±²³×µ¶·÷’”»¼½¾¿"))
G2_CHARS.update({0x50: "\u2015"})
G2_CHARS.update(dict(zip(range(0x60, 0x7F), "ΩÆĐªĦ\ufffdĲĿŁØŒºÞŦŊŉĸæđðħıĳŀłøœßþŧŋ")))
# G0 positions that differ from ASCII (spec ch. 3 fig. 3 and ch. 4 annex 4):
# 6/0 is a horizontal bar, 7/B and 7/D are drawn as vertical bars at the left
# and right cell edge (fallback for braces), 7/E is an overline.
G0_SPECIAL = {0x60: "\u2500", 0x7B: "\u258f", 0x7D: "\u2595", 0x7E: "\u203e"}
COLORS = ["#000", "#f00", "#0f0", "#ff0", "#00f", "#f0f", "#0ff", "#fff"]
ROWS, COLS = 25, 40


def mosaic(c):
    """G1 2x3 block mosaic -> Unicode (Symbols for Legacy Computing)."""
    if c == 0x20:
        return " "
    if c == 0x7F:
        return "█"
    bits = (c & 0x1F) | ((c & 0x40) >> 1)
    if bits == 0:
        return " "
    if bits == 0x3F:
        return "█"
    if bits == 0x15:
        return "▌"
    if bits == 0x2A:
        return "▐"
    idx = bits - 1
    if bits > 0x15:
        idx -= 1
    if bits > 0x2A:
        idx -= 1
    return chr(0x1FB00 + idx)


TRANSPARENT = 8                        # background index for "fond transparent"
LEVEL3 = {"bg": 0, "conceal": False, "dw": False, "dh": False, "lining": False,
          "flash": False, "fg": 7, "inv": False}   # spec ch. 3 sec. 4.1 defaults


class Screen:
    """Attribute levels follow spec ch. 3 sec. 4.1: US and FF fall back to the
    article level (header Y21..Y24), RS to the highest-level defaults."""

    def __init__(self, defaults=None):
        self.level2 = dict(LEVEL3, **{k: v for k, v in (defaults or {}).items()})
        if self.level2["bg"] is None:
            self.level2["bg"] = TRANSPARENT
        d = self.level2
        self.ch = [[" "] * COLS for _ in range(ROWS)]
        self.fg = [[d["fg"]] * COLS for _ in range(ROWS)]
        self.bg = [[d["bg"]] * COLS for _ in range(ROWS)]
        self.fg[0] = [7] * COLS               # row 0 = status row outside the page area
        self.bg[0] = [0] * COLS
        self.dh = [[False] * COLS for _ in range(ROWS)]
        self.sz = [["L"] * COLS for _ in range(ROWS)]
        self.mos = [[None] * COLS for _ in range(ROWS)]
        self.hid = [[False] * COLS for _ in range(ROWS)]
        self.sepm = [[False] * COLS for _ in range(ROWS)]
        self.fl = [[False] * COLS for _ in range(ROWS)]
        self.r, self.c = 1, 0
        self.g1 = False
        self.attrs_reset()
        self.last, self.last_m = " ", None

    def put(self, s, mbits=None):
        if 0 <= self.r < ROWS and 0 <= self.c < COLS:
            self.ch[self.r][self.c] = s
            self.mos[self.r][self.c] = mbits
            self.hid[self.r][self.c] = self.conceal
            self.sepm[self.r][self.c] = self.sep
            self.fl[self.r][self.c] = self.flash
            # size and polarity do not apply to mosaics (spec ch. 2 sec. 5.2)
            self.sz[self.r][self.c] = self.size if mbits is None else "L"
            inv = self.inv and mbits is None
            f, b = (self.cbg, self.cfg) if inv else (self.cfg, self.cbg)
            self.fg[self.r][self.c] = f
            self.bg[self.r][self.c] = b
            self.dh[self.r][self.c] = self.size in "MO" and mbits is None
        self.last, self.last_m = s, mbits
        self.c += 2 if self.size in "NO" and mbits is None else 1
        if self.c >= COLS:
            self.c = 0
            self.r += 1

    def set_level(self, d):
        self.cfg, self.cbg, self.inv = d["fg"], d["bg"], d["inv"]
        self.size = {(False, False): "L", (False, True): "M", (True, False): "N", (True, True): "O"}[
            (d["dw"], d["dh"])]
        self.conceal, self.sep, self.flash = d["conceal"], d["lining"], d["flash"]

    def attrs_reset(self):
        """US / FF: back to the article level (header defaults)."""
        self.set_level(self.level2)
        self.g1 = False                       # US designates G0 (spec ch. 3 sec. 6.4)
        if self.r == 0:
            self.cfg, self.cbg = 7, 0


def render(body, start, defaults=None):
    """Interpret Antiope presentation bytes (parity stripped) on a 25x40 screen.

    defaults: article-level attributes from header_defaults().
    """
    s = Screen(defaults)
    data = [None if x is None else x & 0x7F for x in body[start:]]
    i = 0
    n = len(data)

    def nxt(k=1):
        return data[i + k] if i + k < n else None

    # Antiope sends positioning and mode controls twice (US r c US r c,
    # SO SO, ESC x ESC x); drop the repetition.  REP pairs are real (both count).
    seqlen = {0x1F: 3, 0x1B: 2, 0x19: 2, 0x12: 2}
    out = []
    k = 0
    while k < n:
        x = data[k]
        if x is not None and x < 0x20:
            L = seqlen.get(x, 1)
            seq = data[k:k + L]
            if data[k + L:k + 2 * L] == seq and x not in (0x19, 0x12):
                out.extend(seq)
                k += 2 * L
                continue
            # duplicate with one copy damaged/lost: keep the intact copy
        out.append(x)
        k += 1
    data = out
    n = len(data)
    while i < n:
        x = data[i]
        if x is None:
            s.put("·")
            i += 1
            continue
        if x == 0x1F:                         # US row col
            r, c = nxt(1), nxt(2)
            if r is not None and c is not None and 0x40 <= r < 0x40 + ROWS + 8:
                s.r, s.c = r - 0x40, max(0, c - 0x41)
                s.attrs_reset()
                i += 3
            else:
                i += 1
            continue
        if x == 0x1B:                         # ESC
            e = nxt()
            if e is None:
                break
            if 0x40 <= e <= 0x47:
                s.cfg = e - 0x40
            elif 0x50 <= e <= 0x57:
                s.cbg = e - 0x50
            elif e in (0x4C, 0x4D, 0x4E, 0x4F):  # normal/double height/width/size
                s.size = chr(e)
            elif e == 0x48:                   # flash on
                s.flash = True
            elif e == 0x49:                   # flash off
                s.flash = False
            elif e in (0x4A, 0x4B):           # latent state on/off: TV-mode only, ignored
                pass
            elif e == 0x58:                   # conceal on
                s.conceal = True
            elif e == 0x5F:                   # conceal off
                s.conceal = False
            elif e == 0x59:                   # lining off (underline / contiguous mosaic)
                s.sep = False
            elif e == 0x5A:                   # lining on (underline / separated mosaic)
                s.sep = True
            elif e == 0x5C:                   # positive
                s.inv = False
            elif e == 0x5D:                   # negative
                s.inv = True
            elif e == 0x5E:                   # transparent background
                s.cbg = TRANSPARENT
            i += 2
            continue
        if x == 0x12:                         # REP
            k = nxt()
            if k is not None and 0x40 <= k < 0x80:
                for _ in range(k - 0x40):
                    s.put(s.last, s.last_m)
            i += 2
            continue
        if x == 0x19:                         # SS2 (G2)
            g = nxt()
            if g in G2_ACCENT:
                base = nxt(2)
                if base is not None and 0x20 <= base < 0x7F:
                    import unicodedata
                    s.put(unicodedata.normalize("NFC", chr(base) + G2_ACCENT[g]))
                    i += 3
                    continue
                i += 2
                continue
            if g is not None:
                s.put(G2_CHARS.get(g, "?"))
            i += 2
            continue
        if x == 0x0E:
            s.g1 = True
        elif x == 0x0F:
            s.g1 = False
        elif x == 0x0D:
            s.c = 0
        elif x == 0x0A:
            s.r += 1
        elif x == 0x08:
            s.c = max(0, s.c - 1)
        elif x == 0x09:
            s.c += 1
        elif x == 0x0B:
            s.r = max(0, s.r - 1)
        elif x == 0x0C:                       # FF: clear screen, home
            s.__init__(s.level2)
        elif x == 0x18:                       # CAN: spaces with current attributes to end of row
            r0, c0 = s.r, s.c
            while s.r == r0 and s.c < COLS:
                s.put(" ")
            s.r, s.c = r0, c0
        elif x == 0x1E:                       # RS: home, highest-level defaults
            s.r, s.c = 1, 0
            s.set_level(LEVEL3)
            s.g1 = False
        elif x in (0x03, 0x04):
            pass
        elif x >= 0x20:
            if s.g1:
                # columns 4/5 act like 6/7 in G1 (screenshots: 5/F is a full block)
                if 0x40 <= x < 0x60:
                    x |= 0x20
                s.put(mosaic(x), (x & 0x1F) | ((x & 0x40) >> 1))
            else:
                s.put(G0_SPECIAL.get(x, chr(x)))
        i += 1
    # Double height grows upwards: the row above takes over the background of
    # a row that contains double-height characters (matches the screenshots).
    for r in range(2, ROWS):
        if any(z in "MO" and ch != " " for z, ch in zip(s.sz[r], s.ch[r])):
            for c in range(COLS):
                s.bg[r - 1][c] = s.bg[r][c]
                if s.sz[r][c] in "MO":
                    s.ch[r - 1][c] = " "
                    s.mos[r - 1][c] = None
    return s


def screen_text(s):
    return "\n".join("".join(r).rstrip() for r in s.ch)


def screen_html(s, title):
    """25x40 grid; G1 mosaics drawn as 2x3 blocks, double height grows upwards."""
    cells = []
    for r in range(ROWS):
        for c in range(COLS):
            m = s.mos[r][c]
            # transparent background: no TV picture behind it here, draw black
            fgi, bgi = (0 if v == TRANSPARENT else v for v in (s.fg[r][c], s.bg[r][c]))
            fg, bg = COLORS[fgi], COLORS[bgi]
            sz = s.sz[r][c]
            cls = {"M": " dh", "N": " dw", "O": " ds"}.get(sz, "")
            if s.hid[r][c]:
                cls += " hid"
            if m is not None and s.sepm[r][c]:
                cls += " sep"
            if m is None and s.sepm[r][c]:
                cls += " ul"
            if s.fl[r][c]:
                cls += " fl"
            if m is not None:
                on = fg if not s.hid[r][c] else bg
                inner = "".join('<i style="background:%s"></i>' % (on if m >> k & 1 else bg) for k in range(6))
                cells.append('<b class="m%s" style="background:%s">%s</b>' % (cls, bg, inner))
            else:
                cells.append('<b class="t%s" style="color:%s;background:%s">%s</b>' % (
                    cls, fg, bg, html.escape(s.ch[r][c])))
    return '<div class="page"><h3>%s</h3><div class="scr">%s</div></div>' % (html.escape(title), "".join(cells))


# Benecia (FontStruction by Monique Sammi Brie, CC0) is a 9 x 10 pixel font in
# the style of the Antiope terminals.  At 18 px each font pixel is 2 x 2 screen
# pixels; glyphs it lacks (arrows, frame bars) fall back to a monospace font.
FONT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts", "benecia")
FONT_FILE = os.path.join(FONT_DIR, "benecia.otf.woff2")
FONT_CREDIT = ("Font: Benecia (https://fontstruct.com/fontstructions/show/2950541) "
               "by Monique Sammi Brie, CC0 1.0 Public Domain Dedication")

CSS = """body{background:#222;color:#ddd;font-family:sans-serif}
.page{display:inline-block;margin:10px;vertical-align:top}
h3{font-size:12px;font-weight:normal;max-width:720px;margin:0 0 4px}
.scr{display:grid;grid-template-columns:repeat(40,18px);grid-auto-rows:20px;background:#000;padding:8px;
 font-family:'Benecia','DejaVu Sans Mono','Courier New',monospace;font-size:18px;line-height:20px;
 -webkit-font-smoothing:none;font-smooth:never;font-kerning:none}
.scr b{font-weight:normal;display:block;overflow:visible;text-align:left;white-space:pre}
.scr b.m{display:grid;grid-template-columns:1fr 1fr;grid-template-rows:1fr 1fr 1fr}
.scr b.dh,.scr b.ds{transform-origin:bottom;transform:scaleY(2)}
.scr b.dw{transform-origin:left;transform:scaleX(2)}
.scr b.ds{transform-origin:bottom left;transform:scale(2)}
.scr b.m{overflow:hidden}
.scr b.m i{display:block;min-width:0;min-height:0}
.scr b.sep i{margin:1px}
.scr b.ul{text-decoration:underline}
.scr b.fl{animation:fl 2s step-end infinite}
@keyframes fl{50%{color:transparent}}
.scr b.hid{color:transparent!important}
body.reveal .scr b.hid{color:#888!important}
footer{font-size:11px;color:#999;margin:10px}
"""


def font_css(path=FONT_FILE):
    """@font-face with the font embedded as data URI (output stays one file)."""
    if not path or not os.path.exists(path):
        return ""
    import base64
    data = base64.b64encode(open(path, "rb").read()).decode()
    fmt = "woff2" if path.endswith(".woff2") else ("opentype" if path.endswith(".otf") else "truetype")
    mime = {"woff2": "font/woff2", "opentype": "font/otf", "truetype": "font/ttf"}[fmt]
    return "@font-face{font-family:'Benecia';src:url(data:%s;base64,%s) format('%s')}\n" % (mime, data, fmt)


def content_start(body):
    """Index of first US (row/col) or first printable after the record header."""
    for i in range(2, min(len(body), 40)):
        if body[i] is not None and body[i] & 0x7F == 0x1F:
            return i
    return 2


def page_number(hdr):
    if len(hdr) >= 3 and None not in hdr[:3]:
        return "%X%X%X" % tuple(hdr[:3])
    return "???"


def tc(frame, fps=25.0):
    t = int(frame / fps)
    return "%d:%02d" % (t // 60, t % 60)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("packets")
    ap.add_argument("outdir")
    ap.add_argument("--channel", default=None, help="e.g. 0,7,11 (default: all)")
    ap.add_argument("--font", default=FONT_FILE,
                    help="font for index.html (default: bundled Benecia woff2; empty string = system monospace)")
    ap.add_argument("--min-receptions", type=int, default=2,
                    help="further versions of a page need this many receptions")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    streams = channel_streams(read_packets(a.packets))
    summary = []
    allhtml = []
    for ch, st in sorted(streams.items(), key=lambda kv: -len(kv[1])):
        if a.channel and ch != tuple(int(v) for v in a.channel.split(",")):
            continue
        chname = "%X%X%X" % ch
        with open(os.path.join(a.outdir, "ch%s_stream.bin" % chname), "wb") as fo:
            for seg in flatten(st):
                fo.write(bytes(0 if v is None else v for v in seg[0]))
        bynum = collections.defaultdict(list)
        nrec = 0
        for seg in flatten(st):
            for body, fi in split_pages(seg):
                bynum[page_number(page_header(body))].append((body, fi))
                nrec += 1
        summary.append("Channel %s: %d packets, %d page receptions, %d page numbers" % (
            chname, len(st), nrec, len(bynum)))
        for num in sorted(bynum):
            for vi, c in enumerate(cluster(bynum[num])):
                bodies = [b for b, _ in c["recs"]]
                cnt = len(bodies)
                if cnt < a.min_receptions:
                    continue
                res, conf = vote(bodies)
                # trim at voted ETX EOT
                for i in range(len(res) - 1):
                    if res[i] is not None and res[i + 1] is not None and \
                       res[i] & 0x7F == 0x03 and res[i + 1] & 0x7F == 0x04:
                        res = res[:i + 2]
                        break
                lost = sum(v is None for v in res)
                bad = sum(1 for v in res if v is not None and not odd(v))
                hdr = page_header(res)
                h = parse_header(hdr)
                start = (2 if res and res[0] is not None and res[0] & 0x7F == 0x01 else 1) + h["n"]
                if None in hdr[:h["n"]] or start >= len(res) or res[start] is None or res[start] & 0x7F != 0x1F:
                    start = content_start(res)
                scr = render(res, start, header_defaults(h))
                f0, f1 = c["recs"][0][1], c["recs"][-1][1]
                name = "ch%s_p%s_v%d" % (chname, num, vi)
                title = ("Channel %s page %s version %d: %d receptions (frames %d-%d, %s-%s), %d bytes, "
                         "%d lost, %d parity errors") % (
                    chname, num, vi, cnt, f0, f1, tc(f0), tc(f1), len(res), lost, bad)
                with open(os.path.join(a.outdir, name + ".txt"), "w") as fo:
                    fo.write(title + "\nHeader: %s\n\n" % " ".join("-" if h is None else "%X" % h for h in hdr))
                    fo.write(screen_text(scr) + "\n")
                with open(os.path.join(a.outdir, name + ".bin"), "wb") as fo:
                    fo.write(bytes(0 if v is None else v for v in res))
                allhtml.append(screen_html(scr, title))
                summary.append("  %s" % title)
    with open(os.path.join(a.outdir, "index.html"), "w") as fo:
        font = font_css(a.font) if a.font else ""
        fo.write("<!doctype html><meta charset=utf-8><title>Antiope pages</title><style>%s%s</style>\n"
                 % (font, CSS)
                 + '<p><label><input type=checkbox onchange="document.body.classList.toggle(\'reveal\',this.checked)">'
                 + ' reveal concealed text</label></p>\n'
                 + "\n".join(allhtml)
                 + ("\n<footer>%s</footer>\n" % html.escape(FONT_CREDIT) if font else ""))
    with open(os.path.join(a.outdir, "summary.txt"), "w") as fo:
        fo.write("\n".join(summary) + "\n")
    print("\n".join(summary[:60]))


if __name__ == "__main__":
    main()
