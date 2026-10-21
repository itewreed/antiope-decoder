#!/usr/bin/env python3
"""Measure the colour of every Antiope colour code from screen photos.

Each screenshot is paired with the decoded page (.bin from antiope_pages.py).
The page is rendered, a 40x25 cell grid is fitted onto the photo, and the
pixels of cells with a known foreground/background code are collected.
"""
import argparse
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import antiope_pages as P  # noqa: E402


def load_page(path):
    raw = list(open(path, "rb").read())
    hdr = P.page_header(raw)
    return P.render(raw, P.content_start(raw), P.header_defaults(P.parse_header(hdr)))


def cell_stats(img, x0, y0, cw, ch):
    """Per cell: median colour and 'ink' colour (pixels farthest from the median)."""
    H, W, _ = img.shape
    med = np.zeros((P.ROWS, P.COLS, 3))
    ink = np.zeros((P.ROWS, P.COLS, 3))
    for r in range(P.ROWS):
        for c in range(P.COLS):
            xa, xb = int(x0 + c * cw + 1), int(x0 + (c + 1) * cw - 1)
            ya, yb = int(y0 + r * ch + 2), int(y0 + (r + 1) * ch - 2)
            if xa < 0 or ya < 0 or xb > W or yb > H or xb <= xa or yb <= ya:
                med[r, c] = np.nan
                ink[r, c] = np.nan
                continue
            px = img[ya:yb, xa:xb].reshape(-1, 3)
            m = np.median(px, 0)
            d = np.linalg.norm(px - m, axis=1)
            k = max(3, len(px) // 10)
            med[r, c] = m
            ink[r, c] = px[np.argsort(d)[-k:]].mean(0)
    return med, ink


def cell_means(ii, x0, y0, cw, ch):
    """Mean colour of every cell's interior from an integral image."""
    c = np.arange(P.COLS)
    r = np.arange(P.ROWS)
    xa = np.clip((x0 + c * cw + 2).astype(int), 0, ii.shape[1] - 1)
    xb = np.clip((x0 + (c + 1) * cw - 2).astype(int), 0, ii.shape[1] - 1)
    ya = np.clip((y0 + r * ch + 3).astype(int), 0, ii.shape[0] - 1)
    yb = np.clip((y0 + (r + 1) * ch - 3).astype(int), 0, ii.shape[0] - 1)
    S = ii[yb][:, xb] - ii[ya][:, xb] - ii[yb][:, xa] + ii[ya][:, xa]
    area = np.maximum((yb - ya)[:, None] * (xb - xa)[None, :], 1)
    return S / area[..., None]


def fit_grid(img, scr):
    """Refine grid origin/size so cells with the same bg code look alike."""
    bg = np.array(scr.bg)
    space = np.array([[ch == " " and m is None for ch, m in zip(row, mrow)]
                      for row, mrow in zip(scr.ch, scr.mos)])
    space[0] = False          # status row is drawn outside the page area
    ii = np.zeros((img.shape[0] + 1, img.shape[1] + 1, 3))
    ii[1:, 1:] = img.cumsum(0).cumsum(1)
    best = None
    # Default geometry of these 768x576 captures: rows 1-24 fill y 60..538,
    # columns fill x 54..700 (row 0 is the status row above the page area).
    for cw in (16.05, 16.15, 16.25):
        for chh in (19.82, 19.92, 20.02):
            for x0 in np.arange(51, 57.5, 1.0):
                for y0 in np.arange(37, 43.5, 1.0):
                    if x0 + 40 * cw > img.shape[1] or y0 + 25 * chh > img.shape[0]:
                        continue
                    med = cell_means(ii, x0, y0, cw, chh)
                    err = 0.0
                    for code in np.unique(bg[space]):
                        v = med[space & (bg == code)]
                        err += ((v - v.mean(0)) ** 2).sum()
                    err /= space.sum()
                    if best is None or err < best[0]:
                        best = (err, x0, y0, cw, chh)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pages_dir", help="output of antiope_pages.py")
    ap.add_argument("pairs", nargs="+", help="screenshot.png=PAGE (e.g. shot17.png=062); PAGE may name a file like 300_v9")
    a = ap.parse_args()
    bgs = {k: [] for k in range(8)}
    fgs = {k: [] for k in range(8)}
    for pr in a.pairs:
        shot, page = pr.split("=")
        scr = load_page(os.path.join(a.pages_dir, "ch07B_p%s.bin" % (page if "_v" in page else page + "_v0")))
        img = np.asarray(Image.open(shot).convert("RGB")).astype(float)
        err, x0, y0, cw, chh = fit_grid(img, scr)
        med, ink = cell_stats(img, x0, y0, cw, chh)
        nb = nf = 0
        for r in range(1, P.ROWS):
            for c in range(P.COLS):
                if np.isnan(med[r, c, 0]) or scr.mos[r][c] is not None:
                    continue
                if scr.ch[r][c] == " ":
                    bgs[scr.bg[r][c]].append(med[r, c])
                    nb += 1
                elif scr.ch[r][c].isalnum() and scr.fg[r][c] != scr.bg[r][c]:
                    fgs[scr.fg[r][c]].append((ink[r, c], scr.bg[r][c]))
                    nf += 1
        print("%s page %s: grid x0=%d y0=%d cell=%.1fx%.1f rms=%.1f, %d bg / %d text cells" % (
            os.path.basename(shot), page, x0, y0, cw, chh, err ** 0.5, nb, nf))
    names = ["black", "red", "green", "yellow", "blue", "magenta", "cyan", "white"]
    print("\ncode  name      background (median RGB, n)        text (mean ink RGB, n)")
    for k in range(8):
        b = np.median(np.array(bgs[k]), 0) if bgs[k] else None
        f = np.mean(np.array([v for v, _ in fgs[k]]), 0) if fgs[k] else None
        print("%d     %-8s  %-34s %s" % (
            k, names[k],
            "-" if b is None else "%3d %3d %3d  (%d)" % (*b, len(bgs[k])),
            "-" if f is None else "%3d %3d %3d  (%d, on bg %s)" % (
                *f, len(fgs[k]), sorted(set(names[x] for _, x in fgs[k])))))


if __name__ == "__main__":
    main()
