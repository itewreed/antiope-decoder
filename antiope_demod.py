#!/usr/bin/env python3
"""Didon (Antiope, CCIR 653 system A) demodulator for bttv raw VBI captures.

Input: raw VBI file, 32 lines x 2048 samples per frame at 35.46895 MHz (bt848).
Output: text file, one line per detected Didon packet:
    frame line quality b3 b4 ... b40   (hex bytes, byte 3 = framing code E7)

VHS (LP) reduces the 3.1 MHz fundamental of the 6.203125 Mbit/s signal by
about 18 dB, so bits are recovered with an MLSE (Viterbi) equaliser whose
channel taps are re-estimated decision-directed per chunk of frames.
"""
import argparse
import os
import multiprocessing as mp
import sys

import numpy as np
from scipy.signal import firwin

import mlse

FS = 35.46895e6
SPB_NOM = FS / 6.203125e6        # 5.7179 samples per bit
P0 = 125.0                       # template position of bit 0
T = 330                          # bits sampled per line (2 lead + 16 run-in + 8 FC + 38*8)
LEAD = 2
NB = 38                          # bytes 3..40
RUNIN_LEVEL = 78.5               # template run-in mean above black
M = 6                            # channel memory (taps = M + 1)
# 3.8 MHz low-pass before bit sampling: removes noise and SECAM chroma
# residue above the data band (byte parity 96.8 % -> 98.3 % on this tape)
LPF = firwin(31, 3.8, fs=FS / 1e6)

# Channel taps estimated on this tape (c = black level, g = response to a 1 bit)
C0 = -62.6
G0 = np.r_[13.8, 53.6, 31.2, 13.2, 10.9, 8.1, 3.3, np.zeros(max(0, M - 6))][:M + 1]


def load(path):
    raw = np.memmap(path, dtype=np.uint8, mode="r")
    n = len(raw) // (32 * 2048)
    return raw[:n * 32 * 2048].reshape(n, 32, 2048)


def to_bytes(B, start=2 + 16):
    w = (1 << np.arange(8)).astype(np.int64)
    b = B[:, start:start + NB * 8].reshape(len(B), NB, 8).astype(np.int64)
    return (b * w).sum(-1).astype(np.uint8)


PAR = np.array([bin(i).count("1") & 1 for i in range(256)], np.uint8)


def process_chunk(args):
    path, f0, f1, tmpl = args
    r = load(path)
    ref = tmpl[100:330] - tmpl[100:330].mean()
    rr = np.dot(ref, ref)
    X, L, info, RUN = [], [], [], []
    blk = r[f0:f1].astype(np.float32)
    for fi in range(f1 - f0):
        for li in range(32):
            x = blk[fi, li].astype(float)
            x = np.convolve(x, LPF, "same")
            x -= x[50:112].mean() + 66.0      # black level before run-in -> -66 (template scale)
            cc = np.correlate(x[60:370], ref, "valid")
            k = int(np.argmax(cc))
            if cc[k] < 0.5 * rr:
                continue
            d = 0.0
            if 0 < k < len(cc) - 1:
                den = cc[k - 1] - 2 * cc[k] + cc[k + 1]
                if den != 0:
                    d = (cc[k - 1] - cc[k + 1]) / (2 * den)
            lag = k + d - 40
            a = int(round(135 + lag))
            RUN.append(x[a:a + 77].mean() + 66.0)   # run-in level above black
            X.append(x)
            L.append(lag)
            info.append((f0 + fi, li))
    if not X:
        return []
    # Data amplitude differs between recordings (TF1 tape: ~20 % below the A2
    # tape).  The run-in (alternating bits) sits at mid level, so its median
    # level gives the gain used when the default channel model finds too few
    # framing codes.
    gain = float(np.median(RUN)) / RUNIN_LEVEL
    X = np.array(X)
    L = np.array(L)
    grid = np.arange(2048)
    c, g = C0, G0.copy()

    def run(c, g):
        best_B, best_cost, best_par = None, None, None
        for off in (-0.6, -0.3, 0.0, 0.3, 0.6):
            for spb in (SPB_NOM - 0.012, SPB_NOM - 0.006, SPB_NOM, SPB_NOM + 0.006, SPB_NOM + 0.012):
                pos = P0 + L[:, None] + off + (np.arange(T)[None, :] - LEAD) * spb
                V = np.array([np.interp(p, grid, x) for p, x in zip(pos, X)])
                B, cost = mlse.viterbi_c(V, g, c)
                if best_B is None:
                    best_B, best_cost, best_V = B, cost, V
                else:
                    u = cost < best_cost
                    best_B[u], best_cost[u], best_V[u] = B[u], cost[u], V[u]
        return best_B, best_cost, best_V

    for it in range(2):
        B, cost, V = run(c, g)
        if it == 1:
            break
        by = to_bytes(B)
        ok = by[:, 0] == 0xE7
        if ok.sum() < 0.3 * len(B) and abs(gain - 1) > 0.1:
            # default model does not fit this recording: rescale it to the
            # measured data amplitude and try again
            c2, g2 = -66.0 + (C0 + 66.0) * gain, G0 * gain
            B2, cost2, V2 = run(c2, g2)
            ok2 = to_bytes(B2)[:, 0] == 0xE7
            if ok2.sum() > ok.sum():
                c, g, B, cost, V, ok = c2, g2, B2, cost2, V2, ok2
        if ok.sum() < 50:
            break
        # decision-directed LS refit of channel on good lines
        rows, ys = [], []
        for b, v in zip(B[ok], V[ok]):
            for j in range(M, T - 8):
                rows.append(b[j - np.arange(M + 1)])
                ys.append(v[j])
        A = np.c_[np.ones(len(rows)), np.array(rows, float)]
        sol = np.linalg.lstsq(A, np.array(ys), rcond=None)[0]
        c, g = sol[0], sol[1:]
    by = to_bytes(B)
    out = []
    for (fi, li), b, q in zip(info, by, cost):
        if bin(int(b[0]) ^ 0xE7).count("1") > 1:
            continue
        out.append("%d %d %.0f %s" % (fi, li, q / T, " ".join("%02X" % v for v in b)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vbi")
    ap.add_argument("out")
    ap.add_argument("--template", default=None, help="npy template (default: tmpl.npy next to script)")
    ap.add_argument("--chunk", type=int, default=250)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=None)
    a = ap.parse_args()
    tmpl = np.load(a.template or os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmpl.npy"))
    n = load(a.vbi).shape[0]
    end = min(n, a.end or n)
    jobs = [(a.vbi, f, min(f + a.chunk, end), tmpl) for f in range(a.start, end, a.chunk)]
    os.environ["OMP_NUM_THREADS"] = "1"
    with open(a.out, "w") as fo, mp.Pool(a.jobs) as pool:
        for i, res in enumerate(pool.imap(process_chunk, jobs)):
            fo.write("".join(s + "\n" for s in res))
            fo.flush()
            print("chunk %d/%d: %d packets" % (i + 1, len(jobs), len(res)), file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
