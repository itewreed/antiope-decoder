"""Bit-spaced MLSE (Viterbi) for heavily low-passed NRZ data, vectorised over lines."""
import numpy as np


def viterbi(v, g, c, pre_bits=None):
    """v: (N, T) bit-centre samples; model v[:, j-D] = c + sum_m g[m] * b[j-m].

    g has length M+1 (causal taps, g[0] = newest bit). D = lead (anticausal span)
    is handled by the caller by shifting v. Returns (N, T) bit decisions.
    pre_bits: initial state assumed all zeros (black level before run-in).
    """
    N, T = v.shape
    M = len(g) - 1
    S = 1 << M                       # state = last M bits, bit0 = newest
    states = np.arange(S)
    # expected output for (state, new bit)
    hist = ((states[:, None] >> np.arange(M)[None, :]) & 1).astype(float)  # (S, M) b[j-1..j-M]
    exp0 = c + hist @ g[1:]
    exp1 = exp0 + g[0]
    expv = np.stack([exp0, exp1], 1)            # (S, 2)
    nxt = ((states[:, None] << 1) | np.array([0, 1])[None, :]) & (S - 1)  # (S, 2)
    cost = np.full((N, S), 1e18)
    cost[:, 0] = 0.0
    back = np.zeros((T, N, S), np.int32)
    # predecessor table: for state s', predecessors p with nxt[p, s'&1] == s'
    # s' = ((p << 1) | b) & mask  -> b = s' & 1, p in {(s' >> 1), (s' >> 1) | (1 << (M-1))}
    p0 = states >> 1
    p1 = p0 | (1 << (M - 1))
    bnew = states & 1
    e0 = expv[p0, bnew]                          # (S,)
    e1 = expv[p1, bnew]
    for j in range(T):
        x = v[:, j][:, None]
        c0 = cost[:, p0] + (x - e0[None, :]) ** 2
        c1 = cost[:, p1] + (x - e1[None, :]) ** 2
        pick = c1 < c0
        cost = np.where(pick, c1, c0)
        back[j] = np.where(pick, p1[None, :], p0[None, :])
    s = np.argmin(cost, 1)
    final = cost[np.arange(N), s]
    out = np.zeros((N, T), np.uint8)
    rows = np.arange(N)
    for j in range(T - 1, -1, -1):
        out[:, j] = s & 1
        s = back[j, rows, s]
    return out, final


import ctypes
import os

_lib = None


def viterbi_c(v, g, c):
    """Same as viterbi() but in C (libviterbi.so next to this file)."""
    global _lib
    if _lib is None:
        _lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libviterbi.so"))
    v = np.ascontiguousarray(v, np.float64)
    g = np.ascontiguousarray(g, np.float64)
    N, T = v.shape
    out = np.zeros((N, T), np.uint8)
    cost = np.zeros(N)
    P = ctypes.POINTER
    _lib.viterbi(v.ctypes.data_as(P(ctypes.c_double)), N, T, g.ctypes.data_as(P(ctypes.c_double)),
                 len(g) - 1, ctypes.c_double(c), out.ctypes.data_as(P(ctypes.c_uint8)),
                 cost.ctypes.data_as(P(ctypes.c_double)))
    return out, cost
