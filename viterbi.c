/* Bit-spaced MLSE for Didon lines.  Build: cc -O3 -shared -fPIC -o libviterbi.so viterbi.c */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* v: N*T samples, g: M+1 taps (g[0] = newest bit), c: offset.
   out: N*T bits, cost: N final path metrics.  Start state = all zeros. */
void viterbi(const double *v, int N, int T, const double *g, int M, double c,
             uint8_t *out, double *cost)
{
    int S = 1 << M;
    double *e = malloc(sizeof(double) * S * 2);
    double *pc = malloc(sizeof(double) * S), *nc = malloc(sizeof(double) * S);
    uint8_t *bk = malloc((size_t)T * S);
    for (int s = 0; s < S; s++) {
        double x = c;
        for (int m = 0; m < M; m++)
            if ((s >> m) & 1) x += g[m + 1];
        e[2 * s] = x;
        e[2 * s + 1] = x + g[0];
    }
    for (int n = 0; n < N; n++) {
        const double *vn = v + (size_t)n * T;
        for (int s = 0; s < S; s++) pc[s] = 1e30;
        pc[0] = 0;
        for (int j = 0; j < T; j++) {
            double x = vn[j];
            uint8_t *b = bk + (size_t)j * S;
            for (int s = 0; s < S; s++) {
                int p0 = s >> 1, p1 = p0 | (1 << (M - 1)), bit = s & 1;
                double d0 = x - e[2 * p0 + bit], d1 = x - e[2 * p1 + bit];
                double c0 = pc[p0] + d0 * d0, c1 = pc[p1] + d1 * d1;
                if (c1 < c0) { nc[s] = c1; b[s] = 1; } else { nc[s] = c0; b[s] = 0; }
            }
            double *t = pc; pc = nc; nc = t;
        }
        int s = 0;
        for (int k = 1; k < S; k++) if (pc[k] < pc[s]) s = k;
        cost[n] = pc[s];
        uint8_t *o = out + (size_t)n * T;
        for (int j = T - 1; j >= 0; j--) {
            o[j] = s & 1;
            int p0 = s >> 1;
            s = bk[(size_t)j * S + s] ? (p0 | (1 << (M - 1))) : p0;
        }
    }
    free(e); free(pc); free(nc); free(bk);
}
