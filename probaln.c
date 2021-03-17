/* The MIT License

   Copyright (C) 2003-2006, 2008-2010 by Heng Li <lh3lh3@live.co.uk>
   Copyright (C) 2016-2017, 2020 Genome Research Ltd.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <errno.h>
#include "htslib/hts.h"

/*****************************************
 * Probabilistic banded glocal alignment *
 *****************************************/

#define EI .25
#define EM .33333333333

static float g_qual2prob[256];

#define set_u(u, b, i, k) { int x=(i)-(b); x=x>0?x:0; (u)=((k)-x+1)*3; }

/*
  The topology of the profile HMM:

           /\             /\        /\             /\
           I[1]           I[k-1]    I[k]           I[L]
            ^   \      \    ^    \   ^   \      \   ^
            |    \      \   |     \  |    \      \  |
    M[0]   M[1] -> ... -> M[k-1] -> M[k] -> ... -> M[L]   M[L+1]
                \      \/        \/      \/      /
                 \     /\        /\      /\     /
                       -> D[k-1] -> D[k] ->

   M[0] points to every {M,I}[k] and every {M,I}[k] points to M[L+1].

   On input, ref is the reference sequence and query is the query
   sequence. Both are sequences of 0/1/2/3/4 where 4 stands for an
   ambiguous residue. iqual is the base quality. c sets the gap open
   probability, gap extension probability and band width.

   On output, state and q are arrays of length l_query. The higher 30
   bits give the reference position the query base is matched to and the
   lower two bits can be 0 (an alignment match) or 1 (an
   insertion). q[i] gives the phred scaled posterior probability of
   state[i] being wrong.

   Returns phred-scaled likelihood score, or INT_MIN on failure.
 */
int probaln_glocal(const uint8_t *ref, int l_ref, const uint8_t *query,
                   int l_query, const uint8_t *iqual, const probaln_par_t *c,
                   int *state, uint8_t *q)
{
    double *f = NULL, *b = NULL, *s = NULL, m[9], sI, sM, bI, bM;
    float *qual = NULL;
    int bw, bw2, i, k, is_backward = 1, Pr;

    if ( l_ref<0 || l_query<0 || l_query >= INT_MAX - 2) {
        errno = EINVAL;
        return INT_MIN;
    }
    if (l_ref==0 || l_query==0)
        return 0;  // Is this actually invalid??

    /*** initialization ***/
    is_backward = state && q? 1 : 0;
    bw = l_ref > l_query? l_ref : l_query;
    if (bw > c->bw) bw = c->bw;
    if (bw < abs(l_ref - l_query)) bw = abs(l_ref - l_query);
    bw2 = bw * 2 + 1;
    size_t i_dim = bw2 < l_ref ? (size_t) bw2*3+6 : (size_t) l_ref*3+6;

    // allocate the forward and backward matrices f[][] and b[][]
    // and the scaling array s[]
    // Ideally these callocs would be mallocs + initialisation of
    // the few bits needed.
    if (SIZE_MAX / (l_query+1) / i_dim < sizeof(double)) {
        errno = ENOMEM; // Allocation would fail
        return INT_MIN;
    }
    f = calloc((l_query+1)*i_dim, sizeof(double));
    if (!f) goto fail;
    if (is_backward) {
        b = calloc((l_query+1)*i_dim, sizeof(double));
        if (!b) goto fail;
    }

    // s[] is the scaling factor to avoid underflow
    s = malloc((l_query+2) * sizeof(double));
    if (!s) goto fail;

    // initialize qual
    qual = malloc(l_query * sizeof(float));
    if (!qual) goto fail;
    if (g_qual2prob[0] == 0)
        for (i = 0; i < 256; ++i)
            g_qual2prob[i] = pow(10, -i/10.);
    qual[0] = 0.0; // Should be unused
    for (i = 0; i < l_query; ++i)
        qual[i] = g_qual2prob[iqual? iqual[i] : 30];

    // initialize transition probability
    // the value here seems not to affect results; FIXME: need proof
    sM = sI = 1. / (2 * l_query + 2);
    m[0*3+0] = (1 - c->d - c->d) * (1 - sM);
    m[0*3+1] = m[0*3+2] = c->d * (1 - sM);
    m[1*3+0] = (1 - c->e) * (1 - sI);
    m[1*3+1] = c->e * (1 - sI);
    m[1*3+2] = 0.;
    m[2*3+0] = 1 - c->e;
    m[2*3+1] = 0.;
    m[2*3+2] = c->e;
    bM = (1 - c->d) / l_ref; // (bM+bI)*l_ref==1
    bI = c->d / l_ref;

    /*** forward ***/
    // f[0]
    set_u(k, bw, 0, 0);
    f[0*i_dim+k] = s[0] = 1.;
    { // f[1]
        double *fi = &f[1*i_dim], sum;
        int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1;
        for (k = beg, sum = 0.; k <= end; ++k) {
            int u;
            double e = (ref[k - 1] > 3 || query[0] > 3)? 1. : ref[k - 1] == query[0]? 1. - qual[0] : qual[0] * EM;
            set_u(u, bw, 1, k);
            fi[u+0] = e * bM; fi[u+1] = EI * bI;
            sum += fi[u] + fi[u+1];
        }
        s[1] = sum;
    }
    // f[2..l_query]
    for (i = 2; i <= l_query; ++i) {
        double *fi = &f[i*i_dim], *fi1 = &f[(i-1)*i_dim], sum, qli = qual[i-1];
        int beg = 1, end = l_ref, x;
        uint8_t qyi = query[i - 1];
        x = i - bw; beg = beg > x? beg : x; // band start
        x = i + bw; end = end < x? end : x; // band end

        // NB end-beg is almost always 14 (99.9% of the time)
        // Hence not a large volume to parallelise.
        //
        // Maybe stripe in diagonal doing 14 lines together?
        //
        // Consider rotation? 150x14 vs 14x150 so inner loop
        // takes longer.

        double E[] = {
            qli * EM, // 00
            1. - qli, // 01
            1.,       // 10
            1.,       // 11
        };
        double M = 1./s[i-1];

        // Note this code has the original version listed here (albeit
        // with improved formatting), but we do not compile using
        // -DPROBALN_ORIG.  The purpose of this code is to act as an
        // easier(?) to understand version of the heavily optimised
        // version following it and as an easy validation path in case
        // of any differences in results.
#ifdef PROBALN_ORIG
        for (k = beg, sum = 0.; k <= end; ++k) {
            int u, v11, v01, v10;
            double e;
            e = E[(ref[k - 1] > 3 || qyi > 3)*2 + (ref[k - 1] == qyi)];
            set_u(u, bw, i, k);
            set_u(v11, bw, i-1, k-1);
            set_u(v10, bw, i-1, k);
            set_u(v01, bw, i, k-1);
            fi[u+0] = e * (m[0] * M*fi1[v11+0] + m[3] * M*fi1[v11+1] + m[6] * M*fi1[v11+2]);
            fi[u+1] = EI * (m[1] * M*fi1[v10+0] + m[4] * M*fi1[v10+1]);
            fi[u+2] = m[2] * fi[v01+0] + m[8] * fi[v01+2];
            sum += fi[u] + fi[u+1] + fi[u+2];
        }
#else
        // We use EI*(M*m[1]*? + M*m[4]*?) a lot.  So factor it out here.
        double xm[5];
        xm[0] =    M*m[0];
        xm[1] =    M*m[3];
        xm[2] =    M*m[6];
        xm[3] = EI*M*m[1];
        xm[4] = EI*M*m[4];

        {
            int u, v11;
            set_u(u,   bw, i,   beg);
            set_u(v11, bw, i-1, beg-1);
            // Rather than recompute k->{u,v01,v10,v11} each loop
            // we just increment the pointers.
            double *xi = &fi[u];
            double *yi = &fi1[v11];
            // Derived from xi[0,2] in previous loop iter.
            double  l_x0 = m[2]*xi[0];
            double  l_x2 = m[8]*xi[2];
            for (k = beg, sum = 0.; k <= end; ++k, xi+=3, yi+=3) {
                int cond = (ref[k-1] > 3 || qyi > 3)*2 + (ref[k-1] == qyi);

                double z0 = xm[0]*yi[0];
                double z1 = xm[1]*yi[1];
                double z2 = xm[2]*yi[2];
                double z3 = xm[3]*yi[3];
                double z4 = xm[4]*yi[4];

                xi[0] = E[cond] * (z0+z1+z2);
                xi[1] = z3 + z4;
                xi[2] = l_x0 + l_x2;
                sum  += xi[0] + xi[1] + xi[2];

                l_x0 = m[2]*xi[0];
                l_x2 = m[8]*xi[2];
            }
        }
#endif
        s[i] = sum;
    }

    { // f[l_query+1]
        double sum;
        double M = 1./s[l_query];
        for (k = 1, sum = 0.; k <= l_ref; ++k) {
            int u;
            set_u(u, bw, l_query, k);
            if (u < 3 || u >= i_dim - 3) continue;
            sum += M*f[l_query*i_dim + u+0] * sM + M*f[l_query*i_dim + u+1] * sI;
        }
        s[l_query+1] = sum; // the last scaling factor
    }
    { // compute likelihood
        double p = 1., Pr1 = 0.;
        for (i = 0; i <= l_query + 1; ++i) {
            p *= s[i];
            if (p < 1e-100) Pr1 += -4.343 * log(p), p = 1.;
        }
        Pr1 += -4.343 * log(p * l_ref * l_query);
        Pr = (int)(Pr1 + .499);
        if (!is_backward) { // skip backward and MAP
            free(f); free(s); free(qual);
            return Pr;
        }
    }
    /*** backward ***/
    // b[l_query] (b[l_query+1][0]=1 and thus \tilde{b}[][]=1/s[l_query+1]; this is where s[l_query+1] comes from)
    for (k = 1; k <= l_ref; ++k) {
        int u;
        double *bi = &b[l_query*i_dim];
        set_u(u, bw, l_query, k);
        if (u < 3 || u >= i_dim - 3) continue;
        bi[u+0] = sM / s[l_query] / s[l_query+1]; bi[u+1] = sI / s[l_query] / s[l_query+1];
    }
    // b[l_query-1..1]
    for (i = l_query - 1; i >= 1; --i) {
        int beg = 1, end = l_ref, x;
        double *bi = &b[i*i_dim], *bi1 = &b[(i+1)*i_dim], y = (i > 1), qli1 = qual[i];
        uint8_t qyi1 = query[i];
        x = i - bw; beg = beg > x? beg : x;
        x = i + bw; end = end < x? end : x;
        double E[] = {
            qli1 * EM, //000
            1. - qli1, //001
            1.,        //010
            1.,        //011
            //0,0,0,0    //1xx
        };

#ifdef PROBALN_ORIG
        for (k = end; k >= beg; --k) {
            int u, v11, v01, v10;
            double e;
            set_u(u, bw, i, k);
            set_u(v11, bw, i+1, k+1);
            set_u(v10, bw, i+1, k);
            set_u(v01, bw, i, k+1);
            e = (k>=l_ref)?0 :E[(ref[k] > 3 || qyi1 > 3)*2 + (ref[k] == qyi1)] * bi1[v11];
            bi[u+0] = e * m[0] + EI * m[1] * bi1[v10+1] + m[2] * bi[v01+2]; // bi1[v11] has been foled into e.
            bi[u+1] = e * m[3] + EI * m[4] * bi1[v10+1];
            bi[u+2] = (e * m[6] + m[8] * bi[v01+2]) * y;
//          fprintf(stderr, "B (%d,%d;%d): %lg,%lg,%lg\n", i, k, u, bi[u], bi[u+1], bi[u+2]); // DEBUG
        }
        // rescale
        int _beg, _end;
        set_u(_beg, bw, i, beg); set_u(_end, bw, i, end); _end += 2;
        for (k = _beg, y = 1./s[i]; k <= _end; ++k) bi[k] *= y;
#else
        {
            int u, v10;
            set_u(u,   bw, i,   end);
            set_u(v10, bw, i+1, end);
            // Rather than recompute k->{u,v01,v10,v11} each loop
            // we just increment the pointers.
            double *xi = &bi[u];
            double *yi = &bi1[v10];
            // NB xi[5] is equiv to v01+2.
            double xi_5 = xi[5];
            // Manual loop invariant removal
            double e1 = EI*m[1];
            double e4 = EI*m[4];
            // Do renorm too in the same pass.
            double n = 1./s[i];
            for (k = end; k >= beg; --k, xi -= 3, yi -= 3) {
                double e = (k>=l_ref)
                    ? 0
                    : E[(ref[k]>3 || qyi1>3)*2 + (ref[k] == qyi1)] * yi[3];

                xi[1] =  e * m[3] +   e4 * yi[1];
                xi[0] =  e * m[0] +   e1 * yi[1]  + m[2] * xi_5;
                xi[2] = (e * m[6] + m[8] * xi_5) * y;
                // bi[u+2] from this iter becomes bi[v01+2] in next iter
                xi_5 = xi[2];

                // rescale
                xi[1] *= n;
                xi[0] *= n;
                xi[2] *= n;
            }
        }
#endif
    }
    { // b[0]
        int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1;
        double sum = 0.;
        for (k = end; k >= beg; --k) {
            int u;
            double e = (ref[k - 1] > 3 || query[0] > 3)? 1. : ref[k - 1] == query[0]? 1. - qual[0] : qual[0] * EM;
            set_u(u, bw, 1, k);
            if (u < 3 || u >= i_dim - 3) continue;
            sum += e * b[1*i_dim + u+0] * bM + EI * b[1*i_dim + u+1] * bI;
        }
        set_u(k, bw, 0, 0);
        b[0*i_dim + k] = sum / s[0]; // if everything works as is expected, b[0][k] == 1.0
    }
    /*** MAP ***/
    for (i = 1; i <= l_query; ++i) {
        double sum = 0., *fi = &f[i*i_dim], *bi = &b[i*i_dim], max = 0.;
        int beg = 1, end = l_ref, x, max_k = -1;
        x = i - bw; beg = beg > x? beg : x;
        x = i + bw; end = end < x? end : x;
        double M = 1./s[i];
#ifdef PROBALN_ORIG
        for (k = beg; k <= end; ++k) {
            int u;
            double z;
            set_u(u, bw, i, k);
            z = M*fi[u+0] * bi[u+0];
            if (z > max) max = z, max_k = (k-1)<<2 | 0;
            sum += z;
            z = M*fi[u+1] * bi[u+1];
            if (z > max) max = z, max_k = (k-1)<<2 | 1;
            sum += z;
        }
#else
        {
            int u;
            set_u(u, bw, i, beg);
            for (k = beg; k <= end; ++k, u+=3) {
                double z1, z2;
                z1 = M*fi[u+0] * bi[u+0];
                z2 = M*fi[u+1] * bi[u+1];
                int which = z2 > z1; // strictly z2 >= z1 matches old code
                double zm = which ? z2 : z1;
                if (zm > max) {
                    max   = zm;
                    max_k = (k-1)<<2 | which;
                }
                sum += z1 + z2;
            }
        }
#endif
        max /= sum; sum *= s[i]; // if everything works as is expected, sum == 1.0
        if (state) state[i-1] = max_k;
        if (q) k = (int)(-4.343 * log(1. - max) + .499), q[i-1] = k > 100? 99 : k;
#ifdef PROBALN_MAIN
        k = 0;
        set_u(k, bw, 0, 0);
        fprintf(stderr, "(%.10lg,%.10lg) (%d,%d:%c,%c:%d) %lg\n", b[0][k], sum, i-1, max_k>>2,
                "ACGT"[query[i - 1]], "ACGT"[ref[(max_k>>2)]], max_k&3, max); // DEBUG
#endif
    }

    /*** free ***/
    free(f); free(b); free(s); free(qual);
    return Pr;

 fail:
    free(f); free(b); free(s); free(qual);
    return INT_MIN;
}

#ifdef PROBALN_MAIN
#include <unistd.h>
int main(int argc, char *argv[])
{
    uint8_t conv[256], *iqual, *ref, *query;
    probaln_par_t par = { 0.001, 0.1, 10 };
    int c, l_ref, l_query, i, q = 30, b = 10, P;
    while ((c = getopt(argc, argv, "b:q:")) >= 0) {
        switch (c) {
        case 'b': b = atoi(optarg); break;
        case 'q': q = atoi(optarg); break;
        }
    }
    if (optind + 2 > argc) {
        fprintf(stderr, "Usage: %s [-q %d] [-b %d] <ref> <query>\n", argv[0], q, b); // example: acttc attc
        return 1;
    }
    memset(conv, 4, 256);
    conv['a'] = conv['A'] = 0; conv['c'] = conv['C'] = 1;
    conv['g'] = conv['G'] = 2; conv['t'] = conv['T'] = 3;
    ref = (uint8_t*)argv[optind]; query = (uint8_t*)argv[optind+1];
    l_ref = strlen((char*)ref); l_query = strlen((char*)query);
    for (i = 0; i < l_ref; ++i) ref[i] = conv[ref[i]];
    for (i = 0; i < l_query; ++i) query[i] = conv[query[i]];
    iqual = malloc(l_query);
    memset(iqual, q, l_query);
    par.bw = b;
    P = probaln_glocal(ref, l_ref, query, l_query, iqual, &par, 0, 0);
    fprintf(stderr, "%d\n", P);
    free(iqual);
    return 0;
}
#endif
