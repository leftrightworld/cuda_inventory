/**
 * common.c - Shared utilities for perishable inventory VI
 */

#include "common.h"
#include <math.h>
#include <string.h>

static double factorial(int n) {
    if (n <= 1) return 1.0;
    double f = 1.0;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

double poisson_pmf(int k, double mu) {
    if (k < 0) return 0.0;
    return exp(-mu) * pow(mu, k) / factorial(k);
}

double binomial_pmf(int u, int x, double gamma) {
    if (u < 0 || u > x) return 0.0;
    return factorial(x) / (factorial(u) * factorial(x - u)) * pow(gamma, u) * pow(1 - gamma, x - u);
}

/* State j: lexicographic order (Ia1,Ia2,...,Ib1,Ib2,...) */
int state_to_index(const State *s, int Qa, int Qb) {
    int j = 0;
    int stride = 1;
    for (int r = M - 1; r >= 0; r--) {
        j += s->Ia[r] * stride;
        stride *= (Qa + 1);
    }
    for (int r = M - 1; r >= 0; r--) {
        j += s->Ib[r] * stride;
        stride *= (Qb + 1);
    }
    return j;
}

void index_to_state(int j, int Qa, int Qb, State *s) {
    memset(s, 0, sizeof(State));
    int stride_b = 1;
    for (int r = 0; r < M; r++) stride_b *= (Qb + 1);
    int stride_a = stride_b * (Qa + 1);
    for (int r = M - 1; r >= 0; r--) {
        stride_a /= (Qa + 1);
        s->Ia[r] = j / stride_a;
        j %= stride_a;
    }
    for (int r = M - 1; r >= 0; r--) {
        stride_b /= (Qb + 1);
        s->Ib[r] = j / stride_b;
        j %= stride_b;
    }
}

int nre_for_state(const State *s) {
    int sum_a = 0, sum_b = 0;
    for (int r = 0; r < M; r++) {
        sum_a += s->Ia[r];
        sum_b += s->Ib[r];
    }
    return (sum_a + 1) * (sum_b + 1);
}

/* Stubs - full implementation in Phase 2 */
double expected_sales(int j, int Qa, int Qb, int mu_a, int mu_b, double gamma) {
    (void)j; (void)Qa; (void)Qb; (void)mu_a; (void)mu_b; (void)gamma;
    return 0.0;
}

int transition(int qa, int qb, const State *s, int da, int db, int u, int Qa, int Qb) {
    (void)qa; (void)qb; (void)s; (void)da; (void)db; (void)u; (void)Qa; (void)Qb;
    return 0;
}
