/**
 * common.c - Shared utilities for perishable inventory VI
 * Paper: Ortega et al. (2018), Sect 2-3
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

/* State encoding: (Ia1,Ia2,Ib1,Ib2) lexicographic, Ia1 varies fastest.
 * Table 3: Qa=11 means max inventory 10 (0-indexed). N = Qa^M * Qb^M.
 * Ia[r], Ib[r] in {0..Qa-1}, {0..Qb-1}. j = Ia1 + Ia2*Qa + Ib1*Qa^M + Ib2*Qa^M*Qb. */
int state_to_index(const State *s, int Qa, int Qb) {
    int j = 0;
    int stride = 1;
    for (int r = 0; r < M; r++) {
        j += s->Ia[r] * stride;
        stride *= Qa;
    }
    for (int r = 0; r < M; r++) {
        j += s->Ib[r] * stride;
        stride *= Qb;
    }
    return j;
}

void index_to_state(int j, int Qa, int Qb, State *s) {
    memset(s, 0, sizeof(State));
    int stride_ab = 1;
    for (int r = 0; r < M; r++) stride_ab *= Qa;
    for (int r = 0; r < M; r++) stride_ab *= Qb;
    for (int r = M - 1; r >= 0; r--) {
        stride_ab /= Qb;
        s->Ib[r] = j / stride_ab;
        j %= stride_ab;
    }
    for (int r = M - 1; r >= 0; r--) {
        stride_ab /= Qa;
        s->Ia[r] = j / stride_ab;
        j %= stride_ab;
    }
}

/* N = Qa^M * Qb^M (Table 3: Qa=11 -> max inv 10, 11 values per dim) */
int state_space_size(int Qa, int Qb) {
    int n = 1;
    for (int r = 0; r < M; r++) n *= Qa;
    for (int r = 0; r < M; r++) n *= Qb;
    return n;
}

/* Paper Algorithm 3 line 9: Nre = (sum_r I_ar + 1) * (sum_r I_br + 1) */
int nre_for_state(const State *s) {
    int sum_a = 0, sum_b = 0;
    for (int r = 0; r < M; r++) {
        sum_a += s->Ia[r];
        sum_b += s->Ib[r];
    }
    return (sum_a + 1) * (sum_b + 1);
}

/* pu(u,y) - Eq (12): P(substitution demand = u | db > y)
 * pu(u,y) = sum_{x=u}^{inf} P(db=x+y) * Bin(u,x,gamma). Extended for u=0. */
static double pu_func(int u, int y, double mu_b, double gamma) {
    if (u < 0) return 0.0;
    double sum = 0.0;
    int x_start = (u > 0) ? u : 0;
    for (int x = x_start; x < x_start + 200; x++) {
        double p_db = poisson_pmf(x + y, mu_b);
        if (p_db < 1e-15 && x > y + 5) break;
        sum += p_db * binomial_pmf(u, x, gamma);
    }
    return sum;
}

/* pz(z,y) - Eq (13): P(total demand for a = z | db > y)
 * pz(z,y) = sum_{x=0}^{z} P(da=x) * pu(z-x, y) */
static double pz_func(int z, int y, double mu_a, double mu_b, double gamma) {
    double sum = 0.0;
    for (int x = 0; x <= z; x++) {
        sum += poisson_pmf(x, mu_a) * pu_func(z - x, y, mu_b, gamma);
    }
    return sum;
}

/* Expected sales for two-product state j. Paper: Eq (9) for single product.
 * Two products with substitution (Sect 3.1):
 * - When db <= y: no subst, demand_a = da, sales_a = min(da,Ya), sales_b = min(db,Yb)
 * - When db > y: subst, total demand for a = da+u = z with pz(z,y), sales_a = min(z,Ya), sales_b = Yb
 * Esale = E[sales_a] + E[sales_b], with sa=sb=1 (revenue = quantity) */
double expected_sales(int j, int Qa, int Qb, int mu_a, int mu_b, double gamma) {
    State s;
    index_to_state(j, Qa, Qb, &s);
    int Ya = 0, Yb = 0;
    for (int r = 0; r < M; r++) {
        Ya += s.Ia[r];
        Yb += s.Ib[r];
    }

    double Esale_a = 0.0, Esale_b = 0.0;

    /* Product b: Esale_b = sum_{db=1}^{Yb} db * P(db) + Yb * P(db > Yb) = sum_{db=0}^{inf} min(db,Yb)*P(db) */
    for (int db = 0; db < Yb + 100; db++) {
        double p = poisson_pmf(db, mu_b);
        if (p < 1e-15 && db > Yb) break;
        Esale_b += (db < Yb ? db : Yb) * p;
    }

    /* Product a */
    if (Ya == 0) {
        return Esale_b;  /* no stock of a */
    }

    /* P(db <= Yb) */
    double P_db_le_y = 0.0;
    for (int db = 0; db <= Yb; db++)
        P_db_le_y += poisson_pmf(db, mu_b);

    /* When db <= Yb: no substitution. Esale_a = sum_{da} min(da,Ya)*P(da) */
    double Esale_a_no_subst = 0.0;
    for (int da = 0; da < Ya + 100; da++) {
        double p = poisson_pmf(da, mu_a);
        if (p < 1e-15 && da > Ya) break;
        Esale_a_no_subst += (da < Ya ? da : Ya) * p;
    }

    /* When db > Yb: substitution. Esale_a = sum_z min(z,Ya)*pz(z,Yb) */
    double Esale_a_subst = 0.0;
    double P_db_gt_y = 1.0 - P_db_le_y;
    if (P_db_gt_y > 1e-15) {
        for (int z = 0; z <= Ya + 100; z++) {
            double pz = pz_func(z, Yb, mu_a, mu_b, gamma);
            if (pz < 1e-15 && z > Ya) break;
            Esale_a_subst += (z < Ya ? z : Ya) * pz;
        }
    }

    Esale_a = P_db_le_y * Esale_a_no_subst + P_db_gt_y * Esale_a_subst;
    return Esale_a + Esale_b;
}

/* Transition F(qa,qb,I,da,db,u). Paper Eqs (4),(5) FIFO.
 * I_M,t = Q; I_r,t = (I_{r+1} - (d - sum_{m=1}^r I_m)^+)^+ for r=1..M-1.
 * Two products: d_a_eff = min(da+u, Ya), d_b_eff = min(db, Yb). */
int transition(int qa, int qb, const State *s, int da, int db, int u, int Qa, int Qb) {
    int Ya = 0, Yb = 0;
    for (int r = 0; r < M; r++) {
        Ya += s->Ia[r];
        Yb += s->Ib[r];
    }
    int d_a_eff = da + u;
    if (d_a_eff > Ya) d_a_eff = Ya;
    int d_b_eff = (db > Yb) ? Yb : db;

    State next;
    next.Ia[M - 1] = qa;
    next.Ib[M - 1] = qb;

    /* Product a FIFO: I_a1_new = (I_a2 - (d_a_eff - I_a1)^+)^+ */
    int cum_a = s->Ia[0];
    int prev_a = s->Ia[1];
    int shortfall_a = (d_a_eff > cum_a) ? (d_a_eff - cum_a) : 0;
    next.Ia[0] = (prev_a > shortfall_a) ? (prev_a - shortfall_a) : 0;

    int cum_b = s->Ib[0];
    int prev_b = s->Ib[1];
    int shortfall_b = (d_b_eff > cum_b) ? (d_b_eff - cum_b) : 0;
    next.Ib[0] = (prev_b > shortfall_b) ? (prev_b - shortfall_b) : 0;

    return state_to_index(&next, Qa, Qb);
}

double substitution_future_value(int qa, int qb, const State *s, int Qa, int Qb,
    double mu_a, double mu_b, double gamma, const double *W) {
    int Yb = 0;
    for (int r = 0; r < M; r++) Yb += s->Ib[r];
    double P_db_gt = 0.0;
    for (int db = Yb + 1; db < Yb + 100; db++) {
        double p = poisson_pmf(db, mu_b);
        if (p < 1e-15 && db > Yb + 10) break;
        P_db_gt += p;
    }
    int Ya = 0;
    for (int r = 0; r < M; r++) Ya += s->Ia[r];
    double sum = 0.0;
    for (int z = 0; z <= Ya + 50; z++) {
        double pz = pz_func(z, Yb, mu_a, mu_b, gamma);
        if (pz < 1e-15 && z > Ya) break;
        int d_a = (z < Ya) ? z : Ya;
        int k = transition_eff(qa, qb, s, d_a, Yb, Qa, Qb);
        sum += pz * W[k];
    }
    return P_db_gt * sum;
}

/* Same as transition but with effective demands directly (avoids da,db,u) */
int transition_eff(int qa, int qb, const State *s, int d_a_eff, int d_b_eff, int Qa, int Qb) {
    int Ya = 0, Yb = 0;
    for (int r = 0; r < M; r++) {
        Ya += s->Ia[r];
        Yb += s->Ib[r];
    }
    if (d_a_eff > Ya) d_a_eff = Ya;
    if (d_b_eff > Yb) d_b_eff = Yb;

    State next;
    next.Ia[M - 1] = qa;
    next.Ib[M - 1] = qb;
    int cum_a = s->Ia[0], prev_a = s->Ia[1];
    int shortfall_a = (d_a_eff > cum_a) ? (d_a_eff - cum_a) : 0;
    next.Ia[0] = (prev_a > shortfall_a) ? (prev_a - shortfall_a) : 0;
    int cum_b = s->Ib[0], prev_b = s->Ib[1];
    int shortfall_b = (d_b_eff > cum_b) ? (d_b_eff - cum_b) : 0;
    next.Ib[0] = (prev_b > shortfall_b) ? (prev_b - shortfall_b) : 0;
    return state_to_index(&next, Qa, Qb);
}
