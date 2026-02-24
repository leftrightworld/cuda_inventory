/**
 * common.h - Shared definitions for perishable inventory VI
 * Reproduces Ortega et al. (2018) two-product case
 */

#ifndef COMMON_H
#define COMMON_H

#include <stddef.h>
#include <stdint.h>

/* --- Parameters (Table 3) --- */
#define M 2                    /* shelf life (max remaining periods) */
#define SA 1.0
#define SB 1.0
#define CA 0.5
#define CB 0.5
#define GAMMA 0.5              /* substitution probability */
#define EPSILON 1e-4           /* convergence threshold */
#define MAX_ITER 100           /* max iterations for benchmark */

/* Instance parameters */
typedef struct {
    int mu_a, mu_b;            /* Poisson mean demand */
    int Qa, Qb;                /* max order quantities */
    int N;                     /* state space size */
} Instance;

/* Table 3 instances - N = Qa^M * Qb^M (Qa=11 -> max inv 10, 11 vals/dim) */
#define INSTANCE_P1 {5, 5, 11, 11, 0}   /* N=14641 */
#define INSTANCE_P2 {5, 6, 11, 13, 0}   /* N=20449 */
#define INSTANCE_P3 {6, 6, 13, 13, 0}   /* N=28561 */
#define INSTANCE_P4 {7, 7, 14, 14, 0}   /* N=38416 */

/* --- Inventory state --- */
typedef struct {
    int Ia[M];                 /* I_a1, I_a2, ... */
    int Ib[M];                 /* I_b1, I_b2, ... */
} State;

/* --- Functions (to implement) --- */

/* Poisson PMF: P(X=k) for mean mu */
double poisson_pmf(int k, double mu);

/* Binomial PMF: Bin(u, x, gamma) */
double binomial_pmf(int u, int x, double gamma);

/* State index j <-> State (Ia,Ib) */
int state_to_index(const State *s, int Qa, int Qb);
void index_to_state(int j, int Qa, int Qb, State *s);

/* Expected sales for state j */
double expected_sales(int j, int Qa, int Qb, int mu_a, int mu_b, double gamma);

/* Transition F(qa,qb,I,da,db,u) -> next state index */
int transition(int qa, int qb, const State *s, int da, int db, int u, int Qa, int Qb);

/* Transition by effective demands (for aggregated probabilities) */
int transition_eff(int qa, int qb, const State *s, int d_a_eff, int d_b_eff, int Qa, int Qb);

/* Substitution future value: sum over z of pz(z,Yb)*P(db>Yb)*W[k] for k=F(qa,qb,s,z,Yb) */
double substitution_future_value(int qa, int qb, const State *s, int Qa, int Qb,
    double mu_a, double mu_b, double gamma, const double *W);

/* Nre = (sum Ia + 1) * (sum Ib + 1) - Algorithm 3 line 9 */
int nre_for_state(const State *s);

/* Compute N = Qa^M * Qb^M */
int state_space_size(int Qa, int Qb);

#endif /* COMMON_H */
