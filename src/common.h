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

/* P1: mu_a=5, mu_b=5, Qa=11, Qb=11, N=14461 */
#define INSTANCE_P1 {5, 5, 11, 11, 14461}

/* P2: mu_a=5, mu_b=6, Qa=11, Qb=13, N=20449 */
#define INSTANCE_P2 {5, 6, 11, 13, 20449}

/* P3: mu_a=6, mu_b=6, Qa=13, Qb=13, N=28561 */
#define INSTANCE_P3 {6, 6, 13, 13, 28561}

/* P4: mu_a=7, mu_b=7, Qa=14, Qb=14, N=38416 */
#define INSTANCE_P4 {7, 7, 14, 14, 38416}

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

/* Transition F(qa,qb,I,d_a,d_b,u) -> next state index */
int transition(int qa, int qb, const State *s, int da, int db, int u, int Qa, int Qb);

/* Nre = (sum Ia + 1) * (sum Ib + 1) */
int nre_for_state(const State *s);

#endif /* COMMON_H */
