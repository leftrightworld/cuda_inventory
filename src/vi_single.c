#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>



/* Paper Section 2.3 uses: s=1, c=0.5, mu=5, epsilon=1e-4 */
#define S_PRICE   1.0
#define C_COST    0.5
#define MU_DEMAND  5.0
#define EPS       1e-4
#define MAX_IT    200

/* Simulation settings (paper: 400,000 periods) */
#define SIM_T     400000
#define BURN_IN   2000

/* ---------- utilities ---------- */

static int ipow(int a, int b) {
    int r = 1;
    for (int i = 0; i < b; i++) r *= a;
    return r;
}

/* Encode/decode with I[0] fastest (lowest digit):
   j = I[0] + I[1]*(Q+1) + ... + I[M-1]*(Q+1)^(M-1) */
static void index_to_state(int j, int *I, int Q, int M) {
    for (int r = 0; r < M; r++) {
        I[r] = j % (Q + 1);
        j /= (Q + 1);
    }
}
static int state_to_index(const int *I, int Q, int M) {
    int j = 0;
    int stride = 1;
    for (int r = 0; r < M; r++) {
        j += I[r] * stride;
        stride *= (Q + 1);
    }
    return j;
}
static int state_space(int Q, int M) {
    return ipow(Q + 1, M);
}

/* Poisson pmf table via recursion: p0=e^-mu, pk=pk-1*mu/k */
static void poisson_pmf_table(double mu, int K, double *p) {
    p[0] = exp(-mu);
    for (int k = 1; k <= K; k++) p[k] = p[k - 1] * mu / (double)k;
}

/* Sample Poisson(mu) using Knuth (mu=5 is small; good enough for SIM) */
static int poisson_sample(double mu) {
    double L = exp(-mu);
    int k = 0;
    double p = 1.0;
    do {
        k++;
        p *= (double)rand() / (double)RAND_MAX;
    } while (p > L);
    return k - 1;
}

/* Expected sales: E[min(D,Y)] with D~Pois(mu), Y=sum I.
   Stable: sum_{d=0..Y} d p_d + Y * P(D>Y) */
static double expected_sales(const int *I, int M, double mu) {
    int Y = 0;
    for (int r = 0; r < M; r++) Y += I[r];
    if (Y <= 0) return 0.0;

    double *p = (double *)malloc((size_t)(Y + 1) * sizeof(double));
    assert(p);
    poisson_pmf_table(mu, Y, p);

    double sum = 0.0, cdf = 0.0;
    for (int d = 0; d <= Y; d++) {
        cdf += p[d];
        sum += (double)d * p[d];
    }
    double tail = 1.0 - cdf;
    if (tail < 0.0) tail = 0.0;
    sum += (double)Y * tail;

    free(p);
    return sum;
}

/* FIFO transition F(q,I,d), generic M:
   - Oldest bucket I[0] has 1 period remaining; any leftover expires (not carried).
   - New order arrives as freshest bucket next[M-1] = q.
   If d >= Y: next = (0,0,...,0,q)
   Else:
     next[M-1]=q
     cum = 0
     for r=0..M-2:
        cum += I[r]
        shortfall = max(d - cum, 0)
        next[r] = max(I[r+1] - shortfall, 0)
*/
static int transition(int q, const int *I, int d, int Q, int M) {
    int Y = 0;
    for (int r = 0; r < M; r++) Y += I[r];

    int *next = (int *)alloca((size_t)M * sizeof(int));

    if (d >= Y) {
        for (int r = 0; r < M - 1; r++) next[r] = 0;
        next[M - 1] = q;
        return state_to_index(next, Q, M);
    }

    next[M - 1] = q;
    int cum = 0;
    for (int r = 0; r < M - 1; r++) {
        cum += I[r];
        int shortfall = (d > cum) ? (d - cum) : 0;
        int prev = I[r + 1];
        next[r] = (prev > shortfall) ? (prev - shortfall) : 0;
    }
    return state_to_index(next, Q, M);
}

/* Waste in this model (one-product): leftover of oldest bucket expires end-of-day.
   If demand d consumes FIFO from oldest: waste = max(I[0] - d, 0). */
static int waste_amount(const int *I, int M, int d) {
    (void)M;
    return (I[0] > d) ? (I[0] - d) : 0;
}

/* ---------- Value Iteration (average reward via span) ---------- */

typedef struct {
    int M;
    int Q;
    int N;
    double pi;
    int iters;
    int *best_q;  /* length N */
} VIResult;

/* Compute VI and store best action per state. */
static VIResult run_vi(int M, int Q, double mu, double s, double c, double eps) {
    VIResult res;
    res.M = M;
    res.Q = Q;
    res.N = state_space(Q, M);
    res.pi = 0.0;
    res.iters = 0;
    res.best_q = (int *)malloc((size_t)res.N * sizeof(int));
    assert(res.best_q);

    double *V = (double *)malloc((size_t)res.N * sizeof(double));
    double *W = (double *)malloc((size_t)res.N * sizeof(double));
    assert(V && W);

    int *I = (int *)malloc((size_t)M * sizeof(int));
    assert(I);

    /* init: V_j = s * E[min(D,Y)] */
    for (int j = 0; j < res.N; j++) {
        index_to_state(j, I, Q, M);
        V[j] = s * expected_sales(I, M, mu);
        res.best_q[j] = 0;
    }

    for (int iter = 0; iter < MAX_IT; iter++) {
        memcpy(W, V, (size_t)res.N * sizeof(double));

        for (int j = 0; j < res.N; j++) {
            index_to_state(j, I, Q, M);

            int Y = 0;
            for (int r = 0; r < M; r++) Y += I[r];

            double Esale = s * expected_sales(I, M, mu);

            /* Demand sum cutoff K:
               choose K >= Y so that tail transition is the same (d>=Y => fixed next state). */
            int K = Y + (int)(mu + 30); /* conservative for mu=5 */
            if (K < Y) K = Y;

            double *p = (double *)malloc((size_t)(K + 1) * sizeof(double));
            assert(p);
            poisson_pmf_table(mu, K, p);

            double best_val = -1e100;
            int bestq = 0;

            for (int q = 0; q <= Q; q++) {
                double sum = 0.0, cdf = 0.0;

                for (int d = 0; d <= K; d++) {
                    double pd = p[d];
                    cdf += pd;
                    int k = transition(q, I, d, Q, M);
                    sum += pd * W[k];
                }

                double tail = 1.0 - cdf;
                if (tail < 0.0) tail = 0.0;

                if (tail > 0.0) {
                    /* for d > K (and K>=Y), next state equals the d>=Y case => (0,...,0,q) */
                    int *next_over = (int *)alloca((size_t)M * sizeof(int));
                    for (int r = 0; r < M - 1; r++) next_over[r] = 0;
                    next_over[M - 1] = q;
                    int k_over = state_to_index(next_over, Q, M);
                    sum += tail * W[k_over];
                }

                double val = Esale + sum - c * (double)q;
                if (val > best_val) {
                    best_val = val;
                    bestq = q;
                }
            }

            free(p);

            V[j] = best_val;
            res.best_q[j] = bestq;
        }

        /* span convergence */
        double vmin = 1e100, vmax = -1e100;
        for (int j = 0; j < res.N; j++) {
            double d = V[j] - W[j];
            if (d < vmin) vmin = d;
            if (d > vmax) vmax = d;
        }
        double span = vmax - vmin;
        double pi_est = (vmin + vmax) / 2.0;

        printf("  iter %2d: span=%.6f pi=%.4f %s\n",
               iter + 1, span, pi_est, (span < eps ? "(converged)" : ""));

        if (span < eps) {
            res.iters = iter + 1;
            res.pi = pi_est;
            free(V); free(W); free(I);
            return res;
        }
    }

    /* not converged within MAX_IT, still return last pi estimate */
    double vmin = 1e100, vmax = -1e100;
    for (int j = 0; j < res.N; j++) {
        double d = V[j] - W[j];
        if (d < vmin) vmin = d;
        if (d > vmax) vmax = d;
    }
    res.iters = MAX_IT;
    res.pi = (vmin + vmax) / 2.0;

    free(V); free(W); free(I);
    return res;
}

/* ---------- Policy simulation ---------- */

typedef struct {
    double avg_profit;
    double waste_perc;  /* 100 * waste / total_order */
    double avg_order;
} SimResult;

/* policy_type: 0=base-stock, 1=optimal(VI policy) */
static SimResult simulate_policy(
    int M, int Q, double mu, double s, double c,
    int policy_type,
    const int *best_q, /* used if policy_type==1 */
    int base_S,        /* used if policy_type==0 */
    int T, int burn_in
) {
    SimResult out = {0};

    int *I = (int *)calloc((size_t)M, sizeof(int));
    assert(I);

    double total_profit = 0.0;
    long long total_order = 0;
    long long total_waste = 0;
    long long total_order_for_avg = 0;

    int N = state_space(Q, M);

    for (int t = 0; t < T + burn_in; t++) {
        int Y = 0;
        for (int r = 0; r < M; r++) Y += I[r];

        int q = 0;
        if (policy_type == 0) {
            int need = base_S - Y;
            q = (need > 0) ? need : 0;
            if (q > Q) q = Q;
        } else {
            int j = state_to_index(I, Q, M);
            assert(j >= 0 && j < N);
            q = best_q[j];
            if (q < 0) q = 0;
            if (q > Q) q = Q;
        }

        int d = poisson_sample(mu);
        int sales = (d < Y) ? d : Y;
        int waste = waste_amount(I, M, d);

        double profit = s * (double)sales - c * (double)q;

        int next_index = transition(q, I, d, Q, M);
        index_to_state(next_index, I, Q, M);

        if (t >= burn_in) {
            total_profit += profit;
            total_order += q;
            total_waste += waste;
            total_order_for_avg += q;
        }
    }

    out.avg_profit = total_profit / (double)T;
    out.avg_order = (double)total_order_for_avg / (double)T;
    if (total_order > 0) out.waste_perc = 100.0 * (double)total_waste / (double)total_order;
    else out.waste_perc = 0.0;

    free(I);
    return out;
}

/* ---------- Experiments from paper 2.3 ---------- */

static void run_experiment(int M, int Q, int base_S) {
    printf("\n============================================================\n");
    printf("Experiment: M=%d, Qbar=%d, base-stock S=%d\n", M, Q, base_S);
    printf("Params: s=%.1f c=%.1f mu=%.0f eps=%.1e\n",
           S_PRICE, C_COST, MU_DEMAND, EPS);
    printf("State space N=(Q+1)^M = %d\n\n", state_space(Q, M));

    clock_t t0 = clock();
    VIResult vi = run_vi(M, Q, MU_DEMAND, S_PRICE, C_COST, EPS);
    double sec = (double)(clock() - t0) / CLOCKS_PER_SEC;

    int maxq = 0;
    for (int j = 0; j < vi.N; j++) if (vi.best_q[j] > maxq) maxq = vi.best_q[j];

    printf("\nVI done: iters=%d pi=%.4f time=%.3fs\n", vi.iters, vi.pi, sec);
    printf("Policy check: max optimal q over all states = %d\n", maxq);

    /* Simulate base-stock and optimal policy like paper (400k periods) */
    SimResult base = simulate_policy(M, Q, MU_DEMAND, S_PRICE, C_COST, 0, NULL, base_S, SIM_T, BURN_IN);
    SimResult opt  = simulate_policy(M, Q, MU_DEMAND, S_PRICE, C_COST, 1, vi.best_q, base_S, SIM_T, BURN_IN);

    printf("\nSimulation (%d periods, burn-in %d):\n", SIM_T, BURN_IN);
    printf("  Base-stock:  Pi=%.3f, Wasteperc=%.2f, AvgOrder=%.3f\n",
           base.avg_profit, base.waste_perc, base.avg_order);
    printf("  Optimal   :  Pi=%.3f, Wasteperc=%.2f, AvgOrder=%.3f\n",
           opt.avg_profit, opt.waste_perc, opt.avg_order);

    free(vi.best_q);
}

int main(void) {
    srand(1); /* fixed seed for reproducibility */

    printf("One-product perishable inventory VI + simulation (Paper 2.3)\n");
    printf("Paper reports (M=2): base S=13, Q=9, 12 iters, pi~2.215, wasteperc drops.\n");
    printf("Then (M=3): base S=15, N=4096, 15 iters.\n");
    printf("Then (M=4): Q=20, pi~2.47, base S=16 practically same profit.\n");

    /* Paper 2.3 experiments */
    run_experiment(2,  9, 13);
    run_experiment(3, 15, 15);
    run_experiment(4, 20, 16);

    return 0;
}