/**
 * vi_sequential.c - Algorithm 3: Sequential VI for two-product substitution
 * Paper: Ortega et al. (2018), Algorithm 3
 */

#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

/* Run VI for one instance. Returns iterations to converge, writes pi to *out_pi.
 * name: instance name for progress output (e.g. "P1") */
static int run_vi(const Instance *inst, const char *name, double *out_pi, int max_iter) {
    int N = inst->N, Qa = inst->Qa, Qb = inst->Qb;
    double mu_a = (double)inst->mu_a, mu_b = (double)inst->mu_b;

    double *V = (double *)malloc((size_t)N * sizeof(double));
    double *W = (double *)malloc((size_t)N * sizeof(double));
    assert(V && W);

    /* Algorithm 3 line 1: V_j = Esale_j for j = 0,...,N-1 */
    int init_step = (N >= 5000) ? 2000 : (N >= 500 ? 100 : 50);
    for (int j = 0; j < N; j++) {
        V[j] = expected_sales(j, Qa, Qb, inst->mu_a, inst->mu_b, GAMMA);
        if (j > 0 && j % init_step == 0) { printf("[%s] init %d/%d\n", name, j, N); fflush(stdout); }
    }

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        memcpy(W, V, (size_t)N * sizeof(double));

        /* Line 4-5: j=0 empty state. k = (0,0,qa,0,0,qb). Bellman: V0 = max [W[k]-ca*qa-cb*qb] */
        State s0;
        memset(&s0, 0, sizeof(State));
        double best_v0 = -1e99;
        for (int qa = 0; qa < Qa; qa++) {
            for (int qb = 0; qb < Qb; qb++) {
                s0.Ia[M - 1] = qa;
                s0.Ib[M - 1] = qb;
                int k = state_to_index(&s0, Qa, Qb);
                double val = W[k] - CA * qa - CB * qb;
                if (val > best_v0) best_v0 = val;
            }
        }
        V[0] = best_v0;

        /* Line 6: for j = 0,...,N-1 */
        int state_step = (N >= 5000) ? 2000 : (N >= 500 ? 100 : 50);
        for (int j = 1; j < N; j++) {
            if (j % state_step == 0 || j == N - 1) {
                printf("[%s] iter %d state %d/%d (%.0f%%)\n", name, iter + 1, j, N, 100.0 * j / N);
                fflush(stdout);
            }
            State s;
            index_to_state(j, Qa, Qb, &s);
            int Ya = 0, Yb = 0;
            for (int r = 0; r < M; r++) {
                Ya += s.Ia[r];
                Yb += s.Ib[r];
            }
            double Esale = expected_sales(j, Qa, Qb, inst->mu_a, inst->mu_b, GAMMA);

            if (Ya == 0) {
                /* Line 10-12: No substitution. k = F(qa,qb,(0,0,Ib), db) */
                double best_Fv = -1e99;
                for (int qa = 0; qa < Qa; qa++) {
                    for (int qb = 0; qb < Qb; qb++) {
                        double sum = 0.0;
                        int db_lim = Yb + (int)(mu_b + 20);
                        for (int db = 0; db <= db_lim; db++) {
                            double p = poisson_pmf(db, mu_b);
                            if (p < 1e-12 && db > Yb) break;
                            int k = transition_eff(qa, qb, &s, 0, db, Qa, Qb);
                            sum += p * W[k];
                        }
                        double Fv = Esale + sum - CA * qa - CB * qb;
                        if (Fv > best_Fv) best_Fv = Fv;
                    }
                }
                V[j] = best_Fv;
            } else {
                /* Line 13-19: Substitution. Split: db<=Yb (da,db) loop; db>Yb use pz aggregation */
                double best_Fv = -1e99;
                for (int qa = 0; qa < Qa; qa++) {
                    for (int qb = 0; qb < Qb; qb++) {
                        double sum = 0.0;
                        for (int db = 0; db <= Yb; db++) {
                            double p_db = poisson_pmf(db, mu_b);
                            for (int da = 0; da <= Ya + 30; da++) {
                                double p_da = poisson_pmf(da, mu_a);
                                if (p_da < 1e-12 && da > Ya) break;
                                int d_a = (da < Ya) ? da : Ya;
                                int k = transition_eff(qa, qb, &s, d_a, db, Qa, Qb);
                                sum += p_da * p_db * W[k];
                            }
                        }
                        sum += substitution_future_value(qa, qb, &s, Qa, Qb, mu_a, mu_b, GAMMA, W);
                        double Fv = Esale + sum - CA * qa - CB * qb;
                        if (Fv > best_Fv) best_Fv = Fv;
                    }
                }
                V[j] = best_Fv;
            }
        }

        /* Span check */
        double vmin = 1e99, vmax = -1e99;
        for (int j = 0; j < N; j++) {
            double d = V[j] - W[j];
            if (d < vmin) vmin = d;
            if (d > vmax) vmax = d;
        }
        double span = vmax - vmin;
        double pi_est = (vmin + vmax) / 2.0;
        printf("[%s] iter %d done: span=%.6f pi~%.4f %s\n", name, iter + 1, span, pi_est,
               span < EPSILON ? "(converged)" : "");
        fflush(stdout);
        if (span < EPSILON) {
            *out_pi = pi_est;
            free(V);
            free(W);
            return iter + 1;
        }
    }

    double vmin = 1e99, vmax = -1e99;
    for (int j = 0; j < N; j++) {
        double d = V[j] - W[j];
        if (d < vmin) vmin = d;
        if (d > vmax) vmax = d;
    }
    *out_pi = (vmin + vmax) / 2.0;
    free(V);
    free(W);
    return iter;
}

int main(void) {
    printf("CUDA Perishable Inventory VI - Sequential (Algorithm 3)\n");
    printf("Instances: P1, P2, P3, P4 (checkpoint: results/P*_result.txt)\n\n");

    system("mkdir -p results");

    Instance instances[] = {
        INSTANCE_P1, INSTANCE_P2, INSTANCE_P3, INSTANCE_P4
    };
    const char *names[] = { "P1", "P2", "P3", "P4" };

    for (int i = 0; i < 4; i++) {
        Instance inst = instances[i];
        inst.N = state_space_size(inst.Qa, inst.Qb);

        char resultpath[128];
        snprintf(resultpath, sizeof(resultpath), "results/%s_result.txt", names[i]);

        FILE *exist = fopen(resultpath, "r");
        if (exist) {
            int iters;
            double pi, sec;
            if (fscanf(exist, "iter=%d pi=%lf time=%lf", &iters, &pi, &sec) == 3) {
                fclose(exist);
                printf("[%s] checkpoint: iter=%d pi=%.4f time=%.2fs (skip)\n\n", names[i], iters, pi, sec);
                fflush(stdout);
                continue;
            }
            fclose(exist);
        }

        printf("=== %s (N=%d) ===\n", names[i], inst.N);
        fflush(stdout);
        double pi;
        clock_t t0 = clock();
        int iters = run_vi(&inst, names[i], &pi, MAX_ITER);
        double sec = (double)(clock() - t0) / CLOCKS_PER_SEC;

        FILE *fres = fopen(resultpath, "w");
        if (fres) {
            fprintf(fres, "iter=%d pi=%.6f time=%.2f\n", iters, pi, sec);
            fclose(fres);
        }

        printf("%s: mu_a=%d mu_b=%d Qa=%d Qb=%d N=%d\n", names[i], inst.mu_a, inst.mu_b, inst.Qa, inst.Qb, inst.N);
        printf("  iter=%d pi=%.4f time=%.2fs -> %s\n\n", iters, pi, sec, resultpath);
        fflush(stdout);
    }
    printf("Done. Checkpoints in results/P*_result.txt\n");
    return 0;
}
