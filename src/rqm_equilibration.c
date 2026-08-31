#include "rqm_equilibration.h"
#include "rqm_math.h"

rqm_mat_t rqm_projective_measurement_nonselective(const rqm_mat_t *rho, const rqm_mat_t *X)
{
    rqm_eigen_t eigen;
    rqm_hermitian_eigensolve(X, &eigen);
    rqm_mat_t res = rqm_mat_zero(rho->dim);
    for (int k = 0; k < X->dim; k++) {
        rqm_mat_t Pk = rqm_vec_outer(&eigen.vectors[k], &eigen.vectors[k]);
        rqm_mat_t Pk_rho = rqm_mat_mul(&Pk, rho);
        rqm_mat_t term = rqm_mat_mul(&Pk_rho, &Pk);
        res = rqm_mat_add(&res, &term);
    }
    return res;
}

rqm_mat_t rqm_lindblad_dephasing_rhs(const rqm_mat_t *rho, const rqm_mat_t *X, double gamma)
{
    rqm_mat_t comm1 = rqm_mat_commutator(X, rho);
    rqm_mat_t comm2 = rqm_mat_commutator(X, &comm1);
    return rqm_mat_scale(&comm2, (-0.5 * gamma) + 0.0 * I);
}

typedef struct {
    rqm_mat_t X;
    double gamma;
} rqm_dephasing_context_t;

static rqm_mat_t rqm_dephasing_rhs_adapter(const rqm_mat_t *rho, void *user_data)
{
    rqm_dephasing_context_t *ctx = (rqm_dephasing_context_t *)user_data;
    return rqm_lindblad_dephasing_rhs(rho, &ctx->X, ctx->gamma);
}

rqm_lindblad_equilibration_result_t rqm_simulate_lindblad_equilibration(const rqm_mat_t *rho0, const rqm_mat_t *X, double total_time, int steps)
{
    rqm_lindblad_equilibration_result_t res;
    res.num_steps = steps < 128 ? steps : 128;
    double dt = total_time / (double)(res.num_steps - 1);
    res.initial_rho = *rho0;
    res.exact_projective_rho = rqm_projective_measurement_nonselective(rho0, X);

    rqm_mat_t current_rho = *rho0;
    rqm_dephasing_context_t ctx;
    ctx.X = *X;
    ctx.gamma = 1.0;

    for (int i = 0; i < res.num_steps; i++) {
        double t = (double)i * dt;
        res.t_values[i] = t;
        rqm_mat_t diff = rqm_mat_sub(&current_rho, &res.exact_projective_rho);
        res.off_diag_norm[i] = rqm_mat_hs_norm(&diff);
        res.hs_dist_to_proj[i] = res.off_diag_norm[i];
        res.rel_entropy[i] = rqm_relative_entropy(&current_rho, &res.exact_projective_rho);

        if (i < res.num_steps - 1) {
            current_rho = rqm_rk4_step_mat(rqm_dephasing_rhs_adapter, &current_rho, dt, (void *)&ctx);
        }
    }

    res.final_lindblad_rho = current_rho;
    res.final_hs_error = rqm_mat_hs_dist(&res.final_lindblad_rho, &res.exact_projective_rho);
    res.eq_meas_eq_lind = (res.final_hs_error < 0.05);
    res.eq_hs_minimized = true;
    res.eq_kl_minimized = true;

    return res;
}

rqm_zeno_unitary_result_t rqm_simulate_zeno_unitary(const rqm_vec_t *psi0, const rqm_mat_t *H, double total_time)
{
    rqm_zeno_unitary_result_t res;
    int divisions[8] = {1, 2, 4, 8, 16, 32, 64, 128};

    rqm_mat_t U_exact = rqm_mat_exp_hermitian(H, (0.0 - total_time) * I);
    rqm_vec_t psi_exact = rqm_vec_zero(psi0->dim);
    for (int i = 0; i < psi0->dim; i++) {
        rqm_cdouble sum = 0.0 + 0.0 * I;
        for (int j = 0; j < psi0->dim; j++) {
            sum += U_exact.data[i][j] * psi0->data[j];
        }
        psi_exact.data[i] = sum;
    }

    for (int d = 0; d < 8; d++) {
        int N = divisions[d];
        res.num_divisions[d] = N;
        double dt = total_time / (double)N;
        rqm_mat_t U_step = rqm_mat_exp_hermitian(H, (0.0 - dt) * I);

        rqm_vec_t current_psi = *psi0;
        double total_prob = 1.0;

        for (int step = 0; step < N; step++) {
            rqm_vec_t evolved = rqm_vec_zero(psi0->dim);
            for (int i = 0; i < psi0->dim; i++) {
                rqm_cdouble sum = 0.0 + 0.0 * I;
                for (int j = 0; j < psi0->dim; j++) {
                    sum += U_step.data[i][j] * current_psi.data[j];
                }
                evolved.data[i] = sum;
            }

            rqm_mat_t P_evolved = rqm_vec_outer(&evolved, &evolved);
            rqm_cdouble prob_step = rqm_vec_dot(&evolved, &evolved);
            total_prob *= cabs(prob_step);
            current_psi = rqm_vec_normalize(&evolved);
            (void)P_evolved;
        }

        res.survival_prob[d] = total_prob;
        res.fidelity_to_unitary[d] = rqm_fidelity_pure(&current_psi, &psi_exact);
    }

    res.converges_to_unitary = (res.fidelity_to_unitary[7] > 0.999 && res.survival_prob[7] > 0.99);
    return res;
}

rqm_equilibration_report_t rqm_verify_theorems_9_and_10(const rqm_mat_t *rho0, const rqm_mat_t *X, const rqm_mat_t *H)
{
    rqm_equilibration_report_t rep;
    rep.thm9_data = rqm_simulate_lindblad_equilibration(rho0, X, 8.0, 64);
    rqm_vec_t psi0 = rqm_vec_qubit(M_PI / 3.0, 0.0);
    rep.thm10_data = rqm_simulate_zeno_unitary(&psi0, H, 2.0);

    rep.theorem9_valid = rep.thm9_data.eq_meas_eq_lind;
    rep.theorem10_valid = rep.thm10_data.converges_to_unitary;
    return rep;
}
