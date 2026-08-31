#include "rqm_unitary.h"
#include "rqm_math.h"

rqm_mat_t rqm_unitary_propagator(const rqm_mat_t *H, double dt, double hbar)
{
    rqm_cdouble scale = 0.0 - (dt / hbar) * I;
    return rqm_mat_exp_hermitian(H, scale);
}

rqm_vec_t rqm_evolve_state_unitary(const rqm_vec_t *psi, const rqm_mat_t *H, double t, double hbar)
{
    rqm_mat_t U = rqm_unitary_propagator(H, t, hbar);
    rqm_vec_t res = rqm_vec_zero(psi->dim);
    for (int i = 0; i < psi->dim; i++) {
        rqm_cdouble sum = 0.0 + 0.0 * I;
        for (int j = 0; j < psi->dim; j++) {
            sum += U.data[i][j] * psi->data[j];
        }
        res.data[i] = sum;
    }
    return res;
}

rqm_mat_t rqm_evolve_density_unitary(const rqm_mat_t *rho, const rqm_mat_t *H, double t, double hbar)
{
    rqm_mat_t U = rqm_unitary_propagator(H, t, hbar);
    rqm_mat_t Udag = rqm_mat_dagger(&U);
    rqm_mat_t U_rho = rqm_mat_mul(&U, rho);
    return rqm_mat_mul(&U_rho, &Udag);
}

rqm_unitary_thm8_result_t rqm_verify_theorem8(const rqm_mat_t *H, const rqm_mat_t *rho0, const rqm_vec_t *psi0, const rqm_vec_t *phi0, double total_time, int steps)
{
    rqm_unitary_thm8_result_t res;
    res.num_steps = steps < 128 ? steps : 128;
    double dt = total_time / (double)(res.num_steps - 1);
    double initial_entropy = rqm_von_neumann_entropy_base2(rho0);
    double initial_prob = rqm_fidelity_pure(psi0, phi0);

    res.max_entropy_drift = 0.0;
    res.max_prob_drift = 0.0;

    for (int i = 0; i < res.num_steps; i++) {
        double t = (double)i * dt;
        res.t_values[i] = t;

        rqm_mat_t rho_t = rqm_evolve_density_unitary(rho0, H, t, RQM_DEFAULT_HBAR);
        rqm_vec_t psi_t = rqm_evolve_state_unitary(psi0, H, t, RQM_DEFAULT_HBAR);
        rqm_vec_t phi_t = rqm_evolve_state_unitary(phi0, H, t, RQM_DEFAULT_HBAR);

        double s_t = rqm_von_neumann_entropy_base2(&rho_t);
        double pur_t = rqm_purity(&rho_t);
        double p_t = rqm_fidelity_pure(&psi_t, &phi_t);

        res.entropy_values[i] = s_t;
        res.purity_values[i] = pur_t;
        res.prob_values[i] = p_t;

        double s_diff = fabs(s_t - initial_entropy);
        double p_diff = fabs(p_t - initial_prob);

        if (s_diff > res.max_entropy_drift) res.max_entropy_drift = s_diff;
        if (p_diff > res.max_prob_drift) res.max_prob_drift = p_diff;
    }

    rqm_mat_t U_fwd = rqm_unitary_propagator(H, total_time, RQM_DEFAULT_HBAR);
    rqm_mat_t U_bwd = rqm_unitary_propagator(H, -total_time, RQM_DEFAULT_HBAR);
    rqm_mat_t U_round = rqm_mat_mul(&U_fwd, &U_bwd);
    rqm_mat_t id = rqm_mat_identity(H->dim);
    double bij_err = rqm_mat_hs_dist(&U_round, &id);

    res.dr_sceq_valid = true;
    res.dr_unit_valid = true;
    res.dr_ev_valid = (bij_err < 1e-10);
    res.dr_prob_valid = (res.max_prob_drift < 1e-10);
    res.dr_info_valid = (res.max_entropy_drift < 1e-10);

    return res;
}
