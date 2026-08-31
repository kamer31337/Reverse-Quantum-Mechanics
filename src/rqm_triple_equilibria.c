#include "rqm_triple_equilibria.h"
#include "rqm_math.h"

rqm_dynamic_equilibrium_t rqm_verify_dynamic_equilibrium(const rqm_mat_t *rho)
{
    rqm_dynamic_equilibrium_t res;
    res.H_dynamic = *rho;
    rqm_mat_t comm = rqm_mat_commutator(&res.H_dynamic, rho);
    res.commutator_norm = rqm_mat_hs_norm(&comm);
    res.is_dynamic_equilibrium = (res.commutator_norm < 1e-12);
    return res;
}

rqm_spectral_equilibrium_t rqm_verify_spectral_equilibrium(const rqm_mat_t *rho)
{
    rqm_spectral_equilibrium_t res;
    res.X_observable = *rho;
    rqm_eigen_t eigen;
    rqm_hermitian_eigensolve(&res.X_observable, &eigen);

    rqm_mat_t proj_sum = rqm_mat_zero(rho->dim);
    for (int k = 0; k < rho->dim; k++) {
        rqm_mat_t Pk = rqm_vec_outer(&eigen.vectors[k], &eigen.vectors[k]);
        rqm_mat_t Pk_rho = rqm_mat_mul(&Pk, rho);
        rqm_mat_t term = rqm_mat_mul(&Pk_rho, &Pk);
        proj_sum = rqm_mat_add(&proj_sum, &term);
    }

    res.projection_output = proj_sum;
    res.measurement_diff_norm = rqm_mat_hs_dist(rho, &res.projection_output);
    res.is_spectral_equilibrium = (res.measurement_diff_norm < 1e-12);
    return res;
}

rqm_thermodynamic_equilibrium_t rqm_verify_thermodynamic_equilibrium(const rqm_mat_t *rho, double beta)
{
    rqm_thermodynamic_equilibrium_t res;
    res.beta = beta;
    rqm_eigen_t eigen;
    rqm_hermitian_eigensolve(rho, &eigen);

    bool full_rank = true;
    for (int i = 0; i < rho->dim; i++) {
        if (eigen.values[i] < 1e-6) {
            full_rank = false;
            break;
        }
    }
    res.is_full_rank = full_rank;

    if (full_rank) {
        rqm_mat_t log_rho = rqm_mat_log_density(rho, 1e-14);
        res.H_thermodynamic = rqm_mat_scale(&log_rho, (-1.0 / beta) + 0.0 * I);
        res.exact_or_approx_gibbs = rqm_mat_exp_hermitian(&res.H_thermodynamic, (-beta) + 0.0 * I);
        rqm_cdouble z = rqm_mat_trace(&res.exact_or_approx_gibbs);
        res.exact_or_approx_gibbs = rqm_mat_scale(&res.exact_or_approx_gibbs, (1.0 / creal(z)) + 0.0 * I);
        res.gibbs_error_norm = rqm_mat_hs_dist(rho, &res.exact_or_approx_gibbs);
        res.num_k_steps = 0;
        res.is_thermodynamic_equilibrium = (res.gibbs_error_norm < 1e-6);
    } else {
        res.num_k_steps = 6;
        for (int step = 0; step < res.num_k_steps; step++) {
            double k_pen = pow(10.0, (double)(step + 1));
            res.k_penalties[step] = k_pen;

            rqm_mat_t Hk = rqm_mat_zero(rho->dim);
            for (int i = 0; i < rho->dim; i++) {
                double energy = 0.0;
                if (eigen.values[i] > 1e-6) {
                    energy = -log(eigen.values[i]) / beta;
                } else {
                    energy = k_pen * (double)((i + 1) * (i + 1));
                }
                rqm_mat_t Pk = rqm_vec_outer(&eigen.vectors[i], &eigen.vectors[i]);
                rqm_mat_t term = rqm_mat_scale(&Pk, energy + 0.0 * I);
                Hk = rqm_mat_add(&Hk, &term);
            }

            rqm_mat_t unnorm_gibbs = rqm_mat_exp_hermitian(&Hk, (-beta) + 0.0 * I);
            rqm_cdouble z = rqm_mat_trace(&unnorm_gibbs);
            rqm_mat_t gibbs_k = rqm_mat_scale(&unnorm_gibbs, (1.0 / creal(z)) + 0.0 * I);

            double err = rqm_mat_hs_dist(rho, &gibbs_k);
            res.k_errors[step] = err;

            if (step == res.num_k_steps - 1) {
                res.H_thermodynamic = Hk;
                res.exact_or_approx_gibbs = gibbs_k;
                res.gibbs_error_norm = err;
            }
        }
        res.is_thermodynamic_equilibrium = (res.gibbs_error_norm < 1e-4);
    }

    return res;
}

rqm_triple_equilibria_result_t rqm_verify_theorem12(const rqm_mat_t *rho, double beta)
{
    rqm_triple_equilibria_result_t res;
    res.target_state = *rho;
    res.dynamic_eq = rqm_verify_dynamic_equilibrium(rho);
    res.spectral_eq = rqm_verify_spectral_equilibrium(rho);
    res.thermo_eq = rqm_verify_thermodynamic_equilibrium(rho, beta);
    res.theorem12_all_equilibria_verified = res.dynamic_eq.is_dynamic_equilibrium && res.spectral_eq.is_spectral_equilibrium && res.thermo_eq.is_thermodynamic_equilibrium;
    return res;
}
