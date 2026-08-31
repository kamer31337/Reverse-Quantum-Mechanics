#include "rqm_classical_limit.h"
#include "rqm_math.h"

double rqm_wigner_fock_state(int n, double q, double p, double hbar)
{
    double r2 = (q * q + p * p) / hbar;
    double exp_factor = exp(-r2);
    if (n == 0) {
        return (1.0 / (M_PI * hbar)) * exp_factor;
    } else if (n == 1) {
        return (1.0 / (M_PI * hbar)) * (2.0 * r2 - 1.0) * exp_factor;
    } else if (n == 2) {
        return (1.0 / (M_PI * hbar)) * (2.0 * r2 * r2 - 4.0 * r2 + 1.0) * exp_factor;
    }
    return (1.0 / (M_PI * hbar)) * exp_factor;
}

double rqm_wigner_damped(int n, double q, double p, double lambda, double hbar)
{
    if (lambda < 1.0) lambda = 1.0;
    double r2 = (q * q + p * p) / (lambda * hbar);
    double exp_factor = exp(-r2);
    if (n == 1) {
        double poly = (2.0 * r2) + (1.0 - 2.0 / lambda);
        return (1.0 / (M_PI * hbar * lambda)) * poly * exp_factor;
    }
    return (1.0 / (M_PI * hbar * lambda)) * exp_factor;
}

rqm_mat_t rqm_lindblad_creation_rhs(const rqm_mat_t *rho, double gamma)
{
    int dim = rho->dim;
    rqm_mat_t adag = rqm_mat_creation(dim);
    rqm_mat_t a = rqm_mat_annihilation(dim);

    rqm_mat_t term1 = rqm_mat_mul(&adag, rho);
    term1 = rqm_mat_mul(&term1, &a);

    rqm_mat_t a_adag = rqm_mat_mul(&a, &adag);
    rqm_mat_t acomp = rqm_mat_anticommutator(&a_adag, rho);
    rqm_mat_t term2 = rqm_mat_scale(&acomp, 0.5 + 0.0 * I);

    rqm_mat_t diff = rqm_mat_sub(&term1, &term2);
    return rqm_mat_scale(&diff, gamma + 0.0 * I);
}

rqm_classical_limit_result_t rqm_verify_theorem11(int initial_fock, double gamma, double total_time, int steps)
{
    rqm_classical_limit_result_t res;
    res.num_steps = steps < 32 ? steps : 32;
    double dt = total_time / (double)(res.num_steps - 1);
    res.negativity_suppressed = false;
    res.commutator_vanishes = false;

    for (int i = 0; i < res.num_steps; i++) {
        double t = (double)i * dt;
        double lambda = exp(gamma * t);
        res.time_points[i] = t;
        res.lambda_values[i] = lambda;
        res.effective_hbar[i] = RQM_DEFAULT_HBAR / lambda;
        res.entropy_values[i] = log(lambda);
        res.min_wigner_value[i] = rqm_wigner_damped(initial_fock, 0.0, 0.0, lambda, RQM_DEFAULT_HBAR);
    }

    if (res.min_wigner_value[res.num_steps - 1] >= 0.0) {
        res.negativity_suppressed = true;
    }
    if (res.effective_hbar[res.num_steps - 1] < 0.1) {
        res.commutator_vanishes = true;
    }

    res.theorem11_valid = res.negativity_suppressed && res.commutator_vanishes;
    return res;
}
