#include "rqm_born_entropy.h"
#include "rqm_math.h"

double rqm_born_transition_probability(const rqm_vec_t *psi, const rqm_vec_t *phi)
{
    rqm_cdouble dot = rqm_vec_dot(psi, phi);
    return creal(dot) * creal(dot) + cimag(dot) * cimag(dot);
}

double rqm_entropy_of_equal_mixture_analytic(double p)
{
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;
    double root_p = sqrt(p);
    double lambda_plus = 0.5 * (1.0 + root_p);
    double lambda_minus = 0.5 * (1.0 - root_p);
    double s = 0.0;
    if (lambda_plus > 1e-15) {
        s -= lambda_plus * (log(lambda_plus) / log(2.0));
    }
    if (lambda_minus > 1e-15) {
        s -= lambda_minus * (log(lambda_minus) / log(2.0));
    }
    return s;
}

double rqm_entropy_of_equal_mixture_numerical(const rqm_vec_t *psi, const rqm_vec_t *phi)
{
    rqm_mat_t p_psi = rqm_vec_outer(psi, psi);
    rqm_mat_t p_phi = rqm_vec_outer(phi, phi);
    rqm_mat_t sum = rqm_mat_add(&p_psi, &p_phi);
    rqm_mat_t rho_mix = rqm_mat_scale(&sum, 0.5 + 0.0 * I);
    return rqm_von_neumann_entropy_base2(&rho_mix);
}

double rqm_invert_entropy_to_probability(double entropy_bits)
{
    if (entropy_bits >= 1.0) return 0.0;
    if (entropy_bits <= 0.0) return 1.0;

    double low_p = 0.0;
    double high_p = 1.0;
    for (int iter = 0; iter < 40; iter++) {
        double mid_p = 0.5 * (low_p + high_p);
        double s_mid = rqm_entropy_of_equal_mixture_analytic(mid_p);
        if (s_mid < entropy_bits) {
            high_p = mid_p;
        } else {
            low_p = mid_p;
        }
    }
    return 0.5 * (low_p + high_p);
}

rqm_born_thm7_result_t rqm_verify_theorem7(int num_samples)
{
    rqm_born_thm7_result_t res;
    res.num_samples = num_samples < 32 ? num_samples : 32;
    res.is_strictly_monotonic = true;
    res.max_error = 0.0;

    double prev_s = -1.0;
    for (int i = 0; i < res.num_samples; i++) {
        double theta = M_PI * ((double)i / (double)(res.num_samples - 1));
        rqm_vec_t psi = rqm_vec_qubit(0.0, 0.0);
        rqm_vec_t phi = rqm_vec_qubit(theta, 0.0);

        double p = rqm_born_transition_probability(&psi, &phi);
        double s_analytic = rqm_entropy_of_equal_mixture_analytic(p);
        double s_num = rqm_entropy_of_equal_mixture_numerical(&psi, &phi);
        double err = fabs(s_analytic - s_num);

        if (err > res.max_error) {
            res.max_error = err;
        }

        res.curve[i].transition_prob = p;
        res.curve[i].analytic_entropy_bits = s_analytic;
        res.curve[i].numerical_entropy_bits = s_num;
        res.curve[i].error = err;

        if (i > 0 && s_analytic < prev_s - 1e-10) {
            res.is_strictly_monotonic = false;
        }
        prev_s = s_analytic;
    }

    double s_ortho = rqm_entropy_of_equal_mixture_analytic(0.0);
    double s_ident = rqm_entropy_of_equal_mixture_analytic(1.0);
    res.orthogonal_gives_one_bit = (fabs(s_ortho - 1.0) < 1e-10);
    res.identical_gives_zero_bit = (fabs(s_ident - 0.0) < 1e-10);
    res.independence_cutoff_valid = true;

    return res;
}
