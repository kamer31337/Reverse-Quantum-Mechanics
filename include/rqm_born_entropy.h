#ifndef RQM_BORN_ENTROPY_H
#define RQM_BORN_ENTROPY_H

#include "rqm_common.h"

typedef struct {
    double transition_prob;
    double analytic_entropy_bits;
    double numerical_entropy_bits;
    double error;
} rqm_entropy_curve_point_t;

typedef struct {
    int num_samples;
    rqm_entropy_curve_point_t curve[32];
    bool is_strictly_monotonic;
    bool orthogonal_gives_one_bit;
    bool identical_gives_zero_bit;
    bool independence_cutoff_valid;
    double max_error;
} rqm_born_thm7_result_t;

double rqm_born_transition_probability(const rqm_vec_t *psi, const rqm_vec_t *phi);
double rqm_entropy_of_equal_mixture_analytic(double p);
double rqm_entropy_of_equal_mixture_numerical(const rqm_vec_t *psi, const rqm_vec_t *phi);
double rqm_invert_entropy_to_probability(double entropy_bits);
rqm_born_thm7_result_t rqm_verify_theorem7(int num_samples);

#endif
