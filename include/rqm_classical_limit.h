#ifndef RQM_CLASSICAL_LIMIT_H
#define RQM_CLASSICAL_LIMIT_H

#include "rqm_common.h"

typedef struct {
    int num_steps;
    double time_points[32];
    double lambda_values[32];
    double effective_hbar[32];
    double entropy_values[32];
    double min_wigner_value[32];
    bool negativity_suppressed;
    bool commutator_vanishes;
    bool theorem11_valid;
} rqm_classical_limit_result_t;

double rqm_wigner_fock_state(int n, double q, double p, double hbar);
double rqm_wigner_damped(int n, double q, double p, double lambda, double hbar);
rqm_mat_t rqm_lindblad_creation_rhs(const rqm_mat_t *rho, double gamma);
rqm_classical_limit_result_t rqm_verify_theorem11(int initial_fock, double gamma, double total_time, int steps);

#endif
