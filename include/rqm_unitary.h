#ifndef RQM_UNITARY_H
#define RQM_UNITARY_H

#include "rqm_common.h"

typedef struct {
    int num_steps;
    double t_values[128];
    double entropy_values[128];
    double purity_values[128];
    double prob_values[128];
    bool dr_sceq_valid;
    bool dr_unit_valid;
    bool dr_ev_valid;
    bool dr_prob_valid;
    bool dr_info_valid;
    double max_entropy_drift;
    double max_prob_drift;
} rqm_unitary_thm8_result_t;

rqm_mat_t rqm_unitary_propagator(const rqm_mat_t *H, double dt, double hbar);
rqm_vec_t rqm_evolve_state_unitary(const rqm_vec_t *psi, const rqm_mat_t *H, double t, double hbar);
rqm_mat_t rqm_evolve_density_unitary(const rqm_mat_t *rho, const rqm_mat_t *H, double t, double hbar);
rqm_unitary_thm8_result_t rqm_verify_theorem8(const rqm_mat_t *H, const rqm_mat_t *rho0, const rqm_vec_t *psi0, const rqm_vec_t *phi0, double total_time, int steps);

#endif
