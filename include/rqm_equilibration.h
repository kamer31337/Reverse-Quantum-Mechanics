#ifndef RQM_EQUILIBRATION_H
#define RQM_EQUILIBRATION_H

#include "rqm_common.h"

typedef struct {
    int num_steps;
    double t_values[128];
    double off_diag_norm[128];
    double hs_dist_to_proj[128];
    double rel_entropy[128];
    rqm_mat_t initial_rho;
    rqm_mat_t final_lindblad_rho;
    rqm_mat_t exact_projective_rho;
    bool eq_meas_eq_lind;
    bool eq_hs_minimized;
    bool eq_kl_minimized;
    double final_hs_error;
} rqm_lindblad_equilibration_result_t;

typedef struct {
    int num_divisions[8];
    double survival_prob[8];
    double fidelity_to_unitary[8];
    bool converges_to_unitary;
} rqm_zeno_unitary_result_t;

typedef struct {
    bool theorem9_valid;
    bool theorem10_valid;
    rqm_lindblad_equilibration_result_t thm9_data;
    rqm_zeno_unitary_result_t thm10_data;
} rqm_equilibration_report_t;

rqm_mat_t rqm_projective_measurement_nonselective(const rqm_mat_t *rho, const rqm_mat_t *X);
rqm_mat_t rqm_lindblad_dephasing_rhs(const rqm_mat_t *rho, const rqm_mat_t *X, double gamma);
rqm_lindblad_equilibration_result_t rqm_simulate_lindblad_equilibration(const rqm_mat_t *rho0, const rqm_mat_t *X, double total_time, int steps);
rqm_zeno_unitary_result_t rqm_simulate_zeno_unitary(const rqm_vec_t *psi0, const rqm_mat_t *H, double total_time);
rqm_equilibration_report_t rqm_verify_theorems_9_and_10(const rqm_mat_t *rho0, const rqm_mat_t *X, const rqm_mat_t *H);

#endif
