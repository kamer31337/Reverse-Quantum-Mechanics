#ifndef RQM_TRIPLE_EQUILIBRIA_H
#define RQM_TRIPLE_EQUILIBRIA_H

#include "rqm_common.h"

typedef struct {
    rqm_mat_t H_dynamic;
    double commutator_norm;
    bool is_dynamic_equilibrium;
} rqm_dynamic_equilibrium_t;

typedef struct {
    rqm_mat_t X_observable;
    rqm_mat_t projection_output;
    double measurement_diff_norm;
    bool is_spectral_equilibrium;
} rqm_spectral_equilibrium_t;

typedef struct {
    double beta;
    bool is_full_rank;
    rqm_mat_t H_thermodynamic;
    rqm_mat_t exact_or_approx_gibbs;
    double gibbs_error_norm;
    int num_k_steps;
    double k_penalties[16];
    double k_errors[16];
    bool is_thermodynamic_equilibrium;
} rqm_thermodynamic_equilibrium_t;

typedef struct {
    rqm_mat_t target_state;
    rqm_dynamic_equilibrium_t dynamic_eq;
    rqm_spectral_equilibrium_t spectral_eq;
    rqm_thermodynamic_equilibrium_t thermo_eq;
    bool theorem12_all_equilibria_verified;
} rqm_triple_equilibria_result_t;

rqm_dynamic_equilibrium_t rqm_verify_dynamic_equilibrium(const rqm_mat_t *rho);
rqm_spectral_equilibrium_t rqm_verify_spectral_equilibrium(const rqm_mat_t *rho);
rqm_thermodynamic_equilibrium_t rqm_verify_thermodynamic_equilibrium(const rqm_mat_t *rho, double beta);
rqm_triple_equilibria_result_t rqm_verify_theorem12(const rqm_mat_t *rho, double beta);

#endif
