#ifndef RQM_ENSEMBLES_H
#define RQM_ENSEMBLES_H

#include "rqm_common.h"

typedef struct {
    double x;
    double y;
    double z;
} rqm_bloch_vec_t;

typedef struct {
    int num_pure_states;
    double weights[8];
    rqm_vec_t states[8];
    rqm_mat_t reconstructed_rho;
    double reconstruction_error;
} rqm_ensemble_decomposition_t;

typedef struct {
    bool theorem5_pure_state_ambiguity_valid;
    bool theorem6_ensemble_distinction_valid;
    rqm_bloch_vec_t test_bloch;
    rqm_mat_t target_rho;
    int num_decompositions_found;
    rqm_ensemble_decomposition_t decompositions[4];
} rqm_ensembles_thm5_6_result_t;

rqm_mat_t rqm_bloch_to_density(rqm_bloch_vec_t b);
rqm_bloch_vec_t rqm_density_to_bloch(const rqm_mat_t *rho);
rqm_ensemble_decomposition_t rqm_decompose_spectral(const rqm_mat_t *rho);
rqm_ensemble_decomposition_t rqm_decompose_random_pure(const rqm_mat_t *rho, double angle_offset);
rqm_ensembles_thm5_6_result_t rqm_verify_theorems_5_and_6(rqm_bloch_vec_t target_bloch);

#endif
