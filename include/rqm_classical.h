#ifndef RQM_CLASSICAL_H
#define RQM_CLASSICAL_H

#include "rqm_common.h"

typedef struct {
    double mass;
    double spring_k;
    double lambda_anharm;
} rqm_classical_params_t;

typedef struct {
    double q;
    double p;
} rqm_phase_point_t;

typedef struct {
    int num_points;
    rqm_phase_point_t points[256];
    double weights[256];
} rqm_classical_ensemble_t;

typedef struct {
    bool hm_1d_valid;
    bool di_symp_valid;
    bool di_poi_valid;
    bool dr_div_valid;
    bool dr_den_valid;
    bool dr_ev_valid;
    bool dr_ther_valid;
    bool dr_info_valid;
    double initial_energy;
    double final_energy;
    double initial_entropy;
    double final_entropy;
    double max_divergence;
    double symplectic_det;
} rqm_classical_thm2_result_t;

double rqm_classical_hamiltonian(double q, double p, const rqm_classical_params_t *params);
void rqm_classical_flow(double q, double p, double *dq_dt, double *dp_dt, const rqm_classical_params_t *params);
double rqm_classical_poisson_bracket(double q, double p, double (*f)(double, double), double (*g)(double, double));
double rqm_classical_divergence(double q, double p, const rqm_classical_params_t *params);
rqm_phase_point_t rqm_classical_step(rqm_phase_point_t state, double dt, const rqm_classical_params_t *params);
rqm_classical_ensemble_t rqm_classical_create_gaussian_ensemble(double q0, double p0, double sigma_q, double sigma_p, int num_points);
double rqm_classical_ensemble_entropy(const rqm_classical_ensemble_t *ens);
rqm_classical_thm2_result_t rqm_classical_verify_theorem2(const rqm_classical_params_t *params, double duration, double dt);

#endif
