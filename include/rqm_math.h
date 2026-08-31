#ifndef RQM_MATH_H
#define RQM_MATH_H

#include "rqm_common.h"

rqm_mat_t rqm_mat_zero(int dim);
rqm_mat_t rqm_mat_identity(int dim);
rqm_mat_t rqm_mat_pauli_x(void);
rqm_mat_t rqm_mat_pauli_y(void);
rqm_mat_t rqm_mat_pauli_z(void);
rqm_mat_t rqm_mat_annihilation(int dim);
rqm_mat_t rqm_mat_creation(int dim);
rqm_mat_t rqm_mat_position(int dim);
rqm_mat_t rqm_mat_momentum(int dim);
rqm_mat_t rqm_mat_number(int dim);

rqm_vec_t rqm_vec_zero(int dim);
rqm_vec_t rqm_vec_basis(int dim, int index);
rqm_vec_t rqm_vec_qubit(double theta, double phi);
double rqm_vec_norm(const rqm_vec_t *v);
rqm_vec_t rqm_vec_normalize(const rqm_vec_t *v);
rqm_cdouble rqm_vec_dot(const rqm_vec_t *u, const rqm_vec_t *v);
rqm_mat_t rqm_vec_outer(const rqm_vec_t *u, const rqm_vec_t *v);

rqm_mat_t rqm_mat_add(const rqm_mat_t *A, const rqm_mat_t *B);
rqm_mat_t rqm_mat_sub(const rqm_mat_t *A, const rqm_mat_t *B);
rqm_mat_t rqm_mat_mul(const rqm_mat_t *A, const rqm_mat_t *B);
rqm_mat_t rqm_mat_scale(const rqm_mat_t *A, rqm_cdouble s);
rqm_mat_t rqm_mat_dagger(const rqm_mat_t *A);
rqm_cdouble rqm_mat_trace(const rqm_mat_t *A);
rqm_mat_t rqm_mat_commutator(const rqm_mat_t *A, const rqm_mat_t *B);
rqm_mat_t rqm_mat_anticommutator(const rqm_mat_t *A, const rqm_mat_t *B);
rqm_cdouble rqm_mat_hs_dot(const rqm_mat_t *A, const rqm_mat_t *B);
double rqm_mat_hs_norm(const rqm_mat_t *A);
double rqm_mat_hs_dist(const rqm_mat_t *A, const rqm_mat_t *B);
bool rqm_mat_is_hermitian(const rqm_mat_t *A, double tol);
bool rqm_mat_is_density(const rqm_mat_t *rho, double tol);

rqm_status_t rqm_hermitian_eigensolve(const rqm_mat_t *A, rqm_eigen_t *out_eigen);
rqm_mat_t rqm_mat_exp_hermitian(const rqm_mat_t *H, rqm_cdouble scale);
rqm_mat_t rqm_mat_log_density(const rqm_mat_t *rho, double zero_cutoff);

double rqm_von_neumann_entropy(const rqm_mat_t *rho);
double rqm_von_neumann_entropy_base2(const rqm_mat_t *rho);
double rqm_relative_entropy(const rqm_mat_t *rho, const rqm_mat_t *sigma);
double rqm_purity(const rqm_mat_t *rho);
double rqm_fidelity_pure(const rqm_vec_t *psi, const rqm_vec_t *phi);

typedef rqm_mat_t (*rqm_lindblad_rhs_fn)(const rqm_mat_t *rho, void *user_data);
rqm_mat_t rqm_rk4_step_mat(rqm_lindblad_rhs_fn rhs, const rqm_mat_t *rho, double dt, void *user_data);

typedef void (*rqm_ode_rhs_fn)(double t, const double *y, double *dydt, void *user_data);
void rqm_rk4_step_ode(rqm_ode_rhs_fn rhs, double t, const double *y, double *y_next, int dim, double dt, void *user_data);

#endif
