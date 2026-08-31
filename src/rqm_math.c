#include "rqm_math.h"

rqm_mat_t rqm_mat_zero(int dim)
{
    rqm_mat_t m;
    m.dim = dim;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            m.data[i][j] = 0.0 + 0.0 * I;
        }
    }
    return m;
}

rqm_mat_t rqm_mat_identity(int dim)
{
    rqm_mat_t m = rqm_mat_zero(dim);
    for (int i = 0; i < dim; i++) {
        m.data[i][i] = 1.0 + 0.0 * I;
    }
    return m;
}

rqm_mat_t rqm_mat_pauli_x(void)
{
    rqm_mat_t m = rqm_mat_zero(2);
    m.data[0][1] = 1.0 + 0.0 * I;
    m.data[1][0] = 1.0 + 0.0 * I;
    return m;
}

rqm_mat_t rqm_mat_pauli_y(void)
{
    rqm_mat_t m = rqm_mat_zero(2);
    m.data[0][1] = 0.0 - 1.0 * I;
    m.data[1][0] = 0.0 + 1.0 * I;
    return m;
}

rqm_mat_t rqm_mat_pauli_z(void)
{
    rqm_mat_t m = rqm_mat_zero(2);
    m.data[0][0] = 1.0 + 0.0 * I;
    m.data[1][1] = -1.0 + 0.0 * I;
    return m;
}

rqm_mat_t rqm_mat_annihilation(int dim)
{
    rqm_mat_t a = rqm_mat_zero(dim);
    for (int i = 0; i < dim - 1; i++) {
        a.data[i][i + 1] = sqrt((double)(i + 1)) + 0.0 * I;
    }
    return a;
}

rqm_mat_t rqm_mat_creation(int dim)
{
    rqm_mat_t a = rqm_mat_annihilation(dim);
    return rqm_mat_dagger(&a);
}

rqm_mat_t rqm_mat_position(int dim)
{
    rqm_mat_t a = rqm_mat_annihilation(dim);
    rqm_mat_t adag = rqm_mat_creation(dim);
    rqm_mat_t sum = rqm_mat_add(&a, &adag);
    return rqm_mat_scale(&sum, (1.0 / sqrt(2.0)) + 0.0 * I);
}

rqm_mat_t rqm_mat_momentum(int dim)
{
    rqm_mat_t a = rqm_mat_annihilation(dim);
    rqm_mat_t adag = rqm_mat_creation(dim);
    rqm_mat_t diff = rqm_mat_sub(&a, &adag);
    return rqm_mat_scale(&diff, (0.0 + (1.0 / (I * sqrt(2.0)))) );
}

rqm_mat_t rqm_mat_number(int dim)
{
    rqm_mat_t n = rqm_mat_zero(dim);
    for (int i = 0; i < dim; i++) {
        n.data[i][i] = (double)i + 0.0 * I;
    }
    return n;
}

rqm_vec_t rqm_vec_zero(int dim)
{
    rqm_vec_t v;
    v.dim = dim;
    for (int i = 0; i < dim; i++) {
        v.data[i] = 0.0 + 0.0 * I;
    }
    return v;
}

rqm_vec_t rqm_vec_basis(int dim, int index)
{
    rqm_vec_t v = rqm_vec_zero(dim);
    if (index >= 0 && index < dim) {
        v.data[index] = 1.0 + 0.0 * I;
    }
    return v;
}

rqm_vec_t rqm_vec_qubit(double theta, double phi)
{
    rqm_vec_t v = rqm_vec_zero(2);
    v.data[0] = cos(theta / 2.0) + 0.0 * I;
    v.data[1] = sin(theta / 2.0) * (cos(phi) + sin(phi) * I);
    return v;
}

double rqm_vec_norm(const rqm_vec_t *v)
{
    double sum = 0.0;
    for (int i = 0; i < v->dim; i++) {
        sum += creal(v->data[i]) * creal(v->data[i]) + cimag(v->data[i]) * cimag(v->data[i]);
    }
    return sqrt(sum);
}

rqm_vec_t rqm_vec_normalize(const rqm_vec_t *v)
{
    double n = rqm_vec_norm(v);
    rqm_vec_t res = rqm_vec_zero(v->dim);
    if (n > RQM_EPSILON) {
        for (int i = 0; i < v->dim; i++) {
            res.data[i] = v->data[i] / n;
        }
    }
    return res;
}

rqm_cdouble rqm_vec_dot(const rqm_vec_t *u, const rqm_vec_t *v)
{
    rqm_cdouble sum = 0.0 + 0.0 * I;
    int n = u->dim < v->dim ? u->dim : v->dim;
    for (int i = 0; i < n; i++) {
        sum += conj(u->data[i]) * v->data[i];
    }
    return sum;
}

rqm_mat_t rqm_vec_outer(const rqm_vec_t *u, const rqm_vec_t *v)
{
    int dim = u->dim;
    rqm_mat_t m = rqm_mat_zero(dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            m.data[i][j] = u->data[i] * conj(v->data[j]);
        }
    }
    return m;
}

rqm_mat_t rqm_mat_add(const rqm_mat_t *A, const rqm_mat_t *B)
{
    rqm_mat_t res = rqm_mat_zero(A->dim);
    for (int i = 0; i < A->dim; i++) {
        for (int j = 0; j < A->dim; j++) {
            res.data[i][j] = A->data[i][j] + B->data[i][j];
        }
    }
    return res;
}

rqm_mat_t rqm_mat_sub(const rqm_mat_t *A, const rqm_mat_t *B)
{
    rqm_mat_t res = rqm_mat_zero(A->dim);
    for (int i = 0; i < A->dim; i++) {
        for (int j = 0; j < A->dim; j++) {
            res.data[i][j] = A->data[i][j] - B->data[i][j];
        }
    }
    return res;
}

rqm_mat_t rqm_mat_mul(const rqm_mat_t *A, const rqm_mat_t *B)
{
    rqm_mat_t res = rqm_mat_zero(A->dim);
    for (int i = 0; i < A->dim; i++) {
        for (int j = 0; j < A->dim; j++) {
            rqm_cdouble sum = 0.0 + 0.0 * I;
            for (int k = 0; k < A->dim; k++) {
                sum += A->data[i][k] * B->data[k][j];
            }
            res.data[i][j] = sum;
        }
    }
    return res;
}

rqm_mat_t rqm_mat_scale(const rqm_mat_t *A, rqm_cdouble s)
{
    rqm_mat_t res = rqm_mat_zero(A->dim);
    for (int i = 0; i < A->dim; i++) {
        for (int j = 0; j < A->dim; j++) {
            res.data[i][j] = A->data[i][j] * s;
        }
    }
    return res;
}

rqm_mat_t rqm_mat_dagger(const rqm_mat_t *A)
{
    rqm_mat_t res = rqm_mat_zero(A->dim);
    for (int i = 0; i < A->dim; i++) {
        for (int j = 0; j < A->dim; j++) {
            res.data[i][j] = conj(A->data[j][i]);
        }
    }
    return res;
}

rqm_cdouble rqm_mat_trace(const rqm_mat_t *A)
{
    rqm_cdouble tr = 0.0 + 0.0 * I;
    for (int i = 0; i < A->dim; i++) {
        tr += A->data[i][i];
    }
    return tr;
}

rqm_mat_t rqm_mat_commutator(const rqm_mat_t *A, const rqm_mat_t *B)
{
    rqm_mat_t AB = rqm_mat_mul(A, B);
    rqm_mat_t BA = rqm_mat_mul(B, A);
    return rqm_mat_sub(&AB, &BA);
}

rqm_mat_t rqm_mat_anticommutator(const rqm_mat_t *A, const rqm_mat_t *B)
{
    rqm_mat_t AB = rqm_mat_mul(A, B);
    rqm_mat_t BA = rqm_mat_mul(B, A);
    return rqm_mat_add(&AB, &BA);
}

rqm_cdouble rqm_mat_hs_dot(const rqm_mat_t *A, const rqm_mat_t *B)
{
    rqm_mat_t Adag = rqm_mat_dagger(A);
    rqm_mat_t prod = rqm_mat_mul(&Adag, B);
    return rqm_mat_trace(&prod);
}

double rqm_mat_hs_norm(const rqm_mat_t *A)
{
    rqm_cdouble dot = rqm_mat_hs_dot(A, A);
    return sqrt(fabs(creal(dot)));
}

double rqm_mat_hs_dist(const rqm_mat_t *A, const rqm_mat_t *B)
{
    rqm_mat_t diff = rqm_mat_sub(A, B);
    return rqm_mat_hs_norm(&diff);
}

bool rqm_mat_is_hermitian(const rqm_mat_t *A, double tol)
{
    for (int i = 0; i < A->dim; i++) {
        for (int j = 0; j < A->dim; j++) {
            if (cabs(A->data[i][j] - conj(A->data[j][i])) > tol) {
                return false;
            }
        }
    }
    return true;
}

bool rqm_mat_is_density(const rqm_mat_t *rho, double tol)
{
    if (!rqm_mat_is_hermitian(rho, tol)) {
        return false;
    }
    rqm_cdouble tr = rqm_mat_trace(rho);
    if (fabs(creal(tr) - 1.0) > tol || fabs(cimag(tr)) > tol) {
        return false;
    }
    rqm_eigen_t eigen;
    if (rqm_hermitian_eigensolve(rho, &eigen) != RQM_SUCCESS) {
        return false;
    }
    for (int i = 0; i < rho->dim; i++) {
        if (eigen.values[i] < -tol) {
            return false;
        }
    }
    return true;
}

rqm_status_t rqm_hermitian_eigensolve(const rqm_mat_t *A, rqm_eigen_t *out_eigen)
{
    int n = A->dim;
    out_eigen->dim = n;
    rqm_mat_t H = *A;
    rqm_mat_t V = rqm_mat_identity(n);
    int max_iter = 150;
    double tol = 1e-12;

    for (int iter = 0; iter < max_iter; iter++) {
        double max_off = 0.0;
        int p = 0;
        int q = 1;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double val = cabs(H.data[i][j]);
                if (val > max_off) {
                    max_off = val;
                    p = i;
                    q = j;
                }
            }
        }

        if (max_off < tol) {
            break;
        }

        double hpp = creal(H.data[p][p]);
        double hqq = creal(H.data[q][q]);
        rqm_cdouble hpq = H.data[p][q];
        double abs_hpq = cabs(hpq);

        if (abs_hpq < 1e-15) {
            continue;
        }

        rqm_cdouble phase = hpq / abs_hpq;
        double tau = (hqq - hpp) / (2.0 * abs_hpq);
        double t = 0.0;
        if (tau >= 0.0) {
            t = 1.0 / (tau + sqrt(1.0 + tau * tau));
        } else {
            t = -1.0 / (-tau + sqrt(1.0 + tau * tau));
        }

        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;
        rqm_cdouble s_phase = s * phase;

        for (int i = 0; i < n; i++) {
            if (i != p && i != q) {
                rqm_cdouble hip = H.data[i][p];
                rqm_cdouble hiq = H.data[i][q];
                H.data[i][p] = c * hip - conj(s_phase) * hiq;
                H.data[p][i] = conj(H.data[i][p]);
                H.data[i][q] = s_phase * hip + c * hiq;
                H.data[q][i] = conj(H.data[i][q]);
            }
        }

        H.data[p][p] = (hpp - t * abs_hpq) + 0.0 * I;
        H.data[q][q] = (hqq + t * abs_hpq) + 0.0 * I;
        H.data[p][q] = 0.0 + 0.0 * I;
        H.data[q][p] = 0.0 + 0.0 * I;

        for (int i = 0; i < n; i++) {
            rqm_cdouble vip = V.data[i][p];
            rqm_cdouble viq = V.data[i][q];
            V.data[i][p] = c * vip - conj(s_phase) * viq;
            V.data[i][q] = s_phase * vip + c * viq;
        }
    }

    for (int i = 0; i < n; i++) {
        out_eigen->values[i] = creal(H.data[i][i]);
        out_eigen->vectors[i].dim = n;
        for (int j = 0; j < n; j++) {
            out_eigen->vectors[i].data[j] = V.data[j][i];
        }
    }

    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (out_eigen->values[j] > out_eigen->values[i]) {
                double tmp_val = out_eigen->values[i];
                out_eigen->values[i] = out_eigen->values[j];
                out_eigen->values[j] = tmp_val;
                rqm_vec_t tmp_vec = out_eigen->vectors[i];
                out_eigen->vectors[i] = out_eigen->vectors[j];
                out_eigen->vectors[j] = tmp_vec;
            }
        }
    }

    return RQM_SUCCESS;
}

rqm_mat_t rqm_mat_exp_hermitian(const rqm_mat_t *H, rqm_cdouble scale)
{
    rqm_eigen_t eigen;
    rqm_hermitian_eigensolve(H, &eigen);
    rqm_mat_t res = rqm_mat_zero(H->dim);
    for (int k = 0; k < H->dim; k++) {
        rqm_cdouble factor = cexp(scale * (eigen.values[k] + 0.0 * I));
        rqm_mat_t proj = rqm_vec_outer(&eigen.vectors[k], &eigen.vectors[k]);
        rqm_mat_t term = rqm_mat_scale(&proj, factor);
        res = rqm_mat_add(&res, &term);
    }
    return res;
}

rqm_mat_t rqm_mat_log_density(const rqm_mat_t *rho, double zero_cutoff)
{
    rqm_eigen_t eigen;
    rqm_hermitian_eigensolve(rho, &eigen);
    rqm_mat_t res = rqm_mat_zero(rho->dim);
    for (int k = 0; k < rho->dim; k++) {
        if (eigen.values[k] > zero_cutoff) {
            double log_val = log(eigen.values[k]);
            rqm_mat_t proj = rqm_vec_outer(&eigen.vectors[k], &eigen.vectors[k]);
            rqm_mat_t term = rqm_mat_scale(&proj, log_val + 0.0 * I);
            res = rqm_mat_add(&res, &term);
        }
    }
    return res;
}

double rqm_von_neumann_entropy(const rqm_mat_t *rho)
{
    rqm_eigen_t eigen;
    rqm_hermitian_eigensolve(rho, &eigen);
    double s = 0.0;
    for (int k = 0; k < rho->dim; k++) {
        if (eigen.values[k] > 1e-14) {
            s -= eigen.values[k] * log(eigen.values[k]);
        }
    }
    return s < 0.0 ? 0.0 : s;
}

double rqm_von_neumann_entropy_base2(const rqm_mat_t *rho)
{
    return rqm_von_neumann_entropy(rho) / log(2.0);
}

double rqm_relative_entropy(const rqm_mat_t *rho, const rqm_mat_t *sigma)
{
    rqm_mat_t log_rho = rqm_mat_log_density(rho, 1e-14);
    rqm_mat_t log_sigma = rqm_mat_log_density(sigma, 1e-14);
    rqm_mat_t diff_log = rqm_mat_sub(&log_rho, &log_sigma);
    rqm_mat_t prod = rqm_mat_mul(rho, &diff_log);
    rqm_cdouble tr = rqm_mat_trace(&prod);
    double s = creal(tr);
    return s < 0.0 ? 0.0 : s;
}

double rqm_purity(const rqm_mat_t *rho)
{
    rqm_mat_t sq = rqm_mat_mul(rho, rho);
    return creal(rqm_mat_trace(&sq));
}

double rqm_fidelity_pure(const rqm_vec_t *psi, const rqm_vec_t *phi)
{
    rqm_cdouble dot = rqm_vec_dot(psi, phi);
    return creal(dot) * creal(dot) + cimag(dot) * cimag(dot);
}

rqm_mat_t rqm_rk4_step_mat(rqm_lindblad_rhs_fn rhs, const rqm_mat_t *rho, double dt, void *user_data)
{
    rqm_mat_t k1 = rhs(rho, user_data);
    rqm_mat_t s1 = rqm_mat_scale(&k1, (dt * 0.5) + 0.0 * I);
    rqm_mat_t rho_k1 = rqm_mat_add(rho, &s1);
    rqm_mat_t k2 = rhs(&rho_k1, user_data);
    rqm_mat_t s2 = rqm_mat_scale(&k2, (dt * 0.5) + 0.0 * I);
    rqm_mat_t rho_k2 = rqm_mat_add(rho, &s2);
    rqm_mat_t k3 = rhs(&rho_k2, user_data);
    rqm_mat_t s3 = rqm_mat_scale(&k3, dt + 0.0 * I);
    rqm_mat_t rho_k3 = rqm_mat_add(rho, &s3);
    rqm_mat_t k4 = rhs(&rho_k3, user_data);

    rqm_mat_t k2_2 = rqm_mat_scale(&k2, 2.0 + 0.0 * I);
    rqm_mat_t k3_2 = rqm_mat_scale(&k3, 2.0 + 0.0 * I);
    rqm_mat_t sum = rqm_mat_add(&k1, &k2_2);
    sum = rqm_mat_add(&sum, &k3_2);
    sum = rqm_mat_add(&sum, &k4);

    rqm_mat_t delta = rqm_mat_scale(&sum, (dt / 6.0) + 0.0 * I);
    return rqm_mat_add(rho, &delta);
}

void rqm_rk4_step_ode(rqm_ode_rhs_fn rhs, double t, const double *y, double *y_next, int dim, double dt, void *user_data)
{
    double k1[16];
    double k2[16];
    double k3[16];
    double k4[16];
    double y_tmp[16];

    rhs(t, y, k1, user_data);
    for (int i = 0; i < dim; i++) {
        y_tmp[i] = y[i] + 0.5 * dt * k1[i];
    }

    rhs(t + 0.5 * dt, y_tmp, k2, user_data);
    for (int i = 0; i < dim; i++) {
        y_tmp[i] = y[i] + 0.5 * dt * k2[i];
    }

    rhs(t + 0.5 * dt, y_tmp, k3, user_data);
    for (int i = 0; i < dim; i++) {
        y_tmp[i] = y[i] + dt * k3[i];
    }

    rhs(t + dt, y_tmp, k4, user_data);
    for (int i = 0; i < dim; i++) {
        y_next[i] = y[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}
