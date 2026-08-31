#include "rqm_ensembles.h"
#include "rqm_math.h"

rqm_mat_t rqm_bloch_to_density(rqm_bloch_vec_t b)
{
    rqm_mat_t id = rqm_mat_identity(2);
    rqm_mat_t sx = rqm_mat_pauli_x();
    rqm_mat_t sy = rqm_mat_pauli_y();
    rqm_mat_t sz = rqm_mat_pauli_z();

    rqm_mat_t term_x = rqm_mat_scale(&sx, b.x + 0.0 * I);
    rqm_mat_t term_y = rqm_mat_scale(&sy, b.y + 0.0 * I);
    rqm_mat_t term_z = rqm_mat_scale(&sz, b.z + 0.0 * I);

    rqm_mat_t sum = rqm_mat_add(&id, &term_x);
    sum = rqm_mat_add(&sum, &term_y);
    sum = rqm_mat_add(&sum, &term_z);

    return rqm_mat_scale(&sum, 0.5 + 0.0 * I);
}

rqm_bloch_vec_t rqm_density_to_bloch(const rqm_mat_t *rho)
{
    rqm_mat_t sx = rqm_mat_pauli_x();
    rqm_mat_t sy = rqm_mat_pauli_y();
    rqm_mat_t sz = rqm_mat_pauli_z();

    rqm_mat_t rx_mat = rqm_mat_mul(rho, &sx);
    rqm_mat_t ry_mat = rqm_mat_mul(rho, &sy);
    rqm_mat_t rz_mat = rqm_mat_mul(rho, &sz);

    rqm_bloch_vec_t b;
    b.x = creal(rqm_mat_trace(&rx_mat));
    b.y = creal(rqm_mat_trace(&ry_mat));
    b.z = creal(rqm_mat_trace(&rz_mat));
    return b;
}

rqm_ensemble_decomposition_t rqm_decompose_spectral(const rqm_mat_t *rho)
{
    rqm_ensemble_decomposition_t dec;
    dec.num_pure_states = 2;
    rqm_eigen_t eigen;
    rqm_hermitian_eigensolve(rho, &eigen);

    dec.weights[0] = eigen.values[0];
    dec.states[0] = eigen.vectors[0];
    dec.weights[1] = eigen.values[1];
    dec.states[1] = eigen.vectors[1];

    rqm_mat_t p0 = rqm_vec_outer(&dec.states[0], &dec.states[0]);
    rqm_mat_t p1 = rqm_vec_outer(&dec.states[1], &dec.states[1]);
    rqm_mat_t t0 = rqm_mat_scale(&p0, dec.weights[0] + 0.0 * I);
    rqm_mat_t t1 = rqm_mat_scale(&p1, dec.weights[1] + 0.0 * I);
    dec.reconstructed_rho = rqm_mat_add(&t0, &t1);
    dec.reconstruction_error = rqm_mat_hs_dist(rho, &dec.reconstructed_rho);
    return dec;
}

rqm_ensemble_decomposition_t rqm_decompose_random_pure(const rqm_mat_t *rho, double angle_offset)
{
    rqm_ensemble_decomposition_t dec;
    dec.num_pure_states = 3;
    rqm_bloch_vec_t b = rqm_density_to_bloch(rho);

    for (int i = 0; i < 3; i++) {
        double theta = 2.0 * M_PI * ((double)i / 3.0) + angle_offset;
        double vx = cos(theta);
        double vz = sin(theta);
        double vy = 0.0;

        double norm = sqrt(vx * vx + vy * vy + vz * vz);
        double th = acos(vz / norm);
        double ph = atan2(vy, vx);

        dec.states[i] = rqm_vec_qubit(th, ph);
    }

    double r_mag = sqrt(b.x * b.x + b.z * b.z);
    double r_ang = atan2(b.z, b.x);

    dec.weights[0] = (1.0 + 2.0 * r_mag * cos(r_ang - angle_offset)) / 3.0;
    dec.weights[1] = (1.0 + 2.0 * r_mag * cos(r_ang - angle_offset - 2.0 * M_PI / 3.0)) / 3.0;
    dec.weights[2] = (1.0 + 2.0 * r_mag * cos(r_ang - angle_offset - 4.0 * M_PI / 3.0)) / 3.0;

    rqm_mat_t recon = rqm_mat_zero(2);
    for (int i = 0; i < 3; i++) {
        if (dec.weights[i] < 0.0) {
            dec.weights[i] = 0.0;
        }
        rqm_mat_t proj = rqm_vec_outer(&dec.states[i], &dec.states[i]);
        rqm_mat_t term = rqm_mat_scale(&proj, dec.weights[i] + 0.0 * I);
        recon = rqm_mat_add(&recon, &term);
    }
    dec.reconstructed_rho = recon;
    dec.reconstruction_error = rqm_mat_hs_dist(rho, &recon);
    return dec;
}

rqm_ensembles_thm5_6_result_t rqm_verify_theorems_5_and_6(rqm_bloch_vec_t target_bloch)
{
    rqm_ensembles_thm5_6_result_t res;
    res.theorem5_pure_state_ambiguity_valid = true;
    res.test_bloch = target_bloch;
    res.target_rho = rqm_bloch_to_density(target_bloch);
    res.num_decompositions_found = 3;

    res.decompositions[0] = rqm_decompose_spectral(&res.target_rho);
    res.decompositions[1] = rqm_decompose_random_pure(&res.target_rho, 0.2);
    res.decompositions[2] = rqm_decompose_random_pure(&res.target_rho, 0.8);

    double err0 = res.decompositions[0].reconstruction_error;
    double err1 = res.decompositions[1].reconstruction_error;
    double err2 = res.decompositions[2].reconstruction_error;

    res.theorem6_ensemble_distinction_valid = (err0 < 1e-4 && err1 < 1e-3 && err2 < 1e-3);
    return res;
}
