#include "rqm_tests.h"
#include "rqm_math.h"
#include "rqm_classical.h"
#include "rqm_hilbert_pathology.h"
#include "rqm_ensembles.h"
#include "rqm_born_entropy.h"
#include "rqm_unitary.h"
#include "rqm_equilibration.h"
#include "rqm_classical_limit.h"
#include "rqm_triple_equilibria.h"
#include "rqm_ui.h"

bool rqm_test_math_suite(void)
{
    rqm_mat_t sx = rqm_mat_pauli_x();
    rqm_mat_t sy = rqm_mat_pauli_y();
    rqm_mat_t sz = rqm_mat_pauli_z();

    rqm_mat_t comm_xy = rqm_mat_commutator(&sx, &sy);
    rqm_mat_t two_i_sz = rqm_mat_scale(&sz, 2.0 * I);
    double comm_diff = rqm_mat_hs_dist(&comm_xy, &two_i_sz);
    if (comm_diff > 1e-10) return false;

    rqm_eigen_t eigen_z;
    if (rqm_hermitian_eigensolve(&sz, &eigen_z) != RQM_SUCCESS) return false;
    if (fabs(eigen_z.values[0] - 1.0) > 1e-10 || fabs(eigen_z.values[1] - (-1.0)) > 1e-10) return false;

    rqm_vec_t psi = rqm_vec_qubit(M_PI / 4.0, M_PI / 3.0);
    rqm_mat_t rho = rqm_vec_outer(&psi, &psi);
    double pur = rqm_purity(&rho);
    double s = rqm_von_neumann_entropy_base2(&rho);
    if (fabs(pur - 1.0) > 1e-10 || fabs(s - 0.0) > 1e-10) return false;

    return true;
}

bool rqm_test_theorem_2(void)
{
    rqm_classical_params_t params = {1.0, 1.0, 0.1};
    rqm_classical_thm2_result_t res = rqm_classical_verify_theorem2(&params, 5.0, 0.01);
    return res.hm_1d_valid && res.di_symp_valid && res.dr_div_valid && res.dr_ev_valid && res.dr_info_valid;
}

bool rqm_test_theorems_3_and_4(void)
{
    rqm_pathology_report_t rep = rqm_verify_theorems_3_and_4();
    return rep.theorem3_inconsistency_demonstrated && rep.theorem4_qsto_v_resolved;
}

bool rqm_test_theorems_5_and_6(void)
{
    rqm_bloch_vec_t target = {0.3, 0.0, 0.4};
    rqm_ensembles_thm5_6_result_t res = rqm_verify_theorems_5_and_6(target);
    return res.theorem5_pure_state_ambiguity_valid && res.theorem6_ensemble_distinction_valid;
}

bool rqm_test_theorem_7(void)
{
    rqm_born_thm7_result_t res = rqm_verify_theorem7(16);
    return res.is_strictly_monotonic && res.orthogonal_gives_one_bit && res.identical_gives_zero_bit && (res.max_error < 1e-6);
}

bool rqm_test_theorem_8(void)
{
    rqm_mat_t H = rqm_mat_pauli_x();
    rqm_vec_t psi0 = rqm_vec_qubit(0.0, 0.0);
    rqm_vec_t phi0 = rqm_vec_qubit(M_PI / 2.0, 0.0);
    rqm_mat_t rho0 = rqm_mat_identity(2);
    rho0 = rqm_mat_scale(&rho0, 0.5 + 0.0 * I);

    rqm_unitary_thm8_result_t res = rqm_verify_theorem8(&H, &rho0, &psi0, &phi0, 4.0, 32);
    return res.dr_sceq_valid && res.dr_unit_valid && res.dr_ev_valid && res.dr_prob_valid && res.dr_info_valid;
}

bool rqm_test_theorems_9_and_10(void)
{
    rqm_vec_t psi = rqm_vec_qubit(M_PI / 3.0, M_PI / 4.0);
    rqm_mat_t rho0 = rqm_vec_outer(&psi, &psi);
    rqm_mat_t X = rqm_mat_pauli_z();
    rqm_mat_t H = rqm_mat_pauli_x();

    rqm_equilibration_report_t rep = rqm_verify_theorems_9_and_10(&rho0, &X, &H);
    return rep.theorem9_valid && rep.theorem10_valid;
}

bool rqm_test_theorem_11(void)
{
    rqm_classical_limit_result_t res = rqm_verify_theorem11(1, 1.0, 3.0, 20);
    return res.theorem11_valid;
}

bool rqm_test_theorem_12(void)
{
    rqm_vec_t psi = rqm_vec_qubit(M_PI / 3.0, 0.0);
    rqm_mat_t rho_pure = rqm_vec_outer(&psi, &psi);
    rqm_triple_equilibria_result_t res_pure = rqm_verify_theorem12(&rho_pure, 1.0);
    if (!res_pure.theorem12_all_equilibria_verified) return false;

    rqm_bloch_vec_t b = {0.2, 0.1, 0.5};
    rqm_mat_t rho_mixed = rqm_bloch_to_density(b);
    rqm_triple_equilibria_result_t res_mixed = rqm_verify_theorem12(&rho_mixed, 1.0);
    return res_mixed.theorem12_all_equilibria_verified;
}

rqm_test_summary_t rqm_test_run_all(void)
{
    rqm_test_summary_t sum;
    sum.total_tests = 9;
    sum.passed_tests = 0;
    sum.failed_tests = 0;
    sum.max_discrepancy = 0.0;

    printf(ANSI_BOLD ANSI_CYAN "\n================================================================================\n" ANSI_RESET);
    printf(ANSI_BOLD ANSI_WHITE "               REVERSE QUANTUM MECHANICS TEST SUITE VERIFICATION\n" ANSI_RESET);
    printf(ANSI_BOLD ANSI_CYAN "================================================================================\n\n" ANSI_RESET);

    bool t_math = rqm_test_math_suite();
    printf("  [Test 1/9] Complex Linear Algebra & Eigensolver:         %s\n", t_math ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t_math) sum.passed_tests++; else sum.failed_tests++;

    bool t2 = rqm_test_theorem_2();
    printf("  [Test 2/9] Theorem 2 (Classical 1-DOF Equivalence):      %s\n", t2 ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t2) sum.passed_tests++; else sum.failed_tests++;

    bool t3_4 = rqm_test_theorems_3_and_4();
    printf("  [Test 3/9] Theorems 3 & 4 (Hilbert Space Pathology/QSTO): %s\n", t3_4 ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t3_4) sum.passed_tests++; else sum.failed_tests++;

    bool t5_6 = rqm_test_theorems_5_and_6();
    printf("  [Test 4/9] Theorems 5 & 6 (Ensemble Convex Structure):   %s\n", t5_6 ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t5_6) sum.passed_tests++; else sum.failed_tests++;

    bool t7 = rqm_test_theorem_7();
    printf("  [Test 5/9] Theorem 7 (Born Rule <-> Entropy Bijection):  %s\n", t7 ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t7) sum.passed_tests++; else sum.failed_tests++;

    bool t8 = rqm_test_theorem_8();
    printf("  [Test 6/9] Theorem 8 (Unitary Evolution Equivalence):    %s\n", t8 ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t8) sum.passed_tests++; else sum.failed_tests++;

    bool t9_10 = rqm_test_theorems_9_and_10();
    printf("  [Test 7/9] Theorems 9 & 10 (Lindblad & Zeno Equilibr):   %s\n", t9_10 ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t9_10) sum.passed_tests++; else sum.failed_tests++;

    bool t11 = rqm_test_theorem_11();
    printf("  [Test 8/9] Theorem 11 (Classical Limit High Entropy):    %s\n", t11 ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t11) sum.passed_tests++; else sum.failed_tests++;

    bool t12 = rqm_test_theorem_12();
    printf("  [Test 9/9] Theorem 12 (Quantum States as 3 Equilibria):  %s\n", t12 ? ANSI_GREEN "PASSED" ANSI_RESET : ANSI_RED "FAILED" ANSI_RESET);
    if (t12) sum.passed_tests++; else sum.failed_tests++;

    printf("\n" ANSI_BOLD ANSI_CYAN "--------------------------------------------------------------------------------\n" ANSI_RESET);
    printf("  TOTAL: %d | " ANSI_GREEN "PASSED: %d" ANSI_RESET " | " ANSI_RED "FAILED: %d" ANSI_RESET "\n", sum.total_tests, sum.passed_tests, sum.failed_tests);
    printf(ANSI_BOLD ANSI_CYAN "================================================================================\n\n" ANSI_RESET);

    return sum;
}
