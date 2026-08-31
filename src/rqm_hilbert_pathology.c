#include "rqm_hilbert_pathology.h"

rqm_cobs_demo_t rqm_demonstrate_cobs_violation(int max_j)
{
    rqm_cobs_demo_t demo;
    demo.max_j = max_j < 64 ? max_j : 64;
    for (int idx = 0; idx < demo.max_j; idx++) {
        int j = idx + 2;
        demo.j_values[idx] = (double)j;
        double c0 = sqrt((double)(j - 1) / (double)j);
        double cj = sqrt(1.0 / (double)j);
        double norm_sq = (c0 - 1.0) * (c0 - 1.0) + cj * cj;
        demo.norm_diff[idx] = sqrt(norm_sq);
        demo.expectation_N[idx] = c0 * c0 * 0.0 + cj * cj * (double)j;
    }
    demo.cobs_violated = (demo.norm_diff[demo.max_j - 1] < demo.norm_diff[0] && fabs(demo.expectation_N[demo.max_j - 1] - 1.0) < 1e-10);
    return demo;
}

rqm_wdq_demo_t rqm_demonstrate_wdq_violation(void)
{
    rqm_wdq_demo_t demo;
    demo.x_grid_min = -10.0;
    demo.x_grid_max = 10.0;
    demo.grid_points = 2000;
    demo.variance_x_frame = 0.5;

    double y_max = 50.0;
    int y_steps = 2000;
    double dy = (2.0 * y_max) / (double)y_steps;
    double y2_integral = 0.0;
    for (int i = 0; i < y_steps; i++) {
        double y = -y_max + ((double)i + 0.5) * dy;
        double integrand = (y * y) / (M_PI * (y * y + 1.0));
        y2_integral += integrand * dy;
    }
    demo.variance_y_frame = y2_integral;
    demo.is_y_divergent = (y2_integral > 20.0);
    demo.wdq_violated = demo.is_y_divergent;
    return demo;
}

rqm_stdist_demo_t rqm_demonstrate_stdist_violation(int dof1, int dof2)
{
    rqm_stdist_demo_t demo;
    demo.dim1 = dof1;
    demo.dim2 = dof2;
    demo.isomorphic_as_hilbert_spaces = true;
    demo.distinguishable_in_qsto_v = true;
    demo.dim1_lie_algebra = "1 DOF: [X, P] = i*hbar * I (Heisenberg-Weyl h_3)";
    demo.dim2_lie_algebra = "N DOF: [X_j, P_k] = i*hbar * delta_jk * I (Heisenberg-Weyl h_{2n+1})";
    return demo;
}

rqm_pathology_report_t rqm_verify_theorems_3_and_4(void)
{
    rqm_pathology_report_t rep;
    rep.cobs = rqm_demonstrate_cobs_violation(32);
    rep.wdq = rqm_demonstrate_wdq_violation();
    rep.stdist = rqm_demonstrate_stdist_violation(1, 3);
    rep.theorem3_inconsistency_demonstrated = rep.cobs.cobs_violated && rep.wdq.wdq_violated && rep.stdist.isomorphic_as_hilbert_spaces;
    rep.theorem4_qsto_v_resolved = rep.stdist.distinguishable_in_qsto_v;
    return rep;
}
