#ifndef RQM_HILBERT_PATHOLOGY_H
#define RQM_HILBERT_PATHOLOGY_H

#include "rqm_common.h"

typedef struct {
    int max_j;
    double j_values[64];
    double norm_diff[64];
    double expectation_N[64];
    bool cobs_violated;
} rqm_cobs_demo_t;

typedef struct {
    double x_grid_min;
    double x_grid_max;
    int grid_points;
    double variance_x_frame;
    double variance_y_frame;
    bool is_y_divergent;
    bool wdq_violated;
} rqm_wdq_demo_t;

typedef struct {
    int dim1;
    int dim2;
    bool isomorphic_as_hilbert_spaces;
    bool distinguishable_in_qsto_v;
    const char *dim1_lie_algebra;
    const char *dim2_lie_algebra;
} rqm_stdist_demo_t;

typedef struct {
    bool theorem3_inconsistency_demonstrated;
    bool theorem4_qsto_v_resolved;
    rqm_cobs_demo_t cobs;
    rqm_wdq_demo_t wdq;
    rqm_stdist_demo_t stdist;
} rqm_pathology_report_t;

rqm_cobs_demo_t rqm_demonstrate_cobs_violation(int max_j);
rqm_wdq_demo_t rqm_demonstrate_wdq_violation(void);
rqm_stdist_demo_t rqm_demonstrate_stdist_violation(int dof1, int dof2);
rqm_pathology_report_t rqm_verify_theorems_3_and_4(void);

#endif
