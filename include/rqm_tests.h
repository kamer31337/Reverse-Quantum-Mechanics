#ifndef RQM_TESTS_H
#define RQM_TESTS_H

#include "rqm_common.h"

typedef struct {
    int total_tests;
    int passed_tests;
    int failed_tests;
    double max_discrepancy;
} rqm_test_summary_t;

bool rqm_test_math_suite(void);
bool rqm_test_theorem_2(void);
bool rqm_test_theorems_3_and_4(void);
bool rqm_test_theorems_5_and_6(void);
bool rqm_test_theorem_7(void);
bool rqm_test_theorem_8(void);
bool rqm_test_theorems_9_and_10(void);
bool rqm_test_theorem_11(void);
bool rqm_test_theorem_12(void);
rqm_test_summary_t rqm_test_run_all(void);

#endif
