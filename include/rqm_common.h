#ifndef RQM_COMMON_H
#define RQM_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#define RQM_EPSILON 1e-10
#define RQM_MAX_DIM 16
#define RQM_DEFAULT_HBAR 1.0

typedef double complex rqm_cdouble;

typedef struct {
    double real;
    double imag;
} rqm_complex_t;

typedef enum {
    RQM_SUCCESS = 0,
    RQM_ERROR_INVALID_PARAM = -1,
    RQM_ERROR_CONVERGENCE = -2,
    RQM_ERROR_DIM_MISMATCH = -3,
    RQM_ERROR_NOT_HERMITIAN = -4,
    RQM_ERROR_NOT_NORMALIZED = -5
} rqm_status_t;

typedef struct {
    int dim;
    rqm_cdouble data[RQM_MAX_DIM];
} rqm_vec_t;

typedef struct {
    int dim;
    rqm_cdouble data[RQM_MAX_DIM][RQM_MAX_DIM];
} rqm_mat_t;

typedef struct {
    int dim;
    double values[RQM_MAX_DIM];
    rqm_vec_t vectors[RQM_MAX_DIM];
} rqm_eigen_t;

#endif
