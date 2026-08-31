#ifndef RQM_UI_H
#define RQM_UI_H

#include "rqm_common.h"

#define ANSI_RESET   "\033[0m"
#define ANSI_BOLD    "\033[1m"
#define ANSI_DIM     "\033[2m"
#define ANSI_RED     "\033[31m"
#define ANSI_GREEN   "\033[32m"
#define ANSI_YELLOW  "\033[33m"
#define ANSI_BLUE    "\033[34m"
#define ANSI_MAGENTA "\033[35m"
#define ANSI_CYAN    "\033[36m"
#define ANSI_WHITE   "\033[37m"
#define ANSI_BG_BLUE "\033[44m"
#define ANSI_CLEAR   "\033[2J\033[H"

void rqm_ui_clear_screen(void);
void rqm_ui_print_header(const char *title, const char *subtitle);
void rqm_ui_print_matrix(const char *label, const rqm_mat_t *mat);
void rqm_ui_print_vector(const char *label, const rqm_vec_t *vec);
void rqm_ui_render_bloch_sphere(double x, double y, double z);
void rqm_ui_plot_curve(const char *title, const double *x_vals, const double *y_vals, int count, int width, int height, const char *x_label, const char *y_label);
void rqm_ui_render_phase_space(const double *q_vals, const double *p_vals, int count, int width, int height);
void rqm_ui_render_bar_chart(const char *title, const char **labels, const double *values, int count, double max_val);
void rqm_ui_press_enter_to_continue(void);

#endif
