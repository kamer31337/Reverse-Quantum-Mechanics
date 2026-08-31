#include "rqm_ui.h"

void rqm_ui_clear_screen(void)
{
    printf(ANSI_CLEAR);
}

void rqm_ui_print_header(const char *title, const char *subtitle)
{
    printf(ANSI_CYAN ANSI_BOLD "================================================================================\n" ANSI_RESET);
    printf(ANSI_WHITE ANSI_BOLD "  REVERSE QUANTUM MECHANICS (RQM) - arXiV:2608.27543v1\n" ANSI_RESET);
    printf(ANSI_YELLOW "  %s\n" ANSI_RESET, title);
    if (subtitle && strlen(subtitle) > 0) {
        printf(ANSI_DIM "  %s\n" ANSI_RESET, subtitle);
    }
    printf(ANSI_CYAN ANSI_BOLD "================================================================================\n\n" ANSI_RESET);
}

void rqm_ui_print_matrix(const char *label, const rqm_mat_t *mat)
{
    printf(ANSI_MAGENTA ANSI_BOLD "[Matrix: %s (%dx%d)]\n" ANSI_RESET, label, mat->dim, mat->dim);
    for (int i = 0; i < mat->dim; i++) {
        printf("  [ ");
        for (int j = 0; j < mat->dim; j++) {
            double re = creal(mat->data[i][j]);
            double im = cimag(mat->data[i][j]);
            if (fabs(im) < 1e-6) {
                printf(ANSI_CYAN "%8.4f       " ANSI_RESET, re);
            } else if (im >= 0.0) {
                printf(ANSI_CYAN "%7.4f" ANSI_YELLOW "+%6.4fi " ANSI_RESET, re, im);
            } else {
                printf(ANSI_CYAN "%7.4f" ANSI_YELLOW "-%6.4fi " ANSI_RESET, re, fabs(im));
            }
        }
        printf(" ]\n");
    }
    printf("\n");
}

void rqm_ui_print_vector(const char *label, const rqm_vec_t *vec)
{
    printf(ANSI_MAGENTA ANSI_BOLD "[Vector: %s (dim=%d)]\n" ANSI_RESET, label, vec->dim);
    printf("  [ ");
    for (int i = 0; i < vec->dim; i++) {
        double re = creal(vec->data[i]);
        double im = cimag(vec->data[i]);
        if (fabs(im) < 1e-6) {
            printf(ANSI_CYAN "%7.4f " ANSI_RESET, re);
        } else if (im >= 0.0) {
            printf(ANSI_CYAN "%6.3f" ANSI_YELLOW "+%5.3fi " ANSI_RESET, re, im);
        } else {
            printf(ANSI_CYAN "%6.3f" ANSI_YELLOW "-%5.3fi " ANSI_RESET, re, fabs(im));
        }
    }
    printf(" ]^T\n\n");
}

void rqm_ui_render_bloch_sphere(double x, double y, double z)
{
    double r = sqrt(x * x + y * y + z * z);
    printf(ANSI_BOLD ANSI_CYAN "  Bloch Sphere Representation (|r| = %.4f %s):\n" ANSI_RESET, r, r > 1.0001 ? "(Unphysical)" : (fabs(r - 1.0) < 0.01 ? "Pure State" : "Mixed State"));
    printf("        " ANSI_YELLOW "+z |0>" ANSI_RESET "          r_x = %+.4f\n", x);
    printf("           .---.            r_y = %+.4f\n", y);
    printf("        .-'  |  '-.         r_z = %+.4f\n", z);
    printf("       /     |     \\        Purity Tr(rho^2) = %.4f\n", 0.5 * (1.0 + r * r));
    printf("      |   ---+---   |       Entropy S_2(rho) = %.4f bits\n", (r >= 1.0) ? 0.0 : (-0.5 * (1.0 + r) * (log(0.5 * (1.0 + r)) / log(2.0)) - 0.5 * (1.0 - r) * (log(0.5 * (1.0 - r)) / log(2.0))));
    printf("       \\     |     /\n");
    printf("        '-.  |  .-'\n");
    printf("           '---'\n");
    printf("        " ANSI_YELLOW "-z |1>" ANSI_RESET "\n\n");
}

void rqm_ui_plot_curve(const char *title, const double *x_vals, const double *y_vals, int count, int width, int height, const char *x_label, const char *y_label)
{
    if (count <= 0 || width <= 10 || height <= 4) return;
    double min_x = x_vals[0];
    double max_x = x_vals[0];
    double min_y = y_vals[0];
    double max_y = y_vals[0];

    for (int i = 1; i < count; i++) {
        if (x_vals[i] < min_x) min_x = x_vals[i];
        if (x_vals[i] > max_x) max_x = x_vals[i];
        if (y_vals[i] < min_y) min_y = y_vals[i];
        if (y_vals[i] > max_y) max_y = y_vals[i];
    }
    if (fabs(max_y - min_y) < 1e-10) max_y += 1.0;
    if (fabs(max_x - min_x) < 1e-10) max_x += 1.0;

    printf(ANSI_BOLD ANSI_CYAN "  [Graph: %s]\n" ANSI_RESET, title);
    char grid[32][80];
    int plot_h = height < 30 ? height : 30;
    int plot_w = width < 70 ? width : 70;

    for (int r = 0; r < plot_h; r++) {
        for (int c = 0; c < plot_w; c++) {
            grid[r][c] = ' ';
        }
    }

    for (int i = 0; i < count; i++) {
        int c = (int)(((x_vals[i] - min_x) / (max_x - min_x)) * (double)(plot_w - 1));
        int r = (int)(((max_y - y_vals[i]) / (max_y - min_y)) * (double)(plot_h - 1));
        if (r >= 0 && r < plot_h && c >= 0 && c < plot_w) {
            grid[r][c] = '*';
        }
    }

    for (int r = 0; r < plot_h; r++) {
        double y_val = max_y - ((double)r / (double)(plot_h - 1)) * (max_y - min_y);
        printf(ANSI_DIM "%7.3f | " ANSI_RESET, y_val);
        for (int c = 0; c < plot_w; c++) {
            if (grid[r][c] == '*') {
                printf(ANSI_GREEN ANSI_BOLD "*" ANSI_RESET);
            } else {
                printf(" ");
            }
        }
        printf("\n");
    }

    printf("        +");
    for (int c = 0; c < plot_w; c++) printf("-");
    printf("\n         ");
    printf(ANSI_DIM "%-8.3f" ANSI_RESET, min_x);
    for (int c = 0; c < plot_w - 16; c++) printf(" ");
    printf(ANSI_DIM "%8.3f" ANSI_RESET, max_x);
    printf("\n         " ANSI_YELLOW "%s -> (%s on y-axis)\n\n" ANSI_RESET, x_label, y_label);
}

void rqm_ui_render_phase_space(const double *q_vals, const double *p_vals, int count, int width, int height)
{
    rqm_ui_plot_curve("Phase Space Orbit (q, p)", q_vals, p_vals, count, width, height, "Coordinate q", "Momentum p");
}

void rqm_ui_render_bar_chart(const char *title, const char **labels, const double *values, int count, double max_val)
{
    printf(ANSI_BOLD ANSI_CYAN "  [Bar Chart: %s]\n" ANSI_RESET, title);
    for (int i = 0; i < count; i++) {
        printf("  %-16s | ", labels[i]);
        int bar_len = (int)((values[i] / max_val) * 40.0);
        if (bar_len < 0) bar_len = 0;
        if (bar_len > 40) bar_len = 40;
        printf(ANSI_GREEN);
        for (int b = 0; b < bar_len; b++) printf("=");
        printf(ANSI_RESET " %.4f\n", values[i]);
    }
    printf("\n");
}

void rqm_ui_press_enter_to_continue(void)
{
    printf(ANSI_YELLOW "  Press ENTER to continue..." ANSI_RESET);
    fflush(stdout);
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}
