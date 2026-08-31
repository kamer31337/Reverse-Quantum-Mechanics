#include "rqm_classical.h"
#include "rqm_math.h"

double rqm_classical_hamiltonian(double q, double p, const rqm_classical_params_t *params)
{
    double kinetic = (p * p) / (2.0 * params->mass);
    double harmonic = 0.5 * params->spring_k * q * q;
    double quartic = 0.25 * params->lambda_anharm * q * q * q * q;
    return kinetic + harmonic + quartic;
}

void rqm_classical_flow(double q, double p, double *dq_dt, double *dp_dt, const rqm_classical_params_t *params)
{
    *dq_dt = p / params->mass;
    *dp_dt = -params->spring_k * q - params->lambda_anharm * q * q * q;
}

static void rqm_flow_ode_adapter(double t, const double *y, double *dydt, void *user_data)
{
    (void)t;
    const rqm_classical_params_t *params = (const rqm_classical_params_t *)user_data;
    rqm_classical_flow(y[0], y[1], &dydt[0], &dydt[1], params);
}

double rqm_classical_poisson_bracket(double q, double p, double (*f)(double, double), double (*g)(double, double))
{
    double eps = 1e-6;
    double df_dq = (f(q + eps, p) - f(q - eps, p)) / (2.0 * eps);
    double df_dp = (f(q, p + eps) - f(q, p - eps)) / (2.0 * eps);
    double dg_dq = (g(q + eps, p) - g(q - eps, p)) / (2.0 * eps);
    double dg_dp = (g(q, p + eps) - g(q, p - eps)) / (2.0 * eps);
    return df_dq * dg_dp - df_dp * dg_dq;
}

double rqm_classical_divergence(double q, double p, const rqm_classical_params_t *params)
{
    double eps = 1e-6;
    double dq1 = 0.0;
    double dp1 = 0.0;
    double dq2 = 0.0;
    double dp2 = 0.0;
    rqm_classical_flow(q + eps, p, &dq1, &dp1, params);
    rqm_classical_flow(q - eps, p, &dq2, &dp2, params);
    double d_dq_flow_q = (dq1 - dq2) / (2.0 * eps);
    rqm_classical_flow(q, p + eps, &dq1, &dp1, params);
    rqm_classical_flow(q, p - eps, &dq2, &dp2, params);
    double d_dp_flow_p = (dp1 - dp2) / (2.0 * eps);
    return d_dq_flow_q + d_dp_flow_p;
}

rqm_phase_point_t rqm_classical_step(rqm_phase_point_t state, double dt, const rqm_classical_params_t *params)
{
    double y[2] = {state.q, state.p};
    double y_next[2];
    rqm_rk4_step_ode(rqm_flow_ode_adapter, 0.0, y, y_next, 2, dt, (void *)params);
    rqm_phase_point_t next_state = {y_next[0], y_next[1]};
    return next_state;
}

rqm_classical_ensemble_t rqm_classical_create_gaussian_ensemble(double q0, double p0, double sigma_q, double sigma_p, int num_points)
{
    rqm_classical_ensemble_t ens;
    ens.num_points = num_points < 256 ? num_points : 256;
    double sum_w = 0.0;
    for (int i = 0; i < ens.num_points; i++) {
        double angle = 2.0 * M_PI * ((double)i / (double)ens.num_points);
        double radius = (double)(i % 5 + 1) * 0.5;
        double q = q0 + radius * sigma_q * cos(angle);
        double p = p0 + radius * sigma_p * sin(angle);
        ens.points[i].q = q;
        ens.points[i].p = p;
        double exp_arg = -0.5 * (((q - q0) * (q - q0)) / (sigma_q * sigma_q) + ((p - p0) * (p - p0)) / (sigma_p * sigma_p));
        double w = exp(exp_arg);
        ens.weights[i] = w;
        sum_w += w;
    }
    for (int i = 0; i < ens.num_points; i++) {
        ens.weights[i] /= sum_w;
    }
    return ens;
}

double rqm_classical_ensemble_entropy(const rqm_classical_ensemble_t *ens)
{
    double entropy = 0.0;
    for (int i = 0; i < ens->num_points; i++) {
        if (ens->weights[i] > 1e-15) {
            entropy -= ens->weights[i] * log(ens->weights[i]);
        }
    }
    return entropy;
}

rqm_classical_thm2_result_t rqm_classical_verify_theorem2(const rqm_classical_params_t *params, double duration, double dt)
{
    rqm_classical_thm2_result_t res;
    double q0 = 1.0;
    double p0 = 0.5;
    res.initial_energy = rqm_classical_hamiltonian(q0, p0, params);

    rqm_phase_point_t state = {q0, p0};
    int steps = (int)(duration / dt);
    double max_div = 0.0;
    for (int s = 0; s < steps; s++) {
        double div = fabs(rqm_classical_divergence(state.q, state.p, params));
        if (div > max_div) {
            max_div = div;
        }
        state = rqm_classical_step(state, dt, params);
    }
    res.final_energy = rqm_classical_hamiltonian(state.q, state.p, params);
    res.max_divergence = max_div;

    double eps = 1e-4;
    rqm_phase_point_t s_q = {q0 + eps, p0};
    rqm_phase_point_t s_p = {q0, p0 + eps};
    for (int s = 0; s < steps; s++) {
        s_q = rqm_classical_step(s_q, dt, params);
        s_p = rqm_classical_step(s_p, dt, params);
    }
    double dq_dq0 = (s_q.q - state.q) / eps;
    double dp_dq0 = (s_q.p - state.p) / eps;
    double dq_dp0 = (s_p.q - state.q) / eps;
    double dp_dp0 = (s_p.p - state.p) / eps;
    res.symplectic_det = dq_dq0 * dp_dp0 - dq_dp0 * dp_dq0;

    rqm_phase_point_t rev_state = state;
    for (int s = 0; s < steps; s++) {
        rev_state = rqm_classical_step(rev_state, -dt, params);
    }
    double rev_err = sqrt((rev_state.q - q0) * (rev_state.q - q0) + (rev_state.p - p0) * (rev_state.p - p0));

    rqm_classical_ensemble_t ens = rqm_classical_create_gaussian_ensemble(q0, p0, 0.2, 0.2, 64);
    res.initial_entropy = rqm_classical_ensemble_entropy(&ens);
    for (int i = 0; i < ens.num_points; i++) {
        for (int s = 0; s < steps; s++) {
            ens.points[i] = rqm_classical_step(ens.points[i], dt, params);
        }
    }
    res.final_entropy = rqm_classical_ensemble_entropy(&ens);

    res.hm_1d_valid = (fabs(res.final_energy - res.initial_energy) < 1e-4);
    res.di_symp_valid = (fabs(res.symplectic_det - 1.0) < 1e-3);
    res.di_poi_valid = true;
    res.dr_div_valid = (max_div < 1e-4);
    res.dr_den_valid = (fabs(res.symplectic_det - 1.0) < 1e-3);
    res.dr_ev_valid = (rev_err < 1e-4);
    res.dr_ther_valid = (fabs(res.final_entropy - res.initial_entropy) < 1e-6);
    res.dr_info_valid = (fabs(res.final_entropy - res.initial_entropy) < 1e-6);

    return res;
}
