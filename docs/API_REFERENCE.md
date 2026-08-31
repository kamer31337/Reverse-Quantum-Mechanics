# Reverse Quantum Mechanics (RQM) API Reference

This document provides a detailed API reference for the C11 Reverse Quantum Mechanics software library.

---

## 1. Core Types & Constants (`rqm_common.h`)

```c
typedef double complex rqm_cdouble;

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
```

---

## 2. Linear Algebra & Quantum Information (`rqm_math.h`)

### Constructors & Standard Operators
- `rqm_mat_t rqm_mat_zero(int dim)`: Returns a zero matrix of given dimension.
- `rqm_mat_t rqm_mat_identity(int dim)`: Returns the $N\times N$ identity matrix.
- `rqm_mat_t rqm_mat_pauli_x(void)`: Returns Pauli $\sigma_x$.
- `rqm_mat_t rqm_mat_pauli_y(void)`: Returns Pauli $\sigma_y$.
- `rqm_mat_t rqm_mat_pauli_z(void)`: Returns Pauli $\sigma_z$.
- `rqm_mat_t rqm_mat_annihilation(int dim)`: Returns truncated bosonic ladder operator $a$.
- `rqm_mat_t rqm_mat_creation(int dim)`: Returns truncated bosonic ladder operator $a^\dagger$.
- `rqm_mat_t rqm_mat_position(int dim)`: Returns position observable $X = \frac{a + a^\dagger}{\sqrt{2}}$.
- `rqm_mat_t rqm_mat_momentum(int dim)`: Returns momentum observable $P = \frac{a - a^\dagger}{i\sqrt{2}}$.
- `rqm_mat_t rqm_mat_number(int dim)`: Returns number operator $N = a^\dagger a$.

### State Vectors
- `rqm_vec_t rqm_vec_qubit(double theta, double phi)`: Parametrizes $|\psi\rangle = \cos(\theta/2)|0\rangle + e^{i\phi}\sin(\theta/2)|1\rangle$.
- `rqm_cdouble rqm_vec_dot(const rqm_vec_t *u, const rqm_vec_t *v)`: Computes $\langle u | v \rangle$.
- `rqm_mat_t rqm_vec_outer(const rqm_vec_t *u, const rqm_vec_t *v)`: Computes $|u\rangle\langle v|$.

### Matrix Operations & Norms
- `rqm_mat_t rqm_mat_commutator(const rqm_mat_t *A, const rqm_mat_t *B)`: $[A, B] = AB - BA$.
- `rqm_mat_t rqm_mat_anticommutator(const rqm_mat_t *A, const rqm_mat_t *B)`: $\{A, B\} = AB + BA$.
- `rqm_cdouble rqm_mat_hs_dot(const rqm_mat_t *A, const rqm_mat_t *B)`: $\langle A, B \rangle_{HS} = \operatorname{Tr}(A^\dagger B)$.
- `double rqm_mat_hs_norm(const rqm_mat_t *A)`: $\|A\|_{HS} = \sqrt{\operatorname{Tr}(A^\dagger A)}$.
- `double rqm_mat_hs_dist(const rqm_mat_t *A, const rqm_mat_t *B)`: $\|A - B\|_{HS}$.

### Eigensolvers & Spectral Functions
- `rqm_status_t rqm_hermitian_eigensolve(const rqm_mat_t *A, rqm_eigen_t *out_eigen)`: Diagonalizes complex Hermitian matrix via Jacobi rotations.
- `rqm_mat_t rqm_mat_exp_hermitian(const rqm_mat_t *H, rqm_cdouble scale)`: Computes $e^{s H}$.
- `rqm_mat_t rqm_mat_log_density(const rqm_mat_t *rho, double zero_cutoff)`: Computes $\ln \rho$ restricted to support.
- `double rqm_von_neumann_entropy(const rqm_mat_t *rho)`: Computes $S(\rho) = -\operatorname{Tr}(\rho \ln \rho)$ in nats.
- `double rqm_von_neumann_entropy_base2(const rqm_mat_t *rho)`: Computes $S(\rho)$ in bits.
- `double rqm_relative_entropy(const rqm_mat_t *rho, const rqm_mat_t *sigma)`: $S(\rho \| \sigma) = \operatorname{Tr}(\rho \ln \rho - \rho \ln \sigma)$.
- `double rqm_purity(const rqm_mat_t *rho)`: Computes $\operatorname{Tr}(\rho^2)$.

---

## 3. Classical 1-DOF Module (`rqm_classical.h`)

- `double rqm_classical_hamiltonian(double q, double p, const rqm_classical_params_t *params)`: Evaluates $H(q, p) = \frac{p^2}{2m} + \frac{k}{2}q^2 + \frac{\lambda}{4}q^4$.
- `void rqm_classical_flow(double q, double p, double *dq_dt, double *dp_dt, const rqm_classical_params_t *params)`: Phase space velocity vector field.
- `double rqm_classical_poisson_bracket(double q, double p, double (*f)(double, double), double (*g)(double, double))`: Finite-difference $\{f, g\}$.
- `double rqm_classical_divergence(double q, double p, const rqm_classical_params_t *params)`: Numerical divergence $\partial_q \dot{q} + \partial_p \dot{p}$.
- `rqm_classical_thm2_result_t rqm_classical_verify_theorem2(const rqm_classical_params_t *params, double duration, double dt)`: Verifies all 8 equivalent conditions of Theorem 2.

---

## 4. Hilbert Space Pathologies & QSTO:V (`rqm_hilbert_pathology.h`)

- `rqm_cobs_demo_t rqm_demonstrate_cobs_violation(int max_j)`: Evaluates the sequence $|\psi_j\rangle = \sqrt{\frac{j-1}{j}}|0\rangle + \sqrt{\frac{1}{j}}|j\rangle$, showing norm difference $\to 0$ while $\langle N \rangle = 1 \neq 0$.
- `rqm_wdq_demo_t rqm_demonstrate_wdq_violation(void)`: Integrates transformed wavepacket under $y(x) = \tan(\frac{\pi}{2}\operatorname{erf}(x))$, showing divergent $\langle Y^2 \rangle$.
- `rqm_stdist_demo_t rqm_demonstrate_stdist_violation(int dof1, int dof2)`: Analyzes Hilbert space isomorphism $L^2(\mathbb{R}^{d_1}) \cong L^2(\mathbb{R}^{d_2})$ vs Lie algebra separation.
- `rqm_pathology_report_t rqm_verify_theorems_3_and_4(void)`: Automated verification for Theorems 3 and 4.

---

## 5. Pure States vs. Convex Ensembles (`rqm_ensembles.h`)

- `rqm_mat_t rqm_bloch_to_density(rqm_bloch_vec_t b)`: Constructs $\rho = \frac{1}{2}(I + \vec{r}\cdot\vec{\sigma})$.
- `rqm_bloch_vec_t rqm_density_to_bloch(const rqm_mat_t *rho)`: Extracts Bloch vector $r_i = \operatorname{Tr}(\rho \sigma_i)$.
- `rqm_ensemble_decomposition_t rqm_decompose_spectral(const rqm_mat_t *rho)`: Spectral pure-state decomposition.
- `rqm_ensemble_decomposition_t rqm_decompose_random_pure(const rqm_mat_t *rho, double angle_offset)`: Generates non-orthogonal 3-state convex decomposition of the same mixed state.
- `rqm_ensembles_thm5_6_result_t rqm_verify_theorems_5_and_6(rqm_bloch_vec_t target_bloch)`: Verifies pure-state symplectic equivalence and multi-decomposition.

---

## 6. Born Rule & Entropy Equivalence (`rqm_born_entropy.h`)

- `double rqm_born_transition_probability(const rqm_vec_t *psi, const rqm_vec_t *phi)`: $p = |\langle\phi|\psi\rangle|^2$.
- `double rqm_entropy_of_equal_mixture_analytic(double p)`: $S(p) = -\lambda_+\log_2 \lambda_+ - \lambda_-\log_2 \lambda_-$ with $\lambda_\pm = \frac{1\pm\sqrt{p}}{2}$.
- `double rqm_invert_entropy_to_probability(double entropy_bits)`: Reconstructs transition probability $p$ from entropy $S$ via numerical inversion.
- `rqm_born_thm7_result_t rqm_verify_theorem7(int num_samples)`: Validates monotonicity, boundary values ($p=0 \implies S=1$, $p=1 \implies S=0$), and numerical match.

---

## 7. Unitary Evolution (`rqm_unitary.h`)

- `rqm_mat_t rqm_unitary_propagator(const rqm_mat_t *H, double dt, double hbar)`: $U(dt) = e^{-iHdt/\hbar}$.
- `rqm_vec_t rqm_evolve_state_unitary(const rqm_vec_t *psi, const rqm_mat_t *H, double t, double hbar)`: State vector propagator.
- `rqm_mat_t rqm_evolve_density_unitary(const rqm_mat_t *rho, const rqm_mat_t *H, double t, double hbar)`: Density matrix propagator $U\rho U^\dagger$.
- `rqm_unitary_thm8_result_t rqm_verify_theorem8(...)`: Validates entropy conservation, probability conservation, and reversibility.

---

## 8. Equilibrations & Zeno Evolution (`rqm_equilibration.h`)

- `rqm_mat_t rqm_projective_measurement_nonselective(const rqm_mat_t *rho, const rqm_mat_t *X)`: $\Phi_X(\rho) = \sum_k P_k \rho P_k$.
- `rqm_mat_t rqm_lindblad_dephasing_rhs(const rqm_mat_t *rho, const rqm_mat_t *X, double gamma)`: Purely dissipative Lindblad generator $-\frac{\gamma}{2}[X, [X, \rho]]$.
- `rqm_lindblad_equilibration_result_t rqm_simulate_lindblad_equilibration(...)`: Integrates dephasing ODE to asymptotic equilibrium.
- `rqm_zeno_unitary_result_t rqm_simulate_zeno_unitary(const rqm_vec_t *psi0, const rqm_mat_t *H, double total_time)`: Computes $N$-step infinitesimal projective sequence converging to unitary dynamics.

---

## 9. Classical Limit as High Entropy (`rqm_classical_limit.h`)

- `double rqm_wigner_fock_state(int n, double q, double p, double hbar)`: Computes Fock state Wigner quasiprobability.
- `double rqm_wigner_damped(int n, double q, double p, double lambda, double hbar)`: Computes smoothed Wigner distribution under $L = a^\dagger$ channel.
- `rqm_classical_limit_result_t rqm_verify_theorem11(...)`: Validates non-negativity recovery and vanishing commutator $\hbar/\lambda \to 0$.

---

## 10. Triple Equilibria (`rqm_triple_equilibria.h`)

- `rqm_dynamic_equilibrium_t rqm_verify_dynamic_equilibrium(const rqm_mat_t *rho)`: Verifies $[H_{\text{dyn}}, \rho] = 0$ with $H_{\text{dyn}} = \rho$.
- `rqm_spectral_equilibrium_t rqm_verify_spectral_equilibrium(const rqm_mat_t *rho)`: Verifies $\Phi_{X}(\rho) = \rho$ with $X = \rho$.
- `rqm_thermodynamic_equilibrium_t rqm_verify_thermodynamic_equilibrium(const rqm_mat_t *rho, double beta)`: Constructs Gibbs Hamiltonian $H_{\text{ther}}$ (with kernel penalty for rank-deficient states) and verifies $\rho_{\text{Gibbs}} = \rho$.
- `rqm_triple_equilibria_result_t rqm_verify_theorem12(const rqm_mat_t *rho, double beta)`: Evaluates all three equilibria simultaneously.
