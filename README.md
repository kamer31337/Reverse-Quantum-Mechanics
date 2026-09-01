# Reverse Quantum Mechanics (RQM)

A self-contained C11 scientific simulation engine, theoretical framework, and interactive terminal program formalizing **Reverse Quantum Mechanics** (arXiv:2608.27543v1 [quant-ph], 27 Aug 2026 by Gabriele Carcassi, Tobias Thrien, and Christine A. Aidala).

---

## Highlights & Features

- **Self-Contained C11 Mathematics**: Custom complex linear algebra, Hermitian Jacobi eigensolvers, matrix exponentials/logarithms, Runge-Kutta 4th-order ODE integrators, and von Neumann / quantum relative entropy functions with zero external dependencies.
- **Interactive Terminal UI**: ANSI color-coded dashboard with ASCII visualizers for:
  - 3D Bloch sphere projections ($r_x, r_y, r_z$),
  - Classical phase space orbits $(q, p)$,
  - Lindblad dephasing trajectories ($\|\rho(t) - \Phi_X(\rho)\|_{HS} \to 0$),
  - Monotonic Born rule vs. von Neumann entropy curves ($S(p) \leftrightarrow p$),
  - High-entropy Wigner function negativity suppression under $L = a^\dagger$.
- **Automated Validation Suite**: Numerical proof checks for all 12 theorems presented in the paper.
- **Cross-Platform**: Compiles out-of-the-box with Clang, GCC, and MSVC on Linux, macOS, and Windows.

---

## Quick Start

### 1. Build using Clang or GCC
```bash
# Using Clang
clang -std=c11 -O3 -Wall -Wextra -Iinclude src/*.c -lm -o rqm

# Using GCC
gcc -std=c11 -O3 -Wall -Wextra -Iinclude src/*.c -lm -o rqm

# Or use Make
make
```

### 2. Windows Batch or PowerShell Build
```cmd
# Windows Command Prompt
build.bat

# Windows PowerShell
.\build.ps1
```

### 3. Run the Program

```bash
# Launch interactive terminal UI
./rqm

# Run automated theorem verification test suite
./rqm --test

# Run specific theorem demonstration (e.g., Theorem 7 Born rule / entropy)
./rqm --thm 7
```

---

## Implemented Theorems Overview

| Module / Theorem | Core Finding & Reverse Physics Logical Relationship |
|---|---|
| **Theorem 2** (Classical 1-DOF) | Over $\text{CST-1D}$: $\text{HM-1D} \iff \text{DI-SYMP} \iff \text{DI-POI} \iff \text{DR-DIV} \iff \text{DR-DEN} \iff \text{DR-EV} \iff \text{DR-THER} \iff \text{DR-INFO}$. |
| **Theorems 3 & 4** (Hilbert Space Pathology) | Standard $\text{QSTO:H}$ is inconsistent with base conditions $[B]\text{COBS}$ (discontinuity of $\langle N \rangle$), $[B]\text{WDQ}$ (domain divergence under coordinate transforms), and $[B]\text{STDIST}$ (separable space isomorphism). Resolved by Quantum Vector Space $\text{QSTO:V}$. |
| **Theorems 5 & 6** (Pure States vs Ensembles) | Pure states on $\mathbb{C}P^1 \cong S^2$ are symplectic manifolds and cannot exclude classical mechanics. Quantum systems are distinguished by their convex ensemble space $\mathcal{D}(\mathcal{H})$, allowing infinitely many equivalent pure-state decompositions. |
| **Theorem 7** (Born Rule & Entropy) | $\text{BR-BORN} \iff \text{BR-ENT} \iff \text{BR-ORME}$. The transition probability $p$ is strictly monotonically mapped to the von Neumann entropy of equal mixtures $S(p) = -\lambda_+\log_2\lambda_+ - \lambda_-\log_2\lambda_-$. |
| **Theorem 8** (Unitary Evolution) | $\text{DR-SCEQ} \iff \text{DR-UNIT} \iff \text{DR-EV} \iff \text{DR-PROB} \iff \text{DR-INFO}$. Schrödinger evolution is equivalent to deterministic and reversible information conservation. |
| **Theorems 9 & 10** (Equilibrations & Zeno Limit) | Projective measurements are dissipative Lindblad equilibrations $\frac{d\rho}{dt} = -\frac{1}{2}[X, [X, \rho]]$ minimizing Hilbert-Schmidt distance ($\text{EQ-HS}$) and relative entropy ($\text{EQ-KL}$). Unitary evolution is the continuous limit of infinitesimal state-adaptive projections. |
| **Theorem 11** (Classical Limit) | $\text{CL-HBAR} \iff \text{CL-HENT} \iff \text{CL-CST}$. An open Lindblad creation channel $L = a^\dagger$ stretches phase space ($\lambda = e^{\gamma t} \to \infty$), suppresses Wigner negativity, and deforms the Moyal bracket into the classical Poisson bracket. |
| **Theorem 12** (Triple Equilibria) | Every quantum state $\rho$ (pure or mixed) is simultaneously a **Dynamic Equilibrium** ($H=\rho$), a **Spectral Equilibrium** ($X=\rho$), and a **Thermodynamic Gibbs Equilibrium** ($H_{\text{ther}} = -\frac{1}{\beta}\ln\rho$). |

---

## Directory Structure

```
Reverse-Quantum-Mechanics/
├── CMakeLists.txt                       # CMake build file
├── Makefile                             # POSIX Makefile
├── build.bat                            # Windows Batch compiler script
├── build.ps1                            # Windows PowerShell compiler script
├── 2608.27543v1.pdf                     # Original Research Paper
├── README.md                            # Main project documentation
├── docs/
│   ├── REVERSE_QUANTUM_MECHANICS_THEORY.md  # Complete mathematical derivations
│   └── API_REFERENCE.md                     # C11 API reference
├── include/
│   ├── rqm_common.h                     # Core types, constants, complex structures
│   ├── rqm_math.h                       # Linear algebra, eigensolvers, entropies, RK4
│   ├── rqm_classical.h                  # Classical 1-DOF symplectic invariants
│   ├── rqm_hilbert_pathology.h          # Unphysicality counterexamples & QSTO:V
│   ├── rqm_ensembles.h                  # Convex ensemble spaces & multi-decompositions
│   ├── rqm_born_entropy.h               # Born rule / mutual exclusivity / entropy bijection
│   ├── rqm_unitary.h                    # Unitary group & information conservation
│   ├── rqm_equilibration.h              # Lindblad equilibration & Quantum Zeno
│   ├── rqm_classical_limit.h            # High-entropy Lindblad channel & Wigner functions
│   ├── rqm_triple_equilibria.h          # Quantum states as triple equilibria
│   ├── rqm_ui.h                         # ANSI terminal dashboard & ASCII visualizer
│   └── rqm_tests.h                      # Automated verification test suite
└── src/
    ├── main.c                           # Interactive CLI & main menu entry point
    ├── rqm_math.c                       # Matrix and quantum info implementation
    ├── rqm_classical.c                  # Classical 1-DOF dynamics
    ├── rqm_hilbert_pathology.c          # Pathology proofs & QSTO:V Lie algebra
    ├── rqm_ensembles.c                  # Bloch ball multi-decomposition generator
    ├── rqm_born_entropy.c               # Invertible entropy-probability map
    ├── rqm_unitary.c                    # Schrödinger & Liouville-von Neumann propagator
    ├── rqm_equilibration.c              # Lindblad solver & Zeno projective limit
    ├── rqm_classical_limit.c            # Wigner function smoothing & Moyal bracket limit
    ├── rqm_triple_equilibria.c          # Triple equilibria constructor
    ├── rqm_ui.c                         # Terminal rendering & plotting logic
    └── rqm_tests.c                      # Automated test suite implementation
```

---

## References

1. Carcassi, G.; Thrien, T.; Aidala, C. A. *Reverse Quantum Mechanics*. arXiv:2608.27543v1 [quant-ph] (2026).
2. Carcassi, G.; Aidala, C. A. *Assumptions of Physics*. Michigan Publishing (2021).
3. Carcassi, G.; Calderón, F.; Aidala, C. A. *The Unphysicality of Hilbert Spaces*. Quantum Studies: Mathematics and Foundations 12, 13 (2025).
4. Carcassi, G.; Landini, M.; Aidala, C. A. *Classical mechanics as the high-entropy limit of quantum mechanics*. Physica Scripta 101, 065105 (2026).

## Mentioned Names might be fictional propagated over network not real theoretists
