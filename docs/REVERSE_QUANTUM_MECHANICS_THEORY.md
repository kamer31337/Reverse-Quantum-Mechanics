# Reverse Quantum Mechanics: Theoretical Foundations and Mathematical Analysis

**Based on the paper:**
*Reverse Quantum Mechanics*
Gabriele Carcassi, Tobias Thrien, and Christine A. Aidala
Department of Physics, University of Michigan, Ann Arbor (arXiv:2608.27543v1 [quant-ph], 27 Aug 2026)
Part of the open research program: *Assumptions of Physics*

---

## 1. Introduction to Reverse Physics

Physical theories are conventionally formulated either through **axiomatic reconstructions** or **philosophical interpretations**. Both approaches presuppose that standard mathematical formalisms (e.g., Hilbert spaces $\mathcal{H}$) already map one-to-one onto physical reality, reducing foundational physics to choosing starting points or assigning philosophical meaning.

**Reverse Physics** adopts the methodology of **Reverse Mathematics** (Simpson, Stillwell, Friedman):
Instead of postulating axioms and deducing consequences, Reverse Physics **deconstructs physical theories into distinct mathematical and physical conditions** and investigates their **logical dependencies, equivalences, implications, and independences**.

### Core Taxonomical Concepts
- **Base Conditions $[B]$**: Conditions that must hold in *any well-formed physical theory* (e.g., $[B]\text{COBS}$ continuity of measurement, $[B]\text{WDQ}$ frame-invariance, $[B]\text{STDIST}$ physical distinguishability).
- **Physical Conditions**: Statements describing experimentally accessible facts (e.g., reversibility, mutual exclusivity, conservation of information).
- **Mathematical Conditions**: Specific mathematical machinery (e.g., Hilbert spaces, symplectic forms, unitary groups, Lindblad generators).

---

## 2. Classical Mechanics 1-DOF Equivalence (Theorem 2)

Given **Condition CST-1D** (Classical states for 1 degree of freedom: $(q, p) \in \mathbb{R}^2$, Liouville measure $\mu(U) = \int_U dq dp$, probability density $\rho(q, p)$):

### The 8 Equivalent Formulations
1. **HM-1D** (Hamiltonian Mechanics): $\dot{q} = \frac{\partial H}{\partial p}$, $\dot{p} = -\frac{\partial H}{\partial q}$
2. **DI-SYMP** (Symplectomorphism): The flow preserves the symplectic 2-form $\omega = dq \wedge dp$ ($\det J = 1$).
3. **DI-POI** (Poisson Bracket Invariance): $\{f, g\} = \frac{\partial f}{\partial q}\frac{\partial g}{\partial p} - \frac{\partial f}{\partial p}\frac{\partial g}{\partial q}$ is invariant under evolution.
4. **DR-DIV** (Incompressible Flow): The phase space displacement field $S^a = (\dot{q}, \dot{p})$ is divergence-free: $\nabla \cdot S = \partial_q \dot{q} + \partial_p \dot{p} = 0$.
5. **DR-DEN** (Volume Conservation): The Liouville measure $\mu(U)$ is conserved in time.
6. **DR-EV** (Determinism and Reversibility): Past and future states map bijectively ($1$-to-$1$ and onto).
7. **DR-THER** (Thermodynamic Reversibility): Conserved volume entropy $S_{\text{vol}} = \ln \mu(U)$.
8. **DR-INFO** (Information Conservation): Gibbs-Shannon information entropy $S_{\text{info}} = -\int \rho \ln \rho \, dq dp$ is conserved.

$$\text{HM-1D} \iff \text{DI-SYMP} \iff \text{DI-POI} \iff \text{DR-DIV} \iff \text{DR-DEN} \iff \text{DR-EV} \iff \text{DR-THER} \iff \text{DR-INFO}$$

---

## 3. The Unphysicality of Standard Hilbert Spaces (Theorems 3 & 4)

Standard quantum mechanics formalizes physical states as rays in a separable complex Hilbert space $\mathcal{H}$ (**Condition QSTO:H**). Theorem 3 proves that $\text{QSTO:H}$ is **inconsistent** with three necessary physical base conditions:

### 1. Violation of $[B]\text{COBS}$ (Continuity of Measurable Quantities)
- **Base Requirement**: Small state variations must yield small changes in observable expectations.
- **Counterexample**: Harmonic oscillator energy eigenstates $|j\rangle$. Define the state sequence:
  $$|\psi_j\rangle = \sqrt{\frac{j-1}{j}}|0\rangle + \sqrt{\frac{1}{j}}|j\rangle$$
  As $j \to \infty$:
  $$\|\psi_j - |0\rangle\|_{\mathcal{H}} = \sqrt{2 - 2\sqrt{\frac{j-1}{j}}} \longrightarrow 0$$
  Yet the expectation of the number operator $N = a^\dagger a$ satisfies:
  $$\langle \psi_j | N | \psi_j \rangle = \frac{j-1}{j}\cdot 0 + \frac{1}{j}\cdot j = 1 \neq \langle 0 | N | 0 \rangle = 0$$
  In Hilbert space topology, unbounded operators are discontinuous everywhere.

### 2. Violation of $[B]\text{WDQ}$ (Measurable Quantities Well-Defined in All Frames)
- **Base Requirement**: Physical quantities well-defined in one coordinate frame must remain well-defined across valid coordinate transformations.
- **Counterexample**: Consider coordinate transformation $y(x) = \tan\left(\frac{\pi}{2}\operatorname{erf}(x)\right)$. A standard Gaussian wave packet $\psi(x) \in \mathcal{D}(X^2)$ transforms into $\phi(y) = \frac{1}{\sqrt{\pi(y^2+1)}}$.
  The position expectation in frame $Y$ is:
  $$\langle Y^2 \rangle = \int_{-\infty}^\infty \frac{y^2}{\pi(y^2+1)} dy = \infty$$
  The state is outside the domain of observable $Y$, rendering position undefined in frame $Y$.

### 3. Violation of $[B]\text{STDIST}$ (Distinguishability of Different Systems)
- **Base Requirement**: Distinct physical systems must possess distinct mathematical structures.
- **Counterexample**: $L^2(\mathbb{R})$ (1 particle in 1D) and $L^2(\mathbb{R}^{3N})$ ($N$ particles in 3D) are both infinite-dimensional separable Hilbert spaces and are thus **isomorphic** ($\mathcal{H}_1 \cong \mathcal{H}_{3N}$). The inner product alone cannot distinguish a single particle from the entire universe.

### The Resolution: Quantum Vector Space (Condition QSTO:V / Theorem 4)
Definition: A **Quantum Vector Space** $V$ is a complete second-countable complex topological vector space equipped with:
1. An inner product $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{C}$.
2. A set of everywhere-defined self-adjoint defining operators $\{X_i\}_{i\in I}$ satisfying a Lie algebra $[X_i, X_j] = i\hbar \sum_k C_{ij}^k X_k$.
3. The **coarsest vector space topology** making both the defining observables and the inner product continuous.

$$\text{QSTO:V} \implies [B]\text{COBS} \land [B]\text{WDQ} \land [B]\text{STDIST}$$

---

## 4. Pure States vs. Ensembles (Theorems 5 & 6)

### Pure States Are Not Enough (Theorem 5)
- The state space of a classical magnetic dipole is the sphere $S^2$, which has a canonical symplectic 2-form $\omega = \sin\theta d\theta \wedge d\phi$.
- The pure state space of a qubit is $\mathbb{C}P^1 \cong S^2$, which has the identical Fubini-Study symplectic structure.
- Consequently, pure state geometry cannot exclude classical mechanics: $\text{QP-RAY} + \text{CP}$ is consistent with $\neg\text{CQE}$.

### Ensembles Differentiate Classical from Quantum (Theorem 6)
- **Classical Ensemble Space** ($\text{CE}$): Probability distributions over phase space form a **Choquet simplex**. Every mixed state has a **unique** decomposition into pure states.
- **Quantum Ensemble Space** ($\text{QE-DENS}$): Density matrices $\mathcal{D}(\mathcal{H})$ inside the Bloch ball $|\vec{r}| \le 1$ admit **infinitely many equivalent decompositions** into pure states:
  $$\rho = \sum_{i} p_i |\psi_i\rangle\langle\psi_i| = \sum_{j} q_j |\phi_j\rangle\langle\phi_j|$$
- This non-uniqueness of statistical mixtures is what enables quantum superposition and contextuality.

---

## 5. Born Rule, Mutual Exclusivity, and Entropy (Theorem 7)

Over $\text{QE-DENS} + [B]\text{MUTEX} + [B]\text{PMIX}$, the following three conditions are **equivalent**:
1. **BR-BORN**: Transition probability $p(\psi|\phi) = |\langle\phi|\psi\rangle|^2$.
2. **BR-ENT**: Entropy $S(\rho) = -\operatorname{Tr}(\rho \ln \rho)$.
3. **BR-ORME**: Orthogonality represents mutual exclusivity ($\langle\psi|\phi\rangle = 0 \iff p(\psi|\phi) = 0$).

### The Entropic Invertible Bijection
For any two pure states $|\psi\rangle, |\phi\rangle$ with overlap $p = |\langle\phi|\psi\rangle|^2$, their equal mixture $\rho_{\text{mix}} = \frac{1}{2}|\psi\rangle\langle\psi| + \frac{1}{2}|\phi\rangle\langle\phi|$ has eigenvalues:
$$\lambda_\pm = \frac{1 \pm \sqrt{p}}{2}$$
The von Neumann entropy in bits is:
$$S(p) = -\lambda_+\log_2 \lambda_+ - \lambda_-\log_2 \lambda_-$$
- When $p = 0$ (orthogonal), $S = 1$ bit of information.
- When $p = 1$ (identical), $S = 0$ bits.
- $S(p)$ is strictly monotonic and invertible: the Born rule is equivalent to establishing information entropy on mixtures.

---

## 6. Unitary Evolution as Deterministic and Reversible Processes (Theorem 8)

Over $\text{QSTO:H} + [B]\text{EVMIX} + \text{EVTI}$:
$$\text{DR-SCEQ} \iff \text{DR-UNIT} \iff \text{DR-EV} \iff \text{DR-PROB} \iff \text{DR-INFO}$$
- **DR-SCEQ**: $i\hbar \frac{d\psi}{dt} = H\psi$, $\frac{d\rho}{dt} = -\frac{i}{\hbar}[H, \rho]$
- **DR-UNIT**: $\rho(t) = e^{-iHt/\hbar} \rho(0) e^{iHt/\hbar}$
- **DR-EV**: Bijective map of states ($U^\dagger(t)U(t) = I$)
- **DR-PROB**: Conservation of transition probabilities $|\langle\phi(t)|\psi(t)\rangle|^2 = |\langle\phi(0)|\psi(0)\rangle|^2$
- **DR-INFO**: Conservation of von Neumann entropy $S(\rho(t)) = S(\rho(0))$

---

## 7. Projections as Equilibration Processes (Theorems 9 & 10)

### Theorem 9: Projective Measurements as Lindblad Dissipation
A nonselective projective measurement $\Phi_X(\rho) = \sum_x P_x \rho P_x$ associated with observable $X = \sum x P_x$ is equivalent to:
1. **EQ-LIND**: Asymptotic state under purely dissipative Lindblad evolution with $H=0$ and jump operator $X$:
   $$\frac{d\rho}{dt} = -\frac{1}{2}[X, [X, \rho]], \quad \Phi_X(\rho) = \lim_{t\to\infty} \rho(t)$$
2. **EQ-HS**: The closest commuting density operator in Hilbert-Schmidt distance:
   $$\Phi_X(\rho) = \arg\min_{\sigma \in \mathcal{D}(\mathcal{H}), [\sigma, X]=0} \|\rho - \sigma\|_{HS}^2$$
3. **EQ-KL**: The closest commuting density operator in Quantum Relative Entropy:
   $$\Phi_X(\rho) = \arg\min_{\sigma \in \mathcal{D}(\mathcal{H}), [\sigma, X]=0} S(\rho \| \sigma)$$

### Theorem 10: Unitary Evolution as State-Adaptive Projective Limit (Quantum Zeno)
$$\text{DR-UNIT} \iff \text{DR-PROJ}$$
Unitary evolution is the continuous limit ($N \to \infty$) of $N$ infinitesimal state-adaptive projections along the instantaneous trajectory.

---

## 8. Classical Mechanics as the High-Entropy Limit (Theorem 11)

Over $\text{QSTO:H} + \text{1DOF}$:
$$\text{CL-HBAR} \iff \text{CL-HENT} \iff \text{CL-CST}$$
Under an open Lindblad creation channel ($H=0, L=a^\dagger$):
1. Phase space observables stretch: $\hat{X} = X/\sqrt{\lambda}$, $\hat{P} = P/\sqrt{\lambda}$ with $\lambda = e^{\gamma t} \to \infty$.
2. The commutator vanishes: $[\hat{X}, \hat{P}] = \frac{i\hbar}{\lambda} \to 0$.
3. Negative regions of the Wigner quasiprobability distribution $W(q, p)$ smooth out into strictly non-negative classical probability distributions.
4. The Moyal bracket deformation reduces to the Poisson bracket:
   $$\lim_{\hbar \to 0} \{\{f, g\}\}_{MB} = \{f, g\}_{PB}$$

---

## 9. States as Ensembles in Triple Equilibrium (Theorem 12)

Every quantum density operator $\rho$ (pure or mixed) satisfies three simultaneous equilibria:
1. **EE-DYN (Dynamic Equilibrium)**: $\rho$ is stationary under Hamiltonian $H_{\text{dyn}} = \rho$ ($[H_{\text{dyn}}, \rho] = 0$).
2. **EE-MEAS (Spectral Equilibrium)**: $\rho$ is a fixed point of projective measurement $X_{\text{meas}} = \rho$ ($\Phi_X(\rho) = \rho$).
3. **EE-THER (Thermodynamic Equilibrium)**:
   - Full-rank states: $\rho$ is the exact Gibbs state for $H_{\text{ther}} = -\frac{1}{\beta}\ln \rho$.
   - Non-full-rank states (including pure states): $\rho$ is the $k \to \infty$ limit of Gibbs states with kernel penalty $H_k = -\frac{1}{\beta}\ln(\rho|_{\text{supp}}) \oplus \sum_{j \in \text{ker}} k j^2 |j\rangle\langle j|$.

---

## 10. Summary Matrix of Theorems

| Theorem | Focus | Assumptions | Logical Result |
|---|---|---|---|
| **Thm 1** | Mathematical Logic | ZF Set Theory | $\text{AoC} \iff \text{WOP} \iff \text{ZL}$ |
| **Thm 2** | Classical 1-DOF | $\text{CST-1D}$ | $\text{HM-1D} \iff \text{DI-SYMP} \iff \text{DI-POI} \iff \text{DR-DIV} \iff \text{DR-DEN} \iff \text{DR-EV} \iff \text{DR-THER} \iff \text{DR-INFO}$ |
| **Thm 3** | Hilbert Pathology | $\text{QSTO:H}$ | Inconsistent with $[B]\text{COBS}, [B]\text{WDQ}, [B]\text{STDIST}$ |
| **Thm 4** | Quantum Vector Space | $\text{QSTO:V}$ | $\text{QSTO:V} \implies [B]\text{COBS} \land [B]\text{WDQ} \land [B]\text{STDIST}$ |
| **Thm 5** | Pure State Ambiguity | $\text{QP-RAY} + \text{CP}$ | Consistent with $\neg\text{CQE}$ (Symplectic manifolds) |
| **Thm 6** | Ensemble Distinction | $\text{QE-DENS} + \text{CE}$ | $\implies \text{CQE}$ (Simplex vs Multi-decomposition) |
| **Thm 7** | Born Rule & Entropy | $\text{QE-DENS} + [B]\text{MUTEX} + [B]\text{PMIX}$ | $\text{BR-BORN} \iff \text{BR-ENT} \iff \text{BR-ORME}$ (Independent of $\text{QE-DENS}$) |
| **Thm 8** | Unitary Evolution | $\text{QSTO:H} + [B]\text{EVMIX} + \text{EVTI}$ | $\text{DR-SCEQ} \iff \text{DR-UNIT} \iff \text{DR-EV} \iff \text{DR-PROB} \iff \text{DR-INFO}$ |
| **Thm 9** | Measurement as Lindblad | $\text{QSTO:H} (+ [B]\text{FINENT})$ | $\text{EQ-MEAS} \iff \text{EQ-LIND} \iff \text{EQ-HS} (\iff \text{EQ-KL})$ |
| **Thm 10** | Zeno Unitary Limit | $\text{QSTO:H}$ | $\text{DR-UNIT} \iff \text{DR-PROJ}$ |
| **Thm 11** | Classical Limit | $\text{QSTO:H} + \text{1DOF}$ | $\text{CL-HBAR} \iff \text{CL-HENT} \iff \text{CL-CST}$ |
| **Thm 12** | Triple Equilibria | $\text{QSTO:H}$ | $\text{QSTO:H} \implies \text{EE-DYN} \land \text{EE-MEAS} \land \text{EE-THER}$ |
