# Critical Identifiability Theory (CIT)

**A research programme on what quantum systems can actually identify—not just detect.**

Physical sensitivity can diverge at criticality while information about the underlying physical cause remains zero after accounting for nuisance parameters. CIT studies this distinction using operational quantum response maps, nuisance-projected information geometry, critical finite-size scaling, and nonperturbative asymptotics.

**Status:** working mathematical research notes and computational certificates, **not** a peer-reviewed novelty or experimental-gravity claim. The underlying Ising model, edge modes, Fisher-information projections and resurgence tools have established prior literature; the focus here is their combination in operationally specified identifiability questions.

## Reproduce the open-chain certificate

```bash
python -m pip install numpy scipy
python certificates/ising_tunnelling_certificate.py
```

The script performs independent finite-spin-sector and free-fermion checks, tests the critical-window scaling and boundary quantization conditions, and compares the boundary spectral integral with the Borel-cut/Bessel identity. Assertions are **numerical tests, not computer-assisted proofs of every statement below**.

## A microscopic open Ising testbed

Consider

```text
H_N(β,g) = -J Σ_(j=1)^(N-1) Z_j Z_(j+1)
           -Jg [ β X_1 + Σ_(j=2)^N X_j ]
```

with `J>0`, `g>0`, open boundaries, parity `Π=∏_j X_j`, and a boundary transverse coupling `β`. In the Jordan–Wigner representation, positive excitation energies equal `2J` times singular values of a lower-bidiagonal matrix with diagonal `(βg,g,…,g)` and subdiagonal `-1`.

At `β=0`, a parity doublet is degenerate. The *signed* gap slope with respect to `β` is

```text
c_N(g) = (2J/ℏ) g^N / sqrt(Σ_(k=0)^(N-1) g^(2k)).
```

For fixed `0<g<1`, this decays exponentially in `N`; with `S(g)=-log(g)` the action is `e^(-NS)`. The expansion of that tunnelling overlap in `e^(-NS)` **converges**—it does not itself establish factorial divergence or Stokes phenomena.

### Operational nuisance projection

A parity-referenced Ramsey protocol with two tagged, sign-reversed settings has local phase

```text
φ_s = T (η + s c_N β),   s ∈ {-1,+1}.
```

Against an unknown setting-*even* frequency offset `η`, the identifiable classical Fisher information for `β` is

```text
F_(β|η) = w c_N²,   w = V_N(T)² T².
```

An unconstrained setting-*odd* nuisance eliminates that information. With independent calibration information `C` for that nuisance, it becomes

```text
F_(β|η,κ) = c_N² w C/(w+C).
```

The phase reference, controlled reversal, preparation, and measurement access are part of the protocol—not freely available by assumption.

### Critical crossover

With `g=e^(-x/N)` and `x≥0` held fixed,

```text
lim_(N→∞) N [ℏ c_N/(2J)]² = 2x/(e^(2x)-1).
```

At the critical point `g=1`, `c_N=(2J/ℏ)/sqrt(N)`. Thus fixed-time information in this specified parity channel changes from exponentially to algebraically suppressed as the edge localization length diverges. It is **not** a general precision advantage: coherent estimation of the local boundary coupling obeys `F_(ββ)≤(2JgT/ℏ)²`, which a directly interrogated single spin can attain under the same idealized resource constraints.

### Genuine resurgence—different observable, different information

For the thermodynamic-limit critical boundary Majorana spectral correlation, with `z=4Jτ/ℏ`,

```text
C(z) = (4/π) ∫_0^1 sqrt(1-t²) exp(-zt) dt.
```

Its large-`z` endpoint expansion is factorially divergent. Its Borel transform is `(4/π)sqrt(1-t²)`, with branch points at `t=±1`. In the stated upper-minus-lower convention the lateral continuation discontinuity is

```text
S_+C - S_-C = -(8i/π) K_1(z)/z.
```

But this two-insertion spectral observable is invariant under the transformation that reverses the sign of the parity-referenced tunnelling amplitude. **A valid Borel–Stokes fingerprint is not necessarily a fingerprint of the physical distinction one wants to identify.**

The tunnelling action `N|log g|` and the spectral band-edge action `4Jτ/ℏ` arise from different asymptotic variables. This does **not** establish a novel gravitational transseries.

## Broader mathematical direction

For an observable tangent map `D`, nuisance columns `B`, and appropriate accessible-output information metric `G`, nuisance-adjusted information is the minimum residual in the quotient by the observable nuisance space:

```text
F_id(u|B) = min_b ||u-Bb||_G².
```

Under a finite-dimensional scaling filtration with controlled metric and fixed response directions, a weighted-minor argument yields a quotient scaling exponent from the difference of maximal weighted minors. More general monomial weights—exponential, power, and logarithmic—can be treated through Gram determinants, provided the assumptions on metric comparability and asymptotic minors hold.

**Research questions:** What happens when nuisance spaces rotate with size? Which critical exponents survive the experimentally accessible quotient? When does a nonperturbative sector remain distinguishable under realistic preparations, noise and calibration? Can any such mechanism outperform an optimally resourced noncritical probe?

## Verification and limits

The certificate exercises selected exact formulas using finite-size numerical spectra and numerical integration. Independent review, stronger precision/cross-platform tests, a full bibliography, and assessment of publication priority remain open. In particular, the derivations **do not demonstrate an experimental quantum-gravity effect**.

Contributions, corrections, and counterexamples are welcome.
