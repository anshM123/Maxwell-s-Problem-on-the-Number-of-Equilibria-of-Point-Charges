# Maxwell's problem for three point charges: statement

## The question

Maxwell (Treatise on Electricity and Magnetism, §113, 1873) asserted that k fixed point charges in ℝ³ create at most
(k − 1)² points of equilibrium of their electrostatic field.
- For k = 5 the assertion is false: Arathoon, Ball and Kvalheim, *The Maxwell conjecture is false*, arXiv:2607.27197
  (29 July 2026). They give five positive charges with at least 24 non-degenerate critical points, against 16.
- For k = 3 it was open (the outcome is at the end of this file). The bound would be 4.
  - Best known bound: 6 for three positive charges, any α > 0 and any dimension (Gabrielov, Novikov, Novikov and Shapiro,
    *From 12 to 6*, arXiv:2607.28785, 30 July 2026). It improves the bound 12 of Gabrielov, Novikov and Shapiro, *Mystery
    of point charges*, Proc. LMS 2007.
  - Equal charges: settled, bound 4 (Tsai, Physica D 309, 2015).
  - Collinear charges: trivial.
  - GNNS write that "it remains possible that some configuration of three charges has exactly six equilibrium points for
    special positive values of α".

## Setting (as in GNNS)

- Charges ζ₁, ζ₂, ζ₃ > 0 sit at non-collinear points p₁, p₂, p₃ of ℝⁿ, n ≥ 2.
- The potential is V_α(x) = Σ ζᵢ |x − pᵢ|^(−2α), with α > 0.
- The Coulomb case of ℝ³ is α = 1/2. More generally, the Newtonian case of ℝⁿ is α = (n − 2)/2.
- Equilibria are the critical points of V_α. The count is of non-degenerate critical points.

## Reduction facts (GNNS §1, reproved in code/checks)

- Every critical point lies in the open triangle p₁p₂p₃. Outside the plane, the normal component of ∇V_α equals h·Σ(…) ≠ 0.
  Outside the triangle, moving away from a supporting line decreases every term.
- In the plane, V_α is strictly subharmonic (Δ₂ |x|^(−2α) = 4α² |x|^(−2α−2) > 0), so it has no local maximum there.
- Each critical point is either a planar local minimum or a planar saddle. Poincaré–Hopf on the sphere (the plane plus
  ∞) with the three charges as maxima and ∞ as the minimum gives:

  #saddles − #minima = 2,   #equilibria = 2·m₀ + 2,   m₀ = number of planar local minima.

- Hence Maxwell's bound 4 for three charges ⟺ m₀ ≤ 1, and the bound 6 ⟺ m₀ ≤ 2.

## Target

**Decide whether three positive charges can have 6 non-degenerate equilibria (m₀ = 2).**
- Coulomb case α = 1/2 first, then all α > 0.
- If yes: give an explicit configuration and certify its six equilibria rigorously (interval Newton / Krawczyk for each
  zero of ∇V in the plane, plus exclusion of further zeros in the rest of the triangle).
- If no: prove m₀ ≤ 1 by a computer-assisted argument over the 4-parameter family (shape of the triangle and two charge
  ratios), with charts at the degenerate corners (collinear limit, colliding charges, vanishing charges).

## Sanity tests

- Equilateral triangle with equal charges: 4 equilibria (1 minimum at the centre, 3 saddles).
- Collinear charges: 2 equilibria.
- Every sample satisfies #saddles − #minima = 2. This completeness check is built into the scan.

## Outcome

- **Positive charges, every α > 0.** m₀ ≤ 1, so there are at most 4 critical points (PROOF.md). The proof uses no
  computer.
- **Charges of different signs (MIXED.md).**
  - Coulomb case: at most 4 critical points for non-collinear charges, and at most 4 non-degenerate ones for collinear
    charges.
  - The same holds for every α ∈ [1/4, 2].
  - One step, the convexity of an explicit planar region, is verified by interval arithmetic.
