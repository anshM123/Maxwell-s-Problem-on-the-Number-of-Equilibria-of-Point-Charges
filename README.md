# Maxwell's problem for three point charges: at most four equilibria

**Authors:** Ansh Mishra, Aryan Senthilkumar

## The question

In §113 of his *Treatise on Electricity and Magnetism* (1873), Maxwell asserted that k point charges produce at most
(k − 1)² equilibrium points. For k = 3 that is 4.

Known results, for the potentials Σ qᵢ |x − pᵢ|^(−2α):

- **Gabrielov, Novikov and Shapiro**, *Mystery of point charges*, Proc. LMS 95 (2007): at most 12 non-degenerate
  equilibria for three charges.
- **Gabrielov, Novikov, Novikov and Shapiro**, arXiv:2607.28785 (July 2026): at most 6 for three positive charges.
  They note that 4 was known only for three equal charges (Tsai, Physica D 309, 2015).

For five charges, Maxwell's bound is false (Arathoon, Ball and Kvalheim, arXiv:2607.27197).

## The results

**Theorem 1 (charges of one sign).** Let α > 0, n ≥ 2, q₁, q₂, q₃ > 0, and p₁, p₂, p₃ distinct points of ℝⁿ. Then

  V(x) = q₁|x − p₁|^(−2α) + q₂|x − p₂|^(−2α) + q₃|x − p₃|^(−2α)

has **at most four** critical points.

- If the charges are collinear, there are exactly two critical points.
- Otherwise, at most one critical point is a minimum of V restricted to the plane of the charges, or degenerate. All
  the others are non-degenerate saddles, and there are at most three of them.

The bound is sharp: equal charges at the vertices of an equilateral triangle give one minimum and three saddles.

**Theorem 2 (charges of different signs, Coulomb).** Let q₁, q₂, q₃ be nonzero real numbers, not all of the same sign,
at distinct points p₁, p₂, p₃ of ℝ³, and V(x) = Σ qᵢ |x − pᵢ|^(−1).

- (a) If the pᵢ are not collinear, V has **at most four** critical points, and all of them lie in the plane of the
  charges.
- (b) If the pᵢ are collinear, the critical points off the line form circles about it, and V has at most four
  non-degenerate critical points.

**Corollary (Maxwell's bound for three charges).** Three point charges in ℝ³ with the Coulomb potential (α = 1/2) have
at most four non-degenerate equilibrium points. If the charges are not collinear, they have at most four equilibrium
points altogether.

For collinear charges of different signs, isolated degenerate equilibria can occur, so there the bound is proved for
non-degenerate equilibria only (an example is in the paper, Remark 14.2).

**Theorem 3 (other exponents).** Theorem 2 also holds for V = Σ qᵢ |x − pᵢ|^(−2α) for every α ∈ [1/4, 2]: part (a) in
ℝⁿ for every n ≥ 2, and part (b) in ℝⁿ for n ≥ 3, where the critical points off the line form (n − 2)-spheres about
it. This range contains the Newtonian potentials of ℝ⁴, ℝ⁵ and ℝ⁶ (α = 1, 3/2, 2).

The proof of Theorem 1 is short and does not use the computer. Theorems 2 and 3 have one computer-assisted step, Claim M:
an explicit planar region is convex. It is verified by interval arithmetic.

- [paper/main.pdf](paper/main.pdf) is the paper, *Maxwell's problem for three point charges*; its LaTeX source is
  [paper/main.tex](paper/main.tex). Part I proves Theorem 1 and Part II proves Theorems 2 and 3. In the paper these
  results are Theorems 1.1, 1.2 and 1.4 and Corollary 1.3.
- [PROOF.md](PROOF.md) is Theorem 1 in step-by-step form.
- [MIXED.md](MIXED.md) is Theorem 2 in step-by-step form, with the certificate.
- [STATEMENT.md](STATEMENT.md) has the exact statement and the sanity checks.

## How the proofs work

### Charges of one sign

Let x be a point of the open triangle T = p₁p₂p₃. Let θᵢ be the angle under which x sees the side opposite pᵢ, and
let Aᵢ be the angles of T. Put b = 2α + 1.

1. **Reduction.** x is an equilibrium exactly for the charges qᵢ ∝ λᵢ rᵢ^(b+1), where λᵢ are the barycentric
   coordinates and rᵢ = |x − pᵢ|. It is a minimum in the plane exactly when D(x) = det Σ λᵢ uᵢuᵢᵀ > b/(b+1)², where
   uᵢ = (x − pᵢ)/rᵢ.
2. **Separation of variables.** Inversion at x gives Euler's pedal identity rⱼrₖ sin(θᵢ − Aᵢ) = (R² − |x − O|²) sin Aᵢ.
   Hence
   qᵢ ∝ sin θᵢ · sin^b(θᵢ − Aᵢ) / sin^b Aᵢ.
   The equilibrium equations become the gradient equations of a function of θ = (θ₁, θ₂, θ₃) whose terms each depend on
   one θᵢ.
3. **The region of minima.** D(x) = 1/(W₁W₂ + W₂W₃ + W₃W₁), with Wᵢ = sin Aᵢ/(sin(θᵢ − Aᵢ) sin θᵢ).
4. **Convexity.** Each Wᵢ is a log-convex function of θᵢ, because −log sin is convex. So the region of minima is a convex
   set in the angle coordinates. The potential function is strictly convex on it, and a strictly convex function has
   at most one critical point on a convex set.
5. **Saddles.** A Morse count (#saddles = #minima + 2) and a perturbation argument for degenerate cases give at most
   three saddles.

A geometric reading: x is a minimum for its charges if and only if the hyperbolic distance between the shape of its
pedal triangle (angles θᵢ − Aᵢ) and the shape of the triangle with angles π − θᵢ is less than log(2α + 1).

### Charges of different signs

Replacing V by −V and relabelling gives q₁, q₂ > 0 > q₃.

1. **Location.** Every critical point lies in the plane of the charges. Those in the vertical angle at p₃ are saddles.
   The minima lie beyond the side p₁p₂, outside the circumcircle, around the excentre opposite p₃.
2. **Gradient structure.** In angle coordinates t = (t₁, t₂) the charges separate again: a critical point solves
   ∇H(t) = (fixed vector). Hess H is negative definite on the region M′ = {E < (b+1)²/b} of possible minima, where E is
   explicit.
3. **Convexity of M′.** For one sign this followed from log-convexity. Here one product enters with a minus sign, and
   M′ is not convex for large α. Its boundary has a closed form: M′ = {cosh Δ < β(s)} in the coordinates
   s = t₁ + t₂, Δ = Λ₁(t₁) − Λ₂(t₂). The curvature of the boundary has the sign of an explicit function of the triangle
   and of one boundary parameter.
4. **Claim M.** That function is positive. This is checked by interval arithmetic on three charts that cover all
   triangles. With M′ convex and H concave on it, there is at most one minimum. An index count then gives at most three
   saddles.

## The computer-assisted step

It needs Python 3 and numpy; the symbolic and high-precision checks also need sympy and mpmath.

```
cd code
python mixed_cert_all.py 0.2 0.05            # Claim M for Coulomb: 0 failed boxes, a few seconds
python symbolic_mixed.py                     # the signed identities (sympy)
python symbolic_product.py                   # the product representation of E (sympy)
python mixed_certF.py 0.95 1.05 0.4 0.05     # Claim M for every alpha in [0.95, 1.05], F as an interval variable
python mixed_certF_sliver.py                 # the strip next to A1 + A2 = pi, for all alpha in [1/4, 2]
python check_certF_containment.py 200        # interval-F enclosures vs. an independent 40-60 digit evaluation
python sweep_coverage.py ../cert/runs.txt     # the certified F-boxes of the alpha sweep cover [13/12, 13/5]
python band_recheck.py ../cert/runs.txt       # every band chart of the sweep re-run with the current ia.py
```

The sweep over α ∈ [1/4, 2] uses `mixed_alpha_sweep.py` and takes a few hours of CPU time. Its outputs are in
[cert/](cert/), and [cert/runs.txt](cert/runs.txt) lists the runs. `sweep_coverage.py` checks that the certified F-boxes
cover F ∈ [13/12, 13/5] = [F(1/4), F(2)], where F = (b+1)²/(2b) − 1. Since F is increasing in α, that is all of
α ∈ [1/4, 2].

## Numerical checks (charges of one sign)

The proof of Theorem 1 does not depend on the computer. The scripts in [code/](code/) check every identity numerically,
and compare the theorem's predictions with a direct search for critical points.

```
cd code
python key_identities.py 20000 0.5   # pedal identity, charge formula, D·E = 1, e2(d) = (b+1)^2 − bE, incentre
python separable.py 0.5 60 1.0       # sign(D − c) = sign(e2(d)); no negative-definite region; boundary curvature
python predict_test.py 300 0.5 1     # theorem vs direct critical-point search (exactly 1 min + 3 saddles inside)
python scan.py 2000 0.5              # random charges: never two minima, and #saddles − #minima = 2
```

## Internal audits

The check scripts in audit/ were written separately from code/.

[audit/AUDIT_1.md](audit/AUDIT_1.md) checks PROOF.md (Theorem 1) line by line.
- Every lemma was re-derived, and every identity was checked at 29–41 significant digits.
- Direct searches for two minima found none, for α between 0.05 and 10.
- No gap was found. The wording notes are incorporated.

[audit/AUDIT_2.md](audit/AUDIT_2.md) checks MIXED.md (Theorem 2) and the certificate.
- Every identity was checked at 40 digits. The curvature criterion was re-derived symbolically.
- Every Coulomb box was re-certified in 80-bit `mpmath.iv` arithmetic with separately written formulas.
- A coverage gap of width about 10⁻¹⁶ at A₁ + A₂ = π was found and fixed.
- Searches for counterexamples found none. No gap in the logic was found.

[audit/AUDIT_3.md](audit/AUDIT_3.md) checks the extension to α ∈ [1/4, 2] (Theorem 3) and the paper.
- The reduction was re-derived for general α. Only β > 2, b > 1, K > 4, F > 1 and γ > 1 are used.
- The interval-F enclosures were tested against a third, separately written reference: 0 failures at 7,164 points.
- Eight slices were re-run, with box counts identical to the logs. The coverage of [13/12, 13/5] was re-checked in exact
  arithmetic.
- A direct search at α = 1/4, 1 and 2 found no configuration with two extrema or more than four critical points.
- No gap was found. Its notes are incorporated. In particular, `ia.py` now propagates NaN through `sinc`, `asinc` and
  `asinhc`. These are used only by the band chart, and every band chart of the sweep was re-run with the fixed code
  ([cert/band_recheck.out](cert/band_recheck.out)).

None of this replaces peer review.

## Citation

See [CITATION.cff](CITATION.cff).

## License

MIT; see [LICENSE](LICENSE).
