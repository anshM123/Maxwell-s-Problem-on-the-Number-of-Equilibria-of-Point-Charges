# Three point charges of mixed signs: at most four equilibria (Coulomb, and every α ∈ [1/4, 2])

**Status (2026-10-24).** The proof is complete. One step, Claim M (§6), is computer-assisted: an interval-arithmetic
certificate, code/mixed_cert_all.py, finishes in seconds with 0 failed boxes.

**Internal audit (audit/AUDIT_2.md): no gap in the logic.**
- Every identity was checked at 40 digits.
- Lemma C was re-derived symbolically.
- Every Coulomb box was re-certified in rigorous 80-bit `mpmath.iv` arithmetic, with an independent implementation of
  the formulas, so the libm assumption is removed.
- Coverage issue (minor): a 1e-16 coverage sliver at A₁ + A₂ = π. Now fixed: the driver's upper edge is up(π), and it
  re-certifies with 0 failures.
- Notes N1–N3 (wording) are fixed below.
- Counterexample search: 1,200 forced-minimum configurations all give exactly (1, 3), and winding numbers are ≤ 1 at
  941 shapes.

**Theorem 3 (every α ∈ [1/4, 2], §8) and its audit (audit/AUDIT_3.md): no gap.**
- The reduction was re-derived for general α.
- The interval-F enclosures were tested against an independent reference: 0 failures at 7,164 points.
- Eight slices were re-run, with identical box counts. The coverage was re-checked in exact arithmetic.
- A direct search at α = 1/4, 1 and 2 found nothing.
- Its ten notes are incorporated. In particular, ia._ratio now propagates NaN (N6). All 163 band charts of the sweep were
  re-run with the fixed code, with counts identical to the logs (cert/band_recheck.out).
- The overclaim it found (N9) concerns collinear charges and is now stated precisely (§8).

These are internal checks, not referee reports.

## Theorem 2

Let V = Σ qᵢ |x − pᵢ|^(−1) be the Coulomb potential of three point charges q₁, q₂, q₃ ≠ 0, not all of the same sign,
at distinct points of ℝ³.

- (a) If the charges are not collinear, V has **at most four** critical points. All of them lie in the plane of the
  charges. The same holds in ℝⁿ for every n ≥ 2.
- (b) If the charges are collinear, the off-axis equilibria come in circles about the axis (so they are not
  isolated), and V has at most four non-degenerate equilibria.

Combined with Theorem 1 (PROOF.md: three positive charges, every α > 0), this gives **Maxwell's bound (k − 1)² = 4
for three point charges in ℝ³**. It holds for all critical points when the charges are not collinear, and for
non-degenerate equilibria when they are collinear.

Notation:
- α = ½, β = 3, b = 2, c = b/(b + 1)² = 2/9, K = (b + 1)²/b = 9/2, F = K/2 − 1 = 5/4.
- Replacing V by −V and relabelling reduces every mixed sign pattern to **q₁, q₂ > 0 > q₃**.
- The triangle p₁p₂p₃ is labelled counter-clockwise, with angles Aᵢ and area Δ > 0.
- R and O are the circumradius and circumcentre, and P_x = R² − |x − O|².

## 1. Planarity (non-collinear charges)

Write x = y + h n with y ∈ Π and n ⊥ Π, and put wᵢ = qᵢrᵢ^(−β) ≠ 0. Then

  ∇V = −2α(Σwᵢ(y − pᵢ) + h n Σwᵢ).

If h ≠ 0, then Σwᵢ = 0 and Σwᵢpᵢ = 0, so w = 0 because the pᵢ are affinely independent. That is impossible. Critical
points therefore lie in Π, and are the critical points of V|Π.

## 2. Location and type

1. x ∈ Π ∖ {pᵢ} is critical iff wᵢ = κλᵢ(x) for some κ ≠ 0, where λ are the signed barycentric coordinates.
2. For the sign pattern (+, +, −) this means:
   - x ∈ R₊ = {λ₁, λ₂ > 0 > λ₃} (κ > 0): beyond side p₁p₂, inside the angle at p₃; or
   - x ∈ R₋ = {λ₁, λ₂ < 0 < λ₃} (κ < 0): the vertical angle at p₃.
3. Hess V = 2ακ(βS − I), with S = Σλᵢuᵢuᵢᵀ and tr S = Σλᵢ = 1. As in PROOF.md Lemma 1, x is an extremum iff
   D = det S > c, a minimum if κ > 0 and a maximum if κ < 0. It is a non-degenerate saddle iff D < c.

**Lemma D (all x off the side lines).** D(x) = P_x sin σ₁ sin σ₂ sin σ₃ / (2Δ).
- Here σᵢ = arg((pₖ − x)/(pⱼ − x)) ∈ (−π, π) are signed angles, for cyclic (i, j, k).
- Proof: λᵢ = rⱼrₖ sin σᵢ/(2Δ) and, by Cauchy–Binet, D = Σₖ λᵢλⱼ sin²σₖ. Hence
  D = (r₁r₂r₃ Π sin σ/(2Δ)²) Σₖ rₖ sin σₖ.
- Also r₁r₂r₃ Σₖ rₖ sin σₖ = Σₖ rₖ² · 2Δλₖ = 2Δ Σₖ λₖ|x − pₖ|² = 2Δ P_x, by the Leibniz identity for affine weights.

**No maxima.**
- In R₋, sign(sin σᵢ) = sign(λᵢ), so Π sin σ > 0.
- Also P_x < 0 there: p₃ is a convex combination of x, p₁ and p₂ with positive weights. If x lay in the closed
  circumdisc, p₃ would lie in the open disc, but it is on the circle.
- So D < 0 < c, and every critical point in R₋ is a non-degenerate saddle.

**Minima lie outside the circumcircle.** In R₊, Π sin σ < 0, so D > 0 forces P_x < 0. Minima therefore lie in
R₊ᵒᵘᵗ = R₊ ∩ {P_x < 0}.

## 3. Angle coordinates, separated charges, gradient structure

For x ∈ R₊, segment p₃x crosses p₁p₂, so σ₁ = t₁ ∈ (0, π), σ₂ = t₂ ∈ (0, π) and σ₃ = −s, where s = t₁ + t₂.

**Lemma P (signed pedal identity).** rⱼrₖ sin(σᵢ − Aᵢ) = P_x sin Aᵢ.
- The proof is PROOF.md Lemma 3 with signed angles: inversion at x, the law of sines in the image triangle, and the
  image of the circumcircle, of radius R/|P_x|.
- Checked numerically, relative error ≤ 1e-9.
- On R₊ᵒᵘᵗ, sin(tᵢ − Aᵢ) < 0, so **tᵢ ∈ (0, Aᵢ)**. Put Q = (0, A₁) × (0, A₂).

**Injectivity of x ↦ (t₁, t₂) on R₊ᵒᵘᵗ.** x lies on the circle through p₂, p₃ fixed by t₁ and on the circle through
p₃, p₁ fixed by t₂. These circles are distinct: their common circle would be the circumcircle, and x is outside it. So
x is their second intersection point.

**Separated charges.** Lemma P gives rᵢ = r₁r₂r₃ |sin(σᵢ − Aᵢ)| / (|P_x| sin Aᵢ), so

  λᵢrᵢ^β = C(x) sin σᵢ |sin(σᵢ − Aᵢ)|^b / sin^b Aᵢ,   C(x) = (r₁r₂r₃)^(b+1)/(2Δ|P_x|^b) > 0.

With φᵢ(t) = log sin t + b log sin(Aᵢ − t) (i = 1, 2), φ₃(s) = log sin s + b log sin(s + A₃), and
H(t) = Φ₁(t₁) + Φ₂(t₂) − Φ₃(t₁ + t₂) (Φᵢ′ = φᵢ):

- x ∈ R₊ᵒᵘᵗ is critical for q iff ∇H(t(x)) = ψ_q;
- ψ_q = (log(q₁/|q₃|) + b log(sin A₁/sin A₃), log(q₂/|q₃|) + b log(sin A₂/sin A₃)).

## 4. The minimum region and the sign of Hess H

Put Uᵢ(t) = cot t + cot(aᵢ − t), with a₁ = A₁, a₂ = A₂, a₃ = A₁ + A₂ (U₃ evaluated at s).
- With Wᵢ = cot(σᵢ − Aᵢ) − cot σᵢ, we get W = (−U₁, −U₂, U₃), and e₂(W) = **E := U₁U₂ − U₃(U₁ + U₂) = 1/D**.
  This is PROOF.md Lemma 4 with signed quantities; checked numerically.
- So the minima, together with the degenerate points D = c, lie in **M̄′**, the closure of M′ = {t ∈ Q : E(t) < K}.
- **Hyperbolic form.** x̃ = (cot t₁, cot t₂, −cot s) and ỹ = (cot(A₁ − t₁), cot(A₂ − t₂), −cot(A₁ + A₂ − s)) satisfy
  e₂ = 1, and x̃ + ỹ = (U₁, U₂, −U₃). So E = 2 + 2⟨x̃, ỹ⟩.
  - x̃ and ỹ are the cot-vectors of the genuine triangles (t₁, t₂, π − s) and (A₁ − t₁, A₂ − t₂, A₃ + s), so they lie
    on the same sheet {e₂ = 1, Σ > 0}.
  - Hence E ≥ 4, with equality iff t = t° = (A₁/2, A₂/2), the excentre opposite p₃.
  - Only E(t°) = 4 < K is used (AUDIT_2, N1).
- **Definiteness.** d = (φ₁′, φ₂′, −φ₃′) = x̃ − bỹ, so e₂(d) = (b + 1)² − bE > 0 on M′. Thus Hess H (Gram determinant
  e₂(d)) is definite on M′.
- **Sign.** At t°, d = ((1 − b) cot(A₁/2), (1 − b) cot(A₂/2), (b − 1) cot((A₁ + A₂)/2)), so the form is negative
  definite. Since M′ is connected (§6), H is strictly concave on M′.

## 5. Counting (planar index)

- Around each charge the degree is 1, for either sign.
- At infinity the degree is 1 if Q_tot = Σqᵢ ≠ 0. If Q_tot = 0 it is 2: the dipole term, with moment
  m = Σqᵢpᵢ ≠ 0, gives a field ∝ (1 − γ/2) − (γ/2)e^{2iφ} on large circles, with γ = 2α + 2 > 2.
- Extrema have index +1 and saddles −1. If V|Π is Morse, #saddles = #extrema + 2 (or + 1 when Q_tot = 0).
- Degeneracy: as in PROOF.md §7 Step 3, Hess V_q = 2ακ L∘DΨ with L invertible (x is not on a side line). Sard gives
  dense Morse charges, and saddles persist. For mixed signs, use log|Qᵢ/Q₃|, take Ψ on R₊ ∪ R₋, and allow κ of either
  sign.

## 6. Claim M: M′ is strictly convex (Coulomb)

### 6.1 Explicit description of M′

**Hyperbolic arches.** Put Λᵢ(t) = log(sin t / sin(aᵢ − t)). Then Λᵢ′ = Uᵢ, Uᵢ = 2(cosh Λᵢ + cos aᵢ)/sin aᵢ,
Uᵢ′ = 2 sinh Λᵢ Uᵢ / sin aᵢ, and t = arg(1 + e^{Λ + i a}).

**Product representation** (strip model of the hyperbolic plane, anchored at vertex 3). With
Δ = Λ₁(t₁) − Λ₂(t₂),

  E = (2U₃(s)/sin A₃)(cosh Δ + cos A₃).

Hence, since 1/U₃ = sin s sin(s + A₃)/sin A₃,

  E < K  ⟺  cosh Δ < β(s) := [(F − 1) cos A₃ − (F + 1) cos(2s + A₃)]/2.

The analogous representations at vertices 1 and 2 also hold (all checked to 1e-11):

  E = (2U₁(t₁)/sin A₁)(cosh(Λ₂ − Λ₃) − cos A₁),

and likewise for vertex 2.

**(s, Δ) coordinates.**
- The map t ↦ (s, Δ) is a diffeomorphism of Q onto (0, A₁ + A₂) × ℝ. Its Jacobian has det = −(U₁ + U₂) ≠ 0. For
  fixed s, Δ increases strictly from −∞ to +∞ along t₁ + t₂ = s.
- Put φ = 2s + A₃ = π + ψ. Then β = [(F − 1) cos A₃ + (F + 1) cos ψ]/2, and β ≥ 1 iff |ψ| ≤ ψ_max, where
  cos ψ_max = (2 − (F − 1) cos A₃)/(F + 1).
- ψ_max < A₁ + A₂, so the s-range [A₁₂/2 − ψ_max/2, A₁₂/2 + ψ_max/2] lies inside (0, A₁₂), with A₁₂ = A₁ + A₂.
- So M′ = {|Δ| < B(s)}, B = arccosh β, is a topological disc. Its boundary {cosh Δ = β(s)} is a smooth simple closed
  curve: ∇(cosh Δ − β) = (−β′, sinh Δ) ≠ 0, because β′ ≠ 0 where Δ = 0.
- On M′, for Coulomb (F = 5/4): P := −β″ = 2(F + 1) cos ψ ≥ 3.5, and all |Λᵢ| ≤ log 2 and |Δ| ≤ log 2. These bounds
  depend on F; the certificate does not use them.

### 6.2 The curvature inequality

**Lemma C (bordered Hessian).** Take coordinates X with Ĝ(X) = β(X₁ + X₂) − cosh(Λ̃₁(X₁) − Λ̃₂(X₂)), and write
Ũᵢ = Λ̃ᵢ′ (any sign). Let S = sinh Δ, C = cosh Δ = β, P = −β″, and take the tangent
τ = (β′ + SŨ₂, SŨ₁ − β′) ≠ 0. Then

  −Hess Ĝ(τ, τ) = J := (Ũ₁ + Ũ₂)²(P S² + C β′²) + S[Ũ₁′(β′ + SŨ₂)² − Ũ₂′(SŨ₁ − β′)²].

*Proof.*
1. ∇Δ = (Ũ₁, −Ũ₂) and Hess Δ = diag(Ũ₁′, −Ũ₂′), so Hess Ĝ = β″JJᵀ − C∇Δ∇Δᵀ − S Hess Δ.
2. τ ⊥ ∇Ĝ = β′(1, 1) − S∇Δ.
3. τ₁ + τ₂ = S(Ũ₁ + Ũ₂) and Ũ₁τ₁ − Ũ₂τ₂ = β′(Ũ₁ + Ũ₂). Substituting gives the formula. ∎

The boundary of {Ĝ > 0} has positive curvature at a point iff J > 0 there. In particular {Ĝ > 0} is strictly convex
near such a point.

- **Vertex 3** (X = t, Ũᵢ = Uᵢ > 0). Dividing by (U₁U₂)² gives

    I = (c₁ + c₂)²(P S² + C β′²) + S[a₁(β′c₂ + S)² − a₂(S − β′c₁)²],
    cᵢ = 1/Uᵢ = sin Aᵢ/(2hᵢ),   aᵢ = Uᵢ′/Uᵢ² = sinh Λᵢ / hᵢ,   hᵢ = cosh Λᵢ + cos Aᵢ.

  - I is invariant under the mirror swap 1 ↔ 2, S → −S. So the lower arc of a shape is the upper arc of the mirrored
    shape.
  - Validation:
    - at α = ½ the sign of I matches the bordered Hessian of E (computed directly from E) at 72,000 boundary points;
    - at α = 10, where M′ is genuinely non-convex in places, I < 0 exactly where that bordered Hessian is negative
      (0 disagreements).
- **Vertex 1** (X = (s, −t₂), Λ̃₁ = Λ₃, Λ̃₂(X₂) = Λ₂(−X₂), so Ũ₁ = U₃ and Ũ₂ = −U₂; β₁(t₁) = [(F + 1) cos ψ₁ −
  (F − 1) cos A₁]/2 with ψ₁ = 2t₁ − A₁). Dividing by U₂²:

    J₁/U₂² = (1 − U₃c₂)²(P S² + C β₁′²) + S[U₃′(β₁′c₂ − S)² − a₂(S U₃ − β₁′)²],
    U₃′ = −U₃(cot s + cot(s + A₃)).

  - Validation: its points lie on {E = K} to 2e-10, and at α = 10 there are 0 sign disagreements in 3,943 points.

**Boundary points in closed form.**
- Vertex-3 chart: given (s, Δ = ±B), ζ = e^{Λ₁} is the positive root of

    e^{−Δ} sin(A₁₂ − s) ζ² + [sin(A₁ − s) + e^{−Δ} sin(A₂ − s)] ζ − sin s = 0.

  This comes from arg(1 + ζe^{iA₁}) + arg(1 + ζe^{−Δ}e^{iA₂}) = s.
- Vertex-1 chart: ζ = e^{Λ₂} is the positive root of

    e^{Δ} sin(A₁ − t₁) ζ² + [e^{Δ} sin(A₁₂ − t₁) − sin(A₂ + t₁)] ζ − sin t₁ = 0.

  Then t₂ = arg(1 + ζe^{iA₂}). In the corner chart (A₂ ≤ r₀ < π/2) this equals arctan(ζ sin A₂ / (1 + ζ cos A₂)). The
  arctan form takes the wrong branch when 1 + ζ cos A₂ < 0, which needs A₂ > π/2 (AUDIT_2, N2).

### 6.3 The certificate (code/ia.py, code/mixed_cert.py, code/mixed_charts.py, code/mixed_cert_all.py)

Shape space is covered by three charts. The boundary point is parametrised by ρ ∈ [−1, 1] with ψ = ρψ_max.

| chart | shapes | function certified > 0 | notes |
|---|---|---|---|
| main | A₁ + A₂ ∈ [e₀, π], A₁, A₂ ≤ π − r₀ | I (vertex 3, upper arc) | lower arc: mirror |
| band | A₁ + A₂ = e ∈ [0, e₀] | I/e⁴ in scaled variables | regular at e = 0 |
| corner | A₁ = π − r ≥ π − r₀, A₂ = θr, A₃ = (1 − θ)r | J₁/U₂² (vertex 1), both arcs | A₂ ≥ π − r₀: mirror |

**Band-chart scaling.**
- ψ_max = 2 arcsin(q sin(e/2)) with q = √((F − 1)/(F + 1)), which is 1/3 for Coulomb.
- β − 1 = (F + 1) sin(ψ_max(1 + ρ)/2) sin(ψ_max(1 − ρ)/2) = O(e²).
- S, β′, cᵢ and Δ are each e times an analytic factor, written with sinc, arcsin(x)/x and arcsinh(x)/x.
- (a₁ − a₂)/e is enclosed by the mean value theorem (see mixed_charts.py).

**Arithmetic.**
- IEEE binary64, with every operation rounded outward by 2 ulps and every library transcendental by 4 ulps.
- This assumes libm is accurate to within 2 ulps, the same assumption as Problem 8 (astro/orbit_distance_P8/FLOAT_RIGOUR.md).
- ia.py passes containment tests against 50-digit mpmath for every operation.
- The interval I contains the 40-digit value at random points (0 failures in 300); the corner chart agrees with the
  float formula at 400 points.

**Results** (python mixed_cert_all.py e0 r0 rounds):

| e₀ | r₀ | main boxes | band | corner ± | failed |
|---|---|---|---|---|---|
| 0.2 | 0.05 | 127,403 | 2,048 | 1,024 + 1,024 | **0** |
| 0.1 | 0.1 | 2,027,388 | 2,048 | 1,024 + 1,024 | **0** |
| 0.4 | 0.02 | 21,711 | 2,048 | 1,024 + 1,024 | **0** |

The main chart alone does not reach down to e₀ = 0.05, which is why the band chart exists. Every run with
e₀ ∈ [0.1, 0.4] closes.

- The main chart's A₁ + A₂ range ends at up(π), the double just above π, so triangles with arbitrarily small A₃ > 0
  are covered. This was the AUDIT_2 coverage fix; all three runs above were redone with it, still with 0 failures.
- Independent rigorous re-certification (audit/b4e_mpiv_recert.py, b4f_mpiv_band_corner.py): `mpmath.iv` at
  80 bits, with a separate implementation of the formulas. All Coulomb boxes certify, with no libm assumption.

**Conclusion of §6.** For every triangle, ∂M′ is a smooth simple closed curve with strictly positive curvature, by
Lemma C on both arcs and the mirror symmetry. So M̄′ is strictly convex.

## 7. Proof of Theorem 2

**(a) Non-collinear charges.** Reduce to q₁, q₂ > 0 > q₃. By §1 all critical points lie in Π.
1. Critical points in R₋ are saddles (§2).
2. Critical points with D ≥ c lie in R₊ᵒᵘᵗ (§2) and correspond to points t ∈ M̄′ with ∇H(t) = ψ_q. The map x ↦ t is
   injective (§3).
3. Suppose two such points t ≠ t′ existed. The open segment between them lies in M′ (strict convexity, §6), where
   Hess H ≺ 0 (§4). Then ⟨∇H(t) − ∇H(t′), t − t′⟩ < 0, a contradiction. So there is at most one minimum or degenerate
   point.
4. Saddles: if there were four, a nearby Morse q′ would keep all four (§5). But q′ has at most one extremum, so at most
   1 + 2 = 3 saddles.
5. Total ≤ 4. ∎

**(b) Collinear charges.**
1. Off the axis, ∇V = 0 forces Σwᵢ = 0 and Σwᵢtᵢ = 0. Rotational symmetry about the axis then gives circles of
   equilibria, which are not isolated.
2. A non-degenerate on-axis equilibrium is non-degenerate in a plane Π through the axis: Hess = diag(f″, h, h).
3. Move p₃ slightly inside Π, off the axis. Every such equilibrium survives as a distinct critical point of the
   perturbed, now non-collinear, configuration (implicit function theorem).
4. Part (a) then bounds their number by 4. ∎

Numerically, collinear mixed configurations have 0 or 2 axis equilibria (40,000 samples, code/collinear_mixed.py).

Part (b) is for ℝⁿ with n ≥ 3. In ℝ², the off-axis equilibria of collinear charges are isolated mirror pairs, and (b)
does not cover them (AUDIT_2, N3).

Part (b) counts only non-degenerate equilibria, and isolated degenerate ones occur (AUDIT_3, N9). For charges
(−2, 24, −64) at 0, 1, 3 on a line, the point −1 is critical with zero Hessian. No circle of critical points comes near
it (audit/c12_collinear_degenerate.py). So for collinear charges, Maxwell's bound is proved for non-degenerate equilibria
only.

## 8. Remarks

- **Large α.** Claim M fails for large α. At α = 10, M′ is not convex at the shape (10.6°, 29.3°, 140.1°): I < 0 and the
  bordered Hessian is −6e-3. Yet the winding number of ψ(∂M′) is still 1 there and at 60 random shapes
  (code/mixed_winding.py), so injectivity, and Maxwell's bound, seem to persist. Proving that needs another argument.
- **Where convexity stops** (code/convexity_threshold.py, output cert/convexity_threshold.out; numerical, not rigorous).
  The scale-free minimum of I/(A₁ + A₂)⁴ over triangles and boundary points:

  | α | 1/4 | 2 | 4 | 5 | 5.6 | 5.7 | 5.8 | 6 | 10 |
  |---|---|---|---|---|---|---|---|---|---|
  | min | +0.0066 | +0.073 | +0.076 | +0.040 | +0.011 | +0.0053 | −0.00035 | −0.012 | −0.37 |

  - Beyond α ≈ 5, the minimum sits at nearly degenerate triangles (one angle near 180°, one near 0°).
  - Every value was confirmed with the independent reference at 40–50 digits.
  - So the convexity threshold is about α ≈ 5.75.
  - Between α = 2 and about 5.7, the certified range is limited by computing time rather than by a failure of
    convexity. Near 5.75 the margin tends to 0.
  - As α → 0, K decreases to 4 = min E, so M′ shrinks to the point t° and the charts degenerate. Small α would need
    charts rescaled around t°.
- **Theorem 3 (other exponents, certified 2026-10-24).** Theorem 2 holds for V = Σ qᵢ|x − pᵢ|^(−2α) for every
  α ∈ [1/4, 2]: part (a) in ℝⁿ for n ≥ 2, and part (b) for n ≥ 3, where the off-axis critical points form
  (n − 2)-spheres about the axis. Only Claim M depends on α, and it depends on α only through
  F = (b+1)²/(2b) − 1. F increases from F(1/4) = 13/12 to F(2) = 13/5.
  - code/mixed_certF.py runs the three charts with F as a fourth interval variable, for one α-slice at a time.
  - code/mixed_alpha_sweep.py ran the slices with (e₀, r₀) = (0.4, 0.05). There are 163 certified slices, of widths
    0.05 down to 2⁻⁹, and about 1.3·10⁹ boxes, roughly 4 hours of CPU time. The outputs are cert/alpha_sweep*.out and
    .log; cert/runs.txt lists them.
  - code/mixed_certF_sliver.py covers the main chart for A₁ + A₂ ∈ [π − 10⁻⁶, up(π)] and all F, with output
    cert/mixed_certF_sliver.out. The slices' main chart stops at math.pi < π.
  - code/sweep_coverage.py recomputes each certified slice's F-box. Their union is [1.0833…, 2.6000…], with no gaps. It
    contains [13/12, 13/5], with the ends compared as exact rationals (cert/sweep_coverage.out).
  - code/check_certF_containment.py tests the interval-F enclosures against the independent reference audit/b0_common
    (boundary point by bisection, 40–60 digits). There were 0 failures at 3,965 points (cert/check_certF_containment.out).
  - The α-range run was not repeated in mpmath.iv, so it relies on the libm accuracy assumption of ia.py.
- **Forced-minimum tests** (code/mixed_predict.py). At α = ½, all 300 configurations give exactly (1 minimum,
  0 maxima, 3 saddles).
