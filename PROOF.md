# Three positive point charges have at most four equilibria

**Status (2026-10-24).** The proof below is complete. Every identity in it was also checked numerically
(code/key_identities.py, code/separable.py).

**Internal audit (audit/AUDIT_1.md):** no gap found.
- Every lemma was re-derived, and every identity checked at 29–41 digits with independent code.
- No counterexample came up in forced-minimum searches (Coulomb and α = 0.1–10) or in winding-number shape scans.
- Its wording notes (Step 3; §8 remarks) are fixed below.

This is an internal check, not a referee report; outside expert review is still needed.

## Theorem

Let α > 0, n ≥ 2, q₁, q₂, q₃ > 0, and let p₁, p₂, p₃ be distinct points of ℝⁿ. Then

  V(x) = Σᵢ qᵢ |x − pᵢ|^(−2α)

has at most four critical points in ℝⁿ ∖ {p₁, p₂, p₃}.

- If the pᵢ are collinear, V has exactly two critical points.
- Otherwise, at most one critical point is a planar local minimum or degenerate, and at most three are saddles.

The bound is attained, e.g. by equal charges at the vertices of an equilateral triangle (one minimum, three saddles).

**Context.**
- Maxwell (Treatise, §113) asserted the bound (k − 1)² for k charges, which is 4 for k = 3.
- Gabrielov–Novikov–Shapiro proved 12 for three positive charges.
- Gabrielov–Novikov–Novikov–Shapiro (arXiv:2607.28785, 30 July 2026) improved this to 6, and state that 4 is open
  (known only for equal charges, Tsai 2015).
- Arathoon–Ball–Kvalheim (arXiv:2607.27197) showed the bound fails in general, for five charges.

The theorem settles the three-positive-charge case with the sharp bound 4, for every α > 0.

Charges of mixed signs are treated in MIXED.md, with one computer-assisted step. The first observations (§8):
- the equilibria of three non-collinear charges of any signs still lie in their plane;
- there are no planar maxima;
- the minima lie outside the circumcircle, beyond the side opposite the odd charge;
- but the region of minima is not convex in angle coordinates for large α, so the argument does not extend directly.

## 0. Notation and standard facts

Put β = 2α + 2 and b = β − 1 = 2α + 1 > 1, and c = (β − 1)/β² = b/(b + 1)².

1. At a critical point, ∇V(x) = −2α Σᵢ qᵢ rᵢ^(−β)(x − pᵢ) = 0, with rᵢ = |x − pᵢ|. So x = Σ wᵢpᵢ / Σ wᵢ with
   wᵢ = qᵢ rᵢ^(−β) > 0.
   - Collinear pᵢ: every critical point lies on the open segment spanned by the charges. On each of the two gaps
     between consecutive charges, V is strictly convex (each |t − s|^(−2α) is convex on either side of s) and tends to
     +∞ at both ends. Hence there are exactly two critical points.
   - Non-collinear pᵢ (assumed from now on): every critical point lies in the open triangle T = p₁p₂p₃, in the plane Π
     of the charges. For x ∈ Π the gradient lies in Π, so the critical points of V are exactly the critical points of
     V|Π in T. We work in Π = ℝ².
2. For x ∈ T, write:
   - λᵢ for the barycentric coordinates, uᵢ = (x − pᵢ)/rᵢ, and Δ for the area of T;
   - θᵢ = ∠pⱼ x pₖ for the angle at x subtended by the side opposite pᵢ ({i, j, k} = {1, 2, 3});
   - Aᵢ for the angle of T at pᵢ, aᵢ = |pⱼpₖ|, and R, O for the circumradius and circumcentre;
   - P_x = R² − |x − O|² > 0, the power of x with respect to the circumcircle, negated.

   Then θ₁ + θ₂ + θ₃ = 2π, each θᵢ < π, and λᵢ = rⱼrₖ sin θᵢ / (2Δ), because Δλᵢ is the area of x pⱼ pₖ.
3. S(x) = Σ λᵢ uᵢuᵢᵀ has trace 1. By Cauchy–Binet,

     D(x) := det S = Σ_{i<j} λᵢλⱼ (uᵢ × uⱼ)² = Σₖ λᵢλⱼ sin²θₖ   ({i, j, k} = {1, 2, 3}).

## 1. Reduction

**Lemma 1.** Let x ∈ T. Then x is a critical point of V_q if and only if q ∝ Q(x) := (λᵢ rᵢ^β)ᵢ. In that case:
- if D(x) > c, x is a non-degenerate planar local minimum;
- if D(x) < c, x is a non-degenerate saddle;
- if D(x) = c, x is degenerate.

There are no planar local maxima.

*Proof.*
1. By §0.1, x is critical iff (qᵢ rᵢ^(−β))ᵢ ∝ (λᵢ)ᵢ, by uniqueness of barycentric coordinates.
2. In that case Hess V_q(x) = −2α Σ qᵢ rᵢ^(−β)(I − βuᵢuᵢᵀ) = 2ακ(βS − I) for some κ > 0.
3. The eigenvalues s, 1 − s of S (with s ≤ 1/2) give βs − 1 and β(1 − s) − 1 > 0, since β > 2.
4. So the Hessian is positive definite iff s > 1/β. Because t ↦ t(1 − t) is increasing on [0, 1/2], this holds iff
   D = s(1 − s) > (1/β)(1 − 1/β) = c. It is indefinite iff D < c, and singular iff D = c. ∎

**Index count.** Suppose V_q is Morse in Π. Let m₀ be the number of minima and s₀ the number of saddles. On a large
circle, ∇V points toward the triangle (degree 1); on a small circle around each pᵢ, it points toward pᵢ (degree 1).
Poincaré–Hopf on the disc minus three small discs gives m₀ − s₀ = 1 − 3, so s₀ = m₀ + 2.

## 2. Angle coordinates

**Lemma 2.** The map x ↦ θ(x) = (θ₁, θ₂, θ₃) is injective on T.

*Proof.*
1. The points y ∈ T with ∠p₂ y p₃ = θ₁ lie on one circle Γ₁ through p₂, p₃ (inscribed angles; T lies on one side of
   p₂p₃). Likewise, ∠p₃ y p₁ = θ₂ defines a circle Γ₂ through p₃, p₁.
2. Γ₁ ≠ Γ₂: otherwise this circle would pass through p₁, p₂, p₃, i.e. be the circumcircle, which contains no point of
   the open triangle.
3. Two distinct circles meet in at most two points, one of which is p₃ ∉ T. Hence at most one point of T has the given
   (θ₁, θ₂). ∎

## 3. Separation of variables

**Lemma 3 (pedal identity).** For x ∈ T, with δᵢ := θᵢ − Aᵢ:

  (a) δᵢ > 0 and δ₁ + δ₂ + δ₃ = π;
  (b) rⱼ rₖ sin δᵢ = P_x sin Aᵢ;
  (c) λᵢ rᵢ^β = C(x) · sin θᵢ · sin^b δᵢ / sin^b Aᵢ, where C(x) = (r₁r₂r₃)^(b+1) / (2Δ P_x^b) > 0 does not depend on i.

*Proof.* Let ι be the inversion in the unit circle centred at x, and pᵢ* = ι(pᵢ). Recall three standard properties:
- (i) |ι(y) − ι(z)| = |y − z| / (|y − x||z − x|);
- (ii) the triangles x y z and x ι(z) ι(y) are similar, so ∠x ι(y) ι(z) = ∠x z y;
- (iii) a circle with centre O and radius R, not through x, is mapped to a circle of radius R/|R² − |x − O|²|.

Now the steps.
1. ι preserves rays from x, so ∠pⱼ* x pₖ* = θᵢ. These three angles are each < π and sum to 2π, so x lies inside
   T* = p₁*p₂*p₃*.
2. The angle of T* at pᵢ* is therefore
     ∠pⱼ* pᵢ* x + ∠x pᵢ* pₖ* = ∠x pⱼ pᵢ + ∠x pₖ pᵢ        (by (ii))
                            = (Aⱼ − ∠x pⱼ pₖ) + (Aₖ − ∠x pₖ pⱼ) = Aⱼ + Aₖ − (π − θᵢ) = δᵢ.
   This gives (a).
3. The circumcircle of T does not pass through x (x lies strictly inside it), and it is mapped to the circumcircle of
   T*, whose radius is R/P_x by (iii).
4. The law of sines in T* and (i) give aᵢ/(rⱼrₖ) = |pⱼ*pₖ*| = 2(R/P_x) sin δᵢ. Since aᵢ = 2R sin Aᵢ, this is (b).
5. For (c): by (b), rᵢ = r₁r₂r₃/(rⱼrₖ) = r₁r₂r₃ sin δᵢ / (P_x sin Aᵢ), and λᵢ rᵢ = r₁r₂r₃ sin θᵢ / (2Δ).
   Multiply λᵢrᵢ by rᵢ^b. ∎

(b) is equivalent to Euler's theorem on the area of the pedal triangle: the pedal triangle of x has angles δᵢ and
sides rᵢ sin Aᵢ.

Consequently, with fᵢ(t) := log sin t + b log sin(t − Aᵢ) on (Aᵢ, π),

  log(λᵢ rᵢ^β) = log C(x) − b log sin Aᵢ + fᵢ(θᵢ).        (3.1)

## 4. The minimum region

Put Wᵢ := sin Aᵢ / (sin δᵢ sin θᵢ) = cot δᵢ − cot θᵢ, a function of θᵢ alone, and
E := W₁W₂ + W₂W₃ + W₃W₁.

**Lemma 4.** For every x ∈ T,

  D(x) = P_x sin θ₁ sin θ₂ sin θ₃ / (2Δ) = 1/E(θ(x)).

Hence D(x) > c ⟺ E(θ(x)) < (b + 1)²/b.

*Proof.*
1. By Lemma 3(b), sin Aᵢ/sin δᵢ = rⱼrₖ/P_x, so Wᵢ = rⱼrₖ/(P_x sin θᵢ). Hence

     E = (r₁r₂r₃ / (P_x² Π sin θ)) Σₖ rₖ sin θₖ.

2. From §0.2–3,

     D = Σₖ λᵢλⱼ sin²θₖ = (r₁r₂r₃ Π sin θ / (2Δ)²) Σₖ rₖ sin θₖ.

3. Therefore D·E = (r₁r₂r₃ Σₖ rₖ sin θₖ)² / (2Δ P_x)².
4. Since rᵢrⱼ sin θₖ = 2Δλₖ, we have r₁r₂r₃ Σₖ rₖ sin θₖ = 2Δ Σₖ λₖ rₖ².
5. Finally Σₖ λₖ|x − pₖ|² = |x − O|² − 2⟨x − O, Σₖ λₖ(pₖ − O)⟩ + R² = R² − |x − O|² = P_x.
6. So r₁r₂r₃ Σₖ rₖ sin θₖ = 2Δ P_x. Substituting into steps 1–2 gives D = P_x Π sin θ / (2Δ) and D·E = 1. The last
   claim follows from c = b/(b + 1)². ∎

## 5. Convexity

Let H_A = {θ ∈ ℝ³ : θ₁ + θ₂ + θ₃ = 2π, Aᵢ < θᵢ < π}. This is an open convex polygon in the plane L = {Σθᵢ = 2π},
and θ(T) ⊂ H_A by Lemma 3(a). Note that E and the Wᵢ are defined on all of H_A, by the formulas of §4.

**Lemma 5.** E is strictly convex on H_A (along every segment in L). Hence
M̂ := {θ ∈ H_A : E(θ) < (b + 1)²/b} is convex, and the open segment between any two points of
{E ≤ (b + 1)²/b} lies in M̂.

*Proof.*
1. ℓᵢ(t) := log Wᵢ = log sin Aᵢ − log sin(t − Aᵢ) − log sin t satisfies ℓᵢ″ = csc²(t − Aᵢ) + csc² t > 0.
2. For i ≠ j, WᵢWⱼ = exp(ℓᵢ(θᵢ) + ℓⱼ(θⱼ)) is convex, and Hess(e^f) ⪰ e^f Hess f.
3. A nonzero direction v with Σvᵢ = 0 has at least two nonzero entries, say vᵢ, vⱼ. Then
   ℓᵢ″vᵢ² + ℓⱼ″vⱼ² > 0, so E is strictly convex along v.
4. The claims about the sublevel sets follow. ∎

## 6. Positive definiteness

Let dᵢ(θᵢ) := fᵢ′(θᵢ) = cot θᵢ + b cot δᵢ, and write e₂(v) = v₁v₂ + v₂v₃ + v₃v₁.

**Lemma 6.** On H_A, e₂(d) = (b + 1)² − b E. Moreover, for θ ∈ M̂ the quadratic form v ↦ Σ dᵢ(θᵢ) vᵢ² is positive
definite on {Σvᵢ = 0}.

*Proof.*
1. Two cotangent identities:
   - if three angles in (0, π) sum to π or to 2π, their cotangents C satisfy e₂(C) = 1 (from cot(θ₁ + θ₂) = −cot θ₃);
   - so e₂(X) = 1 for X = (cot δᵢ), and e₂(C) = 1 for C = (cot θᵢ).
2. Set B = Σ_{i≠j} XᵢCⱼ. Then W = X − C and d = C + bX, so

     e₂(W) = e₂(X) + e₂(C) − B = 2 − B,
     e₂(d) = e₂(C) + b² e₂(X) + bB = 1 + b² + b(2 − E) = (b + 1)² − bE.

3. In the basis (e₁ − e₃, e₂ − e₃) the form has Gram matrix [[d₁ + d₃, d₃], [d₃, d₂ + d₃]], with determinant e₂(d).
   This is > 0 on M̂, so the form is definite on M̂.
4. M̂ is connected (Lemma 5), so the sign of definiteness is constant on M̂.
5. At θ° = ((π + Aᵢ)/2)ᵢ ∈ H_A (the angles at the incentre):
   - δᵢ = (π − Aᵢ)/2 and Wᵢ = 2 tan(Aᵢ/2);
   - E = 4 Σ tan(Aᵢ/2) tan(Aⱼ/2) = 4 < (b + 1)²/b (because b ≠ 1), so θ° ∈ M̂;
   - dᵢ = (b − 1) tan(Aᵢ/2) > 0, so the form is positive definite there. ∎

## 7. Proof of the theorem (non-collinear case)

Let F(θ) := (f₁(θ₁), f₂(θ₂), f₃(θ₃)) on H_A, and ψ_q := (log qᵢ + b log sin Aᵢ)ᵢ.

**Step 1: correspondence.** By Lemma 1 and (3.1), if x ∈ T is a critical point of V_q, then
F(θ(x)) − ψ_q ∈ ℝ(1, 1, 1). By Lemmas 1 and 4, x is a planar minimum or degenerate iff D(x) ≥ c iff
E(θ(x)) ≤ (b + 1)²/b.

**Step 2: at most one critical point with D ≥ c.**
1. Suppose x ≠ x′ are two such points. Put θ = θ(x) and θ′ = θ(x′). Then θ ≠ θ′ by Lemma 2.
2. F(θ) − F(θ′) = μ(1, 1, 1) for some μ, by Step 1.
3. Let v = θ − θ′ ≠ 0. Then Σvᵢ = 0, so
     0 = μ Σ vᵢ = Σᵢ (fᵢ(θᵢ) − fᵢ(θ′ᵢ)) vᵢ = ∫₀¹ Σᵢ dᵢ(θ′ᵢ + t vᵢ) vᵢ² dt.
4. The integrand is > 0 for every t ∈ (0, 1), by Lemmas 5 and 6, since θ′ + t v ∈ M̂. This is a contradiction.

**Step 3: at most three saddles.**
1. Saddles (D < c) are non-degenerate, hence isolated.
2. The map Ψ(x) = log Q(x) mod ℝ(1, 1, 1) is smooth. A critical point x of V_q is degenerate iff DΨ(x) is singular.
   Indeed, from ∇V_q(x) = −2ακ(x) Σᵢ λᵢ (e^{εᵢ(x)} − 1)(x − pᵢ) with ε = log q − log Q(x) (normalised so that ε = 0
   at the critical point), one gets Hess V_q = 2ακ L ∘ DΨ. Here L(e) = Σ λᵢ eᵢ (x − pᵢ) is an isomorphism
   ℝ³/ℝ(1, 1, 1) → ℝ², because the vectors x − pᵢ are pairwise non-parallel and Σλᵢ(x − pᵢ) = 0.
3. By Sard's theorem, the q′ for which V_q′ is Morse are dense.
4. Suppose V_q had four distinct saddles x₁, …, x₄. These are non-degenerate with D < c. By the implicit function
   theorem applied to Ψ at the xⱼ, every q′ close enough to q has four distinct critical points near them, still with
   D < c, i.e. four saddles. No finiteness assumption is needed.
5. Choose such a q′ that is Morse. By Step 2, V_q′ has at most one minimum, so by the index count of §1 it has
   s₀(q′) = m₀(q′) + 2 ≤ 3 saddles. This is a contradiction, so V_q has at most three saddles.

**Conclusion.** The total number of critical points is ≤ 1 + 3 = 4. ∎

## 8. Remarks

- **Hyperbolic picture.** Let X = (cot δᵢ) and Y = (cot(π − θᵢ)) = −C. These lie on the hyperboloid e₂ = 1, which is
  the cot-model of the space of triangle shapes, isometric to the hyperbolic plane. Then:
  - W = X + Y and 1 + ⟨X, Y⟩ = e₂(W)/2, so ⟨X, Y⟩ = cosh ρ with ρ = ρ(δ, π − θ);
  - D = 1/(2(1 + cosh ρ)), and the minimum region is {ρ < log b} = {log K < log b}, where K is the dilatation of the
    affine map between the pedal triangle (angles θᵢ − Aᵢ) and the triangle with angles π − θᵢ;
  - the classical |Z| (Z = Σ λᵢ uᵢ², D = (1 − |Z|²)/4) equals tanh(ρ/2).
- **The angle map is a bijection T → H_A.** Its inverse is explicit: λᵢ ∝ sin Aᵢ sin θᵢ / sin(θᵢ − Aᵢ).
  - By (§0.2) and Lemma 3(b), every x ∈ T has λᵢ = P_x sin Aᵢ sin θᵢ / (2Δ sin δᵢ).
  - Given θ ∈ H_A, these weights are positive and define a point x ∈ T.
  - θ(x) = θ, because θ ↦ [sin θᵢ / sin(θᵢ − Aᵢ)] is injective on H_A: each coordinate is decreasing, and Σθᵢ is fixed.
  - (This point was raised in audit/AUDIT_1.md §10.)
- **The incentre is the unique maximiser of D, with D = 1/4.** Each Wᵢ is critical at θᵢ° = (π + Aᵢ)/2. So θ° is a
  critical point of the strictly convex E on H_A, hence its unique minimiser, with E = 4.
- **Non-degenerate minima and four equilibria.** V_q has a non-degenerate planar minimum iff
  ψ_q ∈ F(M̂) + ℝ(1, 1, 1), by Step 1 and the bijection above.
  - In that case every other critical point has D < c (Step 2). So V_q is Morse with exactly four critical points.
  - F is injective on M̂ (Step 2).
  - Whether a degenerate point plus three saddles can occur is not needed for the theorem, and is not claimed.
- **Mixed signs (exploration: code/mixed_scan.py, code/mixed_boundary.py).**
  - For non-collinear charges of any signs, every equilibrium is planar. Off the plane one would need Σwᵢ = 0 and
    Σwᵢpᵢ = 0, i.e. w = 0.
  - With signed barycentric coordinates and signed angles σᵢ, the same algebra gives
    D = P_x Π sin σᵢ / (2Δ) and D·E = 1, and the signed separated form of the charges.
  - For q₃ < 0 < q₁, q₂: in the vertical angle at p₃ we have P_x < 0 and Π sin σ > 0. So D < 0 and there are no
    maxima.
  - Minima need D > c > 0, hence P_x < 0: they lie outside the circumcircle, beyond the side p₁p₂, in angle coordinates
    θᵢ ∈ (0, Aᵢ). Their centre is the excentre opposite p₃, where D = 1/4.
  - The sublevel set {E < (b + 1)²/b} is NOT convex in these coordinates for large α. Checked with mpmath: at α = 10 the
    bordered Hessian is −6·10⁻³. For α = ½ no violation was found. So this case needs another idea.
  - **Structure.** Put δ′ = (t₁, t₂, π − s) and γ′ = (A₁ − t₁, A₂ − t₂, π − A₁ − A₂ + s), with s = t₁ + t₂. Both are
    angle triples of triangles. Then E = e₂(cot δ′ + cot γ′) = 2 + 2 cosh ρ(δ′, γ′), with δ′ + γ′ = (A₁, A₂, π + A₃).
    This is the positive-case formula for the "triangle" A′ = (π − A₁, π − A₂, −A₃), which has a negative angle.
    Then W′₃ < 0, so the log-convexity argument fails. (For pairs with δᵢ + γᵢ > π, convexity of Ξ fails in general;
    see code/hyperbolic.py.)
  - **Numerics** (code/mixed_boundary.py). The boundary curvature of {E < k} is > 0 for α = ½, 1, 2 and 4 over 300
    random shapes (min angle 0.3°). It is < 0 for α = 10 at the shape (10.6°, 29.3°, 140.1°).
  - **Degenerate shapes.** curvature × radius ∝ (smallest angle). After rescaling that coordinate the limit looks
    positive.
  - **Plausible target.** Maxwell's bound for three charges of any signs, for α up to some α* ∈ (4, 10), including
    Coulomb. This needs a quantitative convexity proof, e.g. certified, with a blow-up at degenerate shapes.
  - **Outcome (MIXED.md).**
    - This plan works. The convexity is certified by interval arithmetic, with charts at the degenerate shapes.
    - That proves the Coulomb case (Theorem 2) and every α ∈ [1/4, 2] (Theorem 3).
    - Numerically, the convexity holds up to α ≈ 5.7.
- **Earlier route.** The previous formulation (det DΨ = β²(D − c)/(2Δλ₁λ₂λ₃), code/symbolic_I2.py; closedness of
  Σ log rᵢ dθᵢ) is consistent with, and superseded by, the separated form (3.1).

## 9. Numerical checks (independent of the proof)

- code/key_identities.py checks Lemma 3(b), the charge formula, Wᵢ = sin Aᵢ/(sin δᵢ sin θᵢ), D·E = 1,
  e₂(d) = (b + 1)² − bE, convexity of E, and the incentre facts.
  - Over 20,000 random points, relative errors are ≤ 1.2·10⁻⁹.
- code/separable.py checks the separated charge formula and sign(D − c) = sign(e₂(d)).
  - It also confirms that the negative-definite region is empty and that M̂ is star-shaped with boundary curvature
    > 0.
- code/scan.py and code/scan2.py: direct critical-point searches over more than 12,000 configurations (α = ½ to 3)
  never found two minima.
