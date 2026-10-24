# AUDIT 1: Maxwell's problem for three positive charges

- Scope:
  - **the rewritten PROOF.md**, file time 2026-10-24 10:01, "Status: the proof below is complete";
  - an appendix on the earlier version: the audit was started on it, and part of it is now superseded.
- Auditor: one independent checker.
- All scripts and outputs are in `audit/`. The scripts were written from scratch and import nothing from
  `code/`; the project's own scripts were only re-run, to reproduce their outputs.
- Every run used at most 2 cores. Precision is stated each time: "40 digits" means mpmath with mp.dps = 40.

## Verdict

**I found no gap in the rewritten proof.**
- Every lemma was re-derived by hand.
- Every identity was verified by independent code at 29–41 digits. Lemmas 5–6 were checked on the whole hexagon H_A,
  not only on θ(T).
- The chain "two planar minima (or a minimum and a degenerate point) ⇒ contradiction" was checked step by step. It was
  also attacked numerically: the only apparent violation was float64 rounding, and at 50 digits it is positive.
- Independent searches for a counterexample found none at any α tested (Coulomb and α ∈ {0.05, …, 10}).
- The remaining findings are NOTES about wording and remarks; none affects the theorem.

**Caveat.** This is an internal audit, not a referee report. The rewritten PROOF.md claims the sharp bound 4 for every
α > 0, which GNNS (30 July 2026) list as open. Before any announcement:
- **external expert review**, and
- a **priority check**: the literature after July 2026, including recent preprint releases and arXiv.

The key identities are elementary and would also be cheap to formalise.

| Part of PROOF.md | Verdict |
|---|---|
| §0 facts, collinear case, Lemma 1, index count | **OK** |
| Lemma 2 (θ injective) | **OK** |
| Lemma 3 (pedal identity, separated charge formula) | **OK** |
| Lemma 4 (D = 1/E, Leibniz step) | **OK** |
| Lemma 5 (strict convexity of E on H_A) | **OK** |
| Lemma 6 (e₂(d) = (b+1)² − bE; definiteness; sign at the incentre) | **OK** |
| §7 Step 1 (bookkeeping mod (1,1,1)) | **OK** |
| §7 Step 2 (at most one critical point with D ≥ c) | **OK** |
| §7 Step 3 (at most three saddles: Hess = 2ακ L∘DΨ, Sard, persistence) | **OK**; NOTE (wording) |
| §8 Remarks (not used in the proof) | NOTE: two remarks rely on unproved facts. Both facts are true numerically |
| Independent counterexample search | none found |

---

## 1. §0, the collinear case, Lemma 1 and the index count — OK

Scripts: `a1_lemma1.py` → `a1_lemma1.out`, and `a9_collinear_step2.py` (C) → `a9_collinear_step2.out`.

**Re-derivation.**
- ∇V = −2α Σ qᵢrᵢ^(−β)(x − pᵢ). A critical point is therefore a positive-weight barycentre. So it lies in the open
  triangle, in the plane; for collinear charges it lies on the open hull segment.
- Hess V_q = 2ακ(βS − I), and det(I − βS) = β²(D − c).
- The step "s > 1/β ⇔ s(1−s) > c" needs β ≥ 2, which holds.
- Since β(1−s) − 1 > 0 always, there is no planar maximum.
- The index count: the winding of ∇V is +1 on a large circle and +1 on each small circle around a charge, so the total
  index is −2 and s₀ = m₀ + 2.

**Numerical checks (independent code).**
- Gradient and Hessian against 40-digit finite differences of V in ℝ³: relative errors 8e-16 and 1e-15.
- q = Q(x) makes x critical, to 9e-16.
- 246 critical points found by my own root finder all satisfy Q(x) ∝ q (to 2e-12), and each type equals sign(D − c).
- 2000-point identities (tr S = 1, Cauchy–Binet, D = (1−|Z|²)/4, det(I−βS) = β²(D−c)) hold to ≤ 8e-15.
- 3D root finding: all 60 critical points lie in the plane inside T, and the normal Hessian eigenvalue is −2ακ.
- Planar Laplacian = 4α²Σqᵢrᵢ^(−2α−2), to 9e-16.
- Index bookkeeping, using winding numbers on circles:
  - 120 configurations, 60 of them with a minimum forced (q = Q(x₀), x₀ ∈ M);
  - tallies (0,2)×60 and (1,3)×60, with 0 mismatches.
- Collinear charges in ℝ³ (282 configurations, α from 0.05 to 10): V′ on the line has exactly 2 zeros every time, and
  3D root finding finds 0 critical points off the line.

## 2. Lemma 2 (injectivity of x ↦ θ) — OK

- The circle argument is correct.
  - The level sets lie on two circles through p₃.
  - The circles are distinct, because a common circle would be the circumcircle, which misses the open triangle.
  - Two distinct circles meet in at most two points, one of them p₃.
- Numerically (`a4_consequences.py`), the closed-form inverse reproduces x from θ(x) at 3000 points. The inverse is
  λᵢ ∝ sin Aᵢ sin θᵢ / sin δᵢ, which is exactly what Lemma 3(b) gives.
- Random θ in H_A map back into T with the same angles (0/3000 failures). So **θ(T) = H_A**. The proof only uses
  θ(T) ⊂ H_A.

## 3. Lemma 3 (pedal identity) — OK

Script `a8_new_identities.py` → `a8_new_identities.out`: 400 points of T, half of the triangles with an angle down to
0.5°, a quarter of the points near an edge and a quarter near a vertex, α from 0.05 to 10; 40 digits.

**Proof steps re-derived.**
- The inversion maps rays to themselves, so x is inside T*.
- By similarity (ii), ∠x pᵢ* pⱼ* = ∠x pⱼ pᵢ. Then the angle of T* at pᵢ* is (Aⱼ − ∠xpⱼpₖ) + (Aₖ − ∠xpₖpⱼ) = θᵢ − Aᵢ.
- The circumcircle (not through x) maps to a circle of radius R/P_x.
- The law of sines plus (i) gives aᵢ/(rⱼrₖ) = 2(R/P_x) sin δᵢ.
- Hence (b), and with λᵢrᵢ = r₁r₂r₃ sin θᵢ/(2Δ), (c).

**Numerics.**

| Check | Result |
|---|---|
| δᵢ > 0 | all points |
| \|Σδ − π\| | ≤ 3e-35 |
| (b), relative | ≤ 1.2e-32 |
| (c), relative | ≤ 2.3e-29 |
| Angles of T* = δᵢ | ≤ 3e-34 |
| Circumradius of T* = R/P_x | ≤ 2e-34 |

## 4. Lemma 4 (D = 1/E) — OK

- **Steps 1–6 re-derived.**
  - Wᵢ = rⱼrₖ/(P_x sin θᵢ).
  - E = r₁r₂r₃ Σrₖ sin θₖ / (P_x² Π sin θ).
  - D = r₁r₂r₃ Π sin θ Σrₖ sin θₖ / (2Δ)².
  - r₁r₂r₃ Σ rₖ sin θₖ = 2Δ Σ λₖrₖ².
  - **Leibniz step:** Σλₖ|x − pₖ|² = |x−O|² − 2⟨x−O, Σλₖ(pₖ−O)⟩ + R² = R² − |x−O|² = P_x, because Σλₖ(pₖ−O) = x − O.
    It is correct.
- **Numerics.**
  - Wᵢ = cot δᵢ − cot θᵢ, to 9e-40.
  - Σλₖrₖ² = P_x, to 9e-35 relative.
  - D (computed directly as det S) · E − 1 ≤ 2.3e-29.
  - sign(D − c) ≠ sign(e₂(d)) at 0 points.

## 5. Lemma 5 (strict convexity of E on H_A) — OK

- ℓᵢ″ = csc²(t − Aᵢ) + csc² t > 0, and WᵢWⱼ = exp(ℓᵢ + ℓⱼ) with Hess e^f = e^f(Hess f + ∇f∇fᵀ) ⪰ e^f Hess f. This
  is correct.
- In fact every nonzero vᵢ already contributes WᵢWⱼℓᵢ″vᵢ² > 0. So the "two nonzero entries" remark is not even needed.
- The sublevel-set consequences follow: the open segment between two points of {E ≤ k} lies in {E < k}.
- **Numerics on 400 random points of H_A** (not necessarily images of points of T) and random directions with Σv = 0:
  - ℓᵢ″ matches the formula to 5e-41;
  - E″ > 0 at all points, and E″ ≥ the proof's lower bound at all points;
  - the minimal ratio of E″ to the bound is ≈ 1.
- This analytic convexity explains all of my earlier numerics: θ(M) = M̂ is convex for every shape scanned.

## 6. Lemma 6 (e₂(d) and definiteness) — OK

- **Re-derived.**
  - e₂(cot) = 1 for three angles summing to π or to 2π. No division by zero occurs, since θ₁+θ₂ ∉ πℤ.
  - e₂(X − C) = 2 − B and e₂(C + bX) = 1 + b² + bB.
  - So e₂(d) = (b+1)² − bE.
  - The Gram matrix [[d₁+d₃, d₃], [d₃, d₂+d₃]] has determinant e₂(d).
  - On M̂ the form is definite, and M̂ is convex hence connected.
  - At θ° = ((π+Aᵢ)/2): Wᵢ = 2tan(Aᵢ/2), E = 4 < (b+1)²/b ⇔ (b−1)² > 0 (true since b = 2α+1 > 1), and
    dᵢ = (b−1)tan(Aᵢ/2) > 0.
- **Numerics.**
  - On θ(T) and on random points of H_A: e₂(cot δ) = e₂(cot θ) = 1 to ≤ 1e-28.
  - e₂(d) = (b+1)² − bE to 4e-35 (θ(T)) and 2e-39 (H_A).
  - Gram determinant = e₂(d) to 4e-38.
  - Incentre (100 triangles): θᵢ = (π+Aᵢ)/2, Wᵢ = 2tan(Aᵢ/2), E = 4 and dᵢ = (b−1)tan(Aᵢ/2), all to ≤ 1e-37.
- This also settles the sign question I raised against the earlier version (appendix, A4): with ccw labelling,
  det Hess_θ H has the sign of D − c, and the form is *positive* definite on M̂.

## 7. §7 Step 1 (bookkeeping mod (1,1,1)) — OK

- By (3.1), F(θ(x)) − ψ_q = (log Qᵢ(x) − log qᵢ)ᵢ − log C(x)(1,1,1). At a critical point this lies in ℝ(1,1,1). For
  two critical points of the same V_q, F(θ) − F(θ′) ∈ ℝ(1,1,1), and pairing with v = θ − θ′ (Σvᵢ = 0) kills it.
- Numerically, log Qᵢ − fᵢ(θᵢ) + b log sin Aᵢ − log C(x) = 0 to 2.3e-29, for each i separately.

## 8. §7 Step 2 (at most one critical point with D ≥ c) — OK

**The chain.**
1. x ≠ x′ critical with D ≥ c.
2. Lemma 2 gives θ ≠ θ′.
3. Lemma 4 gives E(θ), E(θ′) ≤ (b+1)²/b, with θ, θ′ ∈ H_A by Lemma 3(a).
4. Lemma 5 puts the open segment in M̂.
5. Lemma 6 makes Σdᵢvᵢ² > 0 there.
6. The fundamental theorem of calculus gives 0 = Σ(fᵢ(θᵢ) − fᵢ(θ′ᵢ))vᵢ = ∫₀¹ Σdᵢ(θ′+tv)vᵢ² dt > 0, a contradiction.

The integrand is continuous on [0, 1], because the closed segment lies in H_A, where the fᵢ are smooth.

**No hidden assumptions.**
- The segment need not lie in θ(T): F, E and d are defined by explicit formulas on all of H_A.
- A degenerate point (D = c, i.e. E = (b+1)²/b) is covered, because the argument uses the *closed* sublevel set and
  the *open* segment.
- Two minima, and a minimum plus a degenerate point, are excluded alike.

**Adversarial numerics** (`a9_collinear_step2.py` (M), then `a9b_step2_mp.py`).
- I sampled 74,850 pairs in {E ≤ k} over 150 shapes and α from 0.05 to 10: interior pairs, boundary pairs, and nearly
  coincident boundary pairs. A Nelder–Mead search then minimised m = Σ(fᵢ(θᵢ) − fᵢ(θ′ᵢ))vᵢ/|v|² over boundary pairs.
- In float64 the minimum came out slightly **negative** (−3.3e-6 of the median). I re-evaluated the minimising pair of
  each of 36 shapes at 50 digits, recomputing the boundary points by 50-digit bisection.
- **All 36 are positive.** The minimum is 5.5e-15 of the median, at separations 2.8e-9 to 8.4e-6. Each value equals the
  tangential limit vᵀdiag(d)v/|v|² of the form.
- So the optimiser collapses pairs onto the points of ∂M̂ where the form's null direction is tangent to ∂M̂ (the cusps
  of the fold image). There inf m = 0 is approached but not attained, which is exactly what Step 2 predicts.

## 9. §7 Step 3 (at most three saddles) — OK, NOTE

- **Hess V_q = 2ακ L∘DΨ.**
  - Write qᵢrᵢ^(−β) = κλᵢe^{εᵢ}, so that ∇V_q = −2ακΣλᵢ(e^{εᵢ} − 1)(x − pᵢ).
  - Differentiating at ε = 0, the common part of ∇ε drops out because L(1,1,1) = Σλᵢ(x − pᵢ) = 0. This gives
    2ακ Σλᵢ(x−pᵢ)⊗∇log Qᵢ.
  - L has kernel exactly ℝ(1,1,1), because the three vectors x − pᵢ are pairwise independent and their only relation is
    λ.
  - Numerically, at 120 critical points (α ∈ {0.25, 0.5, 1, 2.5}): |Hess − 2ακ L DΨ|/|Hess| ≤ 1.2e-29, L(1,1,1) = 0 to
    8e-31, and rank L = 2.
- **Sard.** V_q′ has a degenerate critical point ⇔ [log q′] is a critical value of Ψ: T → ℝ³/ℝ(1,1,1) ≅ ℝ². The
  critical values have measure zero, so Morse q′ are dense. Correct.
- **Persistence.** Non-degenerate saddles persist and stay distinct under small perturbations (implicit function
  theorem), and D < c is an open condition. The index count and Step 2 for q′ give s₀(q′) = m₀(q′) + 2 ≤ 3.
- **NOTE (wording).** "Take a Morse q′ close enough …; each saddle of V_q persists" implicitly treats the saddles of
  V_q as finitely many. Better: "if V_q had four saddles, choose q′ close enough that these four persist". This needs
  no finiteness, and the critical set then turns out finite anyway.
- **Conclusion.** At most 1 + 3 = 4 critical points. This counts *all* critical points, degenerate ones included,
  which is stronger than the non-degenerate count in STATEMENT.md. Correct.

## 10. §8 Remarks — NOTE (not used in the proof)

- "**Exactly** those q with ψ_q ∈ F(M̂) + ℝ(1,1,1) give four equilibria" needs two things.
  - (a) M̂ ⊂ θ(T), i.e. θ(T) = H_A. This is true numerically (§2 above). A proof could start from the explicit inverse
    λᵢ ∝ sin Aᵢ sin θᵢ / sin δᵢ.
  - (b) An argument that a non-Morse q cannot have four critical points (a degenerate point plus three saddles).
    - I believe (b) holds. The degenerate point is the unique point with D ≥ c, and its local degree is ≤ 0: degree
      +1 together with the negative sheet that exists across ∂M (∇D ≠ 0 there, because E is strictly convex) would
      force two minima for a nearby q, contradicting Step 2.
    - So the index count gives at most 2 + 0 saddles. This is a sketch; the remark should either cite such an argument
      or drop "exactly".
- "The incentre is the unique maximiser of D." This needs E ≥ 4 on H_A, with equality only at θ°. That follows from
  the hyperbolic identity E = 2(1 + cosh ρ) in the remark, or from D ≤ 1/4 together with θ(T) = H_A.
  - Numerically it holds. `a8` gives E = 4 (so D = 1/4 and Z = 0) at the incentre of 100 random triangles. The shape
    scans find exactly one zero of Z per shape. For thin isosceles shapes that zero has barycentrics
    (¼, ¼, ½) = (a₁ : a₂ : a₃), the incentre.
- The other remarks (hyperbolic picture, |Z| = tanh(ρ/2), {ρ < log b}) are consistent with D = (1 − |Z|²)/4 and with
  E < (b+1)²/b ⇔ cosh ρ < cosh log b. They are not used.

## 11. Independent counterexample search — none found

| Search | Configurations | Result |
|---|---|---|
| (a) Direct, Coulomb (`a6_direct_search.py`, no reduction used) | 4000, q = Q(x₀) with x₀ ∈ M (half near the fold), angles ≥ 3° | all (1 min, 3 saddles); 0 index-count failures; 0 with ≥ 2 minima |
| (b) Direct, α = 0.1 / 2.5 / 10 | 800 each | α = 0.1: all (1,3) |
| | | α = 2.5 and 10: 2 and 12 first-pass finder misses. Re-solving exactly those (`a6d_failures.py`: denser seeds, log-barycentric solver) gives (1,3) in all 14 |
| (c) All charges at once, per shape (`a67_shape_scan.py`) | Coulomb: 719 shapes on a 2° grid (A₁ ≤ A₂ ≤ A₃, min 0.5°) plus thin extras down to 0.2° | every shape: one zero of Z (the incentre), one maximiser of D above c, θ(M) star-shaped, D monotone on rays, boundary turning one-signed, max winding number 1 |
| | α ∈ {0.05, 0.25, 1, 2.5, 10}: 99 shapes each | same result |

Notes on these searches:
- In (b), 10 of the 12 saddles missed at α = 10 lie within about 1e-10 (in barycentric coordinates) of an edge. The
  other two, and the two misses at α = 2.5, are close to an edge or a vertex.
- How (c) works: Ψ preserves orientation on M, so for every q off the fold image the number of minima equals the
  winding number of Γ = Ψ(∂M). It is evaluated at 10,800 interior points per shape.
- A few thin shapes show float64 "self-intersections" of Γ, with segments of length ~1e-10 near the symmetric cusp.
  `a6c_thin_mp.py` (40 digits) shows they are rounding: Ψ₁ − Ψ₂ has the sign of θ₁ − θ₂ at every boundary point, with
  values down to 1e-21.

## 12. Precision warning for anyone recomputing near degenerate shapes

- In θ-coordinates, θ₃ − A₃ is computed by cancellation for thin triangles.
- My own float64 θ-space code gave **false non-convexity** for thin isosceles shapes at A ≤ 0.1° (and at 0.2° for
  α = 0.05). At 50 digits (`a7b_thin_mp.py`, `a7c_thin_general_mp.py`) the same shapes are convex.
  - The relative margin is ≈ 0.46 for Coulomb, independent of the shape as A → 0 (an elliptic limit profile).
  - The margin is 0.50 / 0.49 / 0.30 / 0.14 at α = 0.05 / 0.25 / 2.5 / 10.
- Fold-image cusps are extremely sharp (Ψ₁ − Ψ₂ ~ 1e-15 to 1e-21 at A = 0.2°).
- This no longer matters for the proof, since Lemma 5 is analytic. It does matter for any numerical cross-check.

---

## Appendix: findings on the earlier PROOF.md (superseded)

The earlier version had: Lemma 2 (det identity), Lemma 3 (closed 1-form), an unproved convexity Lemma 4, and a
deduction.

- **A1. Lemma 1:** OK, as in §1.
- **A2. det DΨ = β²(D − c)/(2Δλ₁λ₂λ₃): OK.**
  - `code/symbolic_I2.py` re-run gives 0. Its algebra was checked line by line.
  - My own sympy re-derivation (6 free vertex symbols, Cramer barycentrics, D = det S) expands to 0. 50-digit
    numerics agree to 4e-48.
  - NOTE: Δ must be the *signed* area (ccw labelling). For cw labelling the sign flips.
- **A3. Σ log rᵢ dθᵢ closed: OK.**
  - θᵢ = φₖ − φⱼ (ccw), the Re(f′ḡ′) formula and the cancellation are all correct. The curl is 0 to 30 digits, and
    loop integrals are ≤ 1.6e-15.
  - Σ log sin θᵢ dθᵢ = d(−ΣЛ(θᵢ)).
  - DΨ·A⁻¹ is symmetric to 6e-30.
  - The explicit Lobachevsky form holds with linear part −Σ log(sᵢ)θᵢ (residual 7e-16).
  - Lemma 3 of the new version makes this transparent: Ψ ≡ F(θ) mod (1,1,1) with F diagonal.
- **A4. ISSUE (minor), now resolved.**
  - "det Hess_θ H has the sign ±(D − c) … therefore definite on M" required the sign to be +, and this was not shown.
    It is + (det A = 2·area of the inverted triangle > 0 for ccw; checked at 20,000 points).
  - The new Lemma 6 (incentre evaluation) settles it.
- **A4. NOTE, now resolved.** Injectivity of x ↦ θ was not argued; it is now Lemma 2.
- **A5. Deduction given old Lemma 4: OK**, with NOTES that are now all handled by the new §7:
  - the density of Morse q (Sard);
  - degenerate minima (Step 2 uses D ≥ c);
  - the observation that this route needs convexity in θ specifically.
- **A7. Old Lemma 4 numerics.** Reproduced: `run_extreme_shapes_a0.5.out` and `run_convexity_test_a0.5.out` give
  convex shapes with turn ratio 0.125–0.162.
  - Their ray method presupposed connectivity and star-shapedness of M, which those scripts did not check; my scans
    confirm both.
  - The "turn ratio" is not a convexity margin.
  - The question is now moot: θ(M) = M̂ is a sublevel set of the strictly convex E (Lemma 5).
  - The project's earlier `scan_bal_a20.out` (1466/3000 parity failures) and `scan_bal_a0.5.out` (8 failures) show
    that its direct finder breaks down at large α and for near-degenerate shapes. Those scans are not used by the
    proof.

## Files (all in `audit/`)

| Script | Output | Checks |
|---|---|---|
| `a1_lemma1.py` | `a1_lemma1.out` | Lemma 1, standard facts, index count |
| `a2_lemma2.py` | `a2_lemma2.out` | old det identity; plus `run_symbolic_I2.out` |
| `a3_lemma3.py` | `a3_lemma3.out` | old closed-form lemma |
| `a4_consequences.py` | `a4_consequences.out` | θ injective, θ(T) = H_A, signs |
| `a8_new_identities.py` | `a8_new_identities.out` | new Lemmas 3–6, Step 1, Hess = 2ακ L∘DΨ |
| `a9_collinear_step2.py` | `a9_collinear_step2.out` | collinear case; Step-2 stress test |
| `a9b_step2_mp.py` | `a9b_step2_mp.out` | 50-digit resolution of the float64 negatives |
| `a6_direct_search.py` | `a6_direct_search_a{0.5,0.1,2.5,10}.out` | direct counterexample search |
| `a6d_failures.py` | `a6d_failures_a{2.5,10}.out` | robust re-solve of the large-α finder misses |
| `a67_shape_scan.py` | `a67_scan_a0.5_step2.out`, `a67_scan_a{0.05,0.25,1,2.5,10}_step6.out` | per-shape winding-number scans |
| `a67_thin_trend.py` | `a67_thin_trend_a0.5.out` | float64 thin-shape trend |
| `a6b_thin_diag.py` | `a6b_thin_diag.out` | fold-image self-intersection diagnostics |
| `a6c_thin_mp.py` | `a6c_thin_mp.out` | 40-digit fold-image check |
| `a7b_thin_mp.py`, `a7c_thin_general_mp.py` | `a7b_thin_mp_A*.out`, `a7c_thin_general_mp.out` | 50-digit thin-shape convexity |
| project scripts (re-run only) | `run_extreme_shapes_a0.5.out`, `run_convexity_test_a0.5.out` | reproduction of the project's outputs |
