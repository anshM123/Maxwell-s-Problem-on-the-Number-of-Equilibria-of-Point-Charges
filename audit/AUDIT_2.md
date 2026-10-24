# AUDIT 2: Three point charges of mixed signs, Coulomb (MIXED.md, Theorem 2)

**Scope.**
- `MIXED.md` (file time 2026-10-24 12:21).
- The certificate: `code/ia.py`, `code/mixed_cert.py`, `code/mixed_charts.py`, `code/mixed_cert_all.py`.
- This continues AUDIT_1 (positive charges, PROOF.md).

**Rules followed.**
- All verification code is my own, in `audit/b*.py`; outputs are in `b*.out`.
- The project's modules are imported only where they are *the object under test*: the certificate re-run, the `ia.py`
  unit tests and the chart-containment tests. No existing file was modified.
- Every run used at most 2 cores.
- The certificate was re-run from a scratch working directory, because running it from `code/` would overwrite the
  existing `cert/claimM_e0_0.2_r0_0.05.npz`.

## Verdict

**No gap was found in the logic of Theorem 2.**
- Every identity in §§1–6 was re-derived by hand and verified by independent 40–50-digit code.
- Lemma C was re-derived symbolically for arbitrary functions.
- The certificate is sound as implemented, and its runs reproduce MIXED.md's table exactly. More importantly, I
  **re-certified all of its boxes independently in rigorous interval arithmetic** (`mpmath.iv`, 80 bits), using my
  own implementation of the formulas. That removes the "libm within 2 ulps" assumption.
- A direct counterexample search (no reduction used, 50-digit verification) and a winding-number scan over 941 shapes
  found nothing.

**Findings.**
- **One ISSUE (minor):** the driver leaves a 1.2e-16 coverage sliver below A₁₂ = π. The one-line fix has been verified
  to certify, in their code and in mine.
- **Three NOTES:** wording only.

**Caveat.** This is an internal audit, not a referee report. Claim M is computer-assisted and holds for Coulomb only;
MIXED.md itself shows that M′ is non-convex at α = 10. Together with Theorem 1 this would settle Maxwell's problem for
three charges in ℝ³, which GNNS (July 2026) list as open. **External expert review and a literature/priority check**
(including recent preprint collections and arXiv) are needed before any announcement.

| Item | Verdict |
|---|---|
| 1. §§1–5: planarity, signed D / pedal / separated charges, R₊/R₋, no maxima, minima in R₊ᵒᵘᵗ, injectivity, E = 1/D, e₂(d), sign at the excentre, index count incl. Q_tot = 0 | **OK** (NOTE N1) |
| 2. §6.1: product representation, (s, Δ) diffeomorphism, M′ = {cosh Δ < β(s)}, ψ_max < A₁₂, ζ quadratics (vertex 3 and 1) | **OK** (NOTE N2) |
| 3. §6.2 Lemma C: bordered Hessian J, reductions to I and J₁/U₂², coordinate change X = (s, −t₂), signs of Ũ, mirror | **OK** |
| 4. §6.3 certificate: ia.py soundness, chart formulas (incl. band scaling and mean-value enclosure), coverage, re-run, high-precision spot checks | **OK**, with **ISSUE (minor)**: the π-sliver in coverage; fix verified |
| 5. §7 assembly: degenerate points, saddles (Hess = 2ακ L∘DΨ, Sard, persistence), collinear perturbation | **OK** (NOTE N3) |
| 6. Independent counterexample search (Coulomb, mixed signs, incl. near-degenerate triangles) | **none found** |

---

## 1. §§1–5 — OK

Script `b1_sections1to5.py` → `b1_sections1to5.out` (40 digits; my own geometry: signed angles via arg, barycentrics
via cross products, D = det S from the vectors). The inverse t → x is computed *geometrically*: inversion at p₃ turns
both angle loci into lines. So Lemma P is not assumed anywhere.

**Re-derived by hand.**
- §1 planarity: h ≠ 0 forces Σwᵢ = 0 and Σwᵢpᵢ = 0, hence w = 0. This works in every ℝⁿ.
- §2: w = κλ with κ = Σw ≠ 0, which gives the regions R± from the sign pattern.
- Hess V = 2ακ(βS − I), and extremum ⇔ (βs₁−1)(βs₂−1) > 0 ⇔ D > c. Here S is not PSD, but the determinant identity
  alone suffices.
- Lemma D (with the Leibniz step for signed weights).
- No maxima: in R₋, p₃ is a strict convex combination of x, p₁, p₂, so P_x < 0.
- Lemma P and the separated charges; ∇H = ψ_q.
- W = (−U₁, −U₂, U₃), E = e₂(W) = 1/D, d = x̃ − bỹ, and e₂(d) = (b+1)² − bE.
- The trace at t° is (1−b)(cot(A₁/2) + cot(A₂/2) − 2cot(A₁₂/2)) < 0, so the form there is negative definite.

**Numerics.**

| Check | Result |
|---|---|
| Planarity | 77 critical points from 3D root finding (off-plane seeds), 0 off the plane |
| Critical-point locations | 120 critical points: 73 in R₊, 47 in R₋, 0 elsewhere |
| w = κλ | wᵢ/λᵢ constant to 2e-13; κ sign matches the region every time |
| Hess V = 2ακ(βS − I) | to 4e-14 |
| Lemma D | to 1e-35 |
| Lemma P (all regions) | to 3e-37 |
| Separated charges | to 6e-37 |
| ∇H(t(x)) = ψ_q | to 6e-37 |
| Signs | R₋: Π sin σ > 0, P_x < 0, D < 0 at every sample; R₊: Π sin σ < 0; R₊ᵒᵘᵗ: tᵢ ∈ (0, Aᵢ) and s = t₁ + t₂ < A₁₂ (0 failures) |
| Injectivity and image | the geometric inverse reproduces x to 1e-37; 1,500 random t ∈ Q all map into R₊ᵒᵘᵗ with t(x(t)) = t to 1e-38, so t is a bijection R₊ᵒᵘᵗ → Q |
| W, E | W = (−U₁, −U₂, U₃) to 2e-39; E·D − 1 to 5e-35 |
| Hyperbolic form | e₂(x̃) = e₂(ỹ) = 1 to 2e-35; E = 2 + 2⟨x̃, ỹ⟩ to 5e-38; E ≥ 4 on Q (0 violations); t(excentre₃) = (A₁/2, A₂/2) to 3e-39; E(t°) = 4 |
| d and Hess H | d = (φ₁′, φ₂′, −φ₃′) to 3e-37; e₂(d) = (b+1)² − bE to 1e-37; Hess H (from ∇H) = Gram(d) to 1e-40; negative definite at t° for all 200 shapes |

**Index count.** Over 80 configurations, the winding is 1 on each small circle for either sign. On a large circle it is
1 if Q_tot ≠ 0, giving (0 extrema, 2 saddles) in 39 cases, and 2 if Q_tot = 0 (the dipole), giving (0, 1) in 40 cases.
- One Q_tot ≠ 0 case showed winding 2 on my r = 400 circle. `b1b_far_saddle.py` shows the cause: for small Q_tot the
  winding changes from 2 to 1 across a saddle at radius ≈ 2|m|/|Q_tot| (27, 274 and 2,746 for Q_tot = 1e-1, 1e-2,
  1e-3). This is consistent with s₀ = m₀ + 2.

**NOTE N1 (exposition).** §4 justifies "x̃ and ỹ lie in the same sheet" by "E = 1/D > 0". That only covers t ∈ t(R₊ᵒᵘᵗ).
This is all of Q (verified above), but a direct argument is cleaner: x̃ and ỹ are the cot-vectors of the genuine
triangles with angles (t₁, t₂, π − s) and (A₁ − t₁, A₂ − t₂, A₃ + s). The proof does not need E ≥ 4 anyway: it uses
E(t°) = 4 < K, and §6.1's product form gives E > 0 on Q.

## 2. §6.1 — OK

Script `b2_section61.py` → `b2_section61.out` (40 digits).

| Check | Result |
|---|---|
| Λᵢ′ = Uᵢ; Uᵢ = 2(cosh Λᵢ + cos aᵢ)/sin aᵢ; Uᵢ′ = 2 sinh Λᵢ Uᵢ/sin aᵢ; t = arg(1 + e^{Λ+ia}) | to ≤ 5e-40 |
| E = (2U₃/sin A₃)(cosh Δ + cos A₃) | to 2e-38 |
| K − E = (2U₃/sin A₃)(β(s) − cosh Δ), i.e. E < K ⇔ cosh Δ < β(s) | to 2e-38 |
| Vertex-1 and vertex-2 representations | to 3.5e-37 |
| (s, Δ) Jacobian det = −(U₁ + U₂) | to 3e-41 |
| Δ along s-lines | strictly increasing on 200 lines × 400 points; tends to ∓∞ at the ends |
| ψ_max < A₁₂ | 0 failures in 3,000 shapes, angles down to 0.01° |
| β ≥ 1 ⇔ \|ψ\| ≤ ψ_max | 0 failures |
| β − 1 = (F+1) sin(·) sin(·) | to 4e-41 |
| Bounds on M′ (3,600 sampled points, 300 shapes) | P ≥ 3.502 ≥ 3.5; max \|Λ₁\|, \|Λ₂\|, \|Λ₃\| = 0.688 < log 2 = 0.693; max \|Δ\| = 0.674; no sampled point of {\|Δ\| < B} has E ≥ K |
| Vertex-3 quadratic | exactly one positive root (400 cases); t₁ + t₂ = s and E = K to 4e-38; agrees with an independent bisection to 2e-40 |
| Vertex-1 quadratic | exactly one positive root; ψ₁,max < A₁; E = K to 3e-37; Λ₃ − Λ₂ = D to 1e-40 |

**NOTE N2.** "t₂ = arctan(ζ sin A₂/(1 + ζ cos A₂))" takes the wrong branch (t₂ − π) when 1 + ζ cos A₂ < 0. This
happened for 5 of 400 random shapes, all with A₂ > π/2. It should read t₂ = arg(1 + ζe^{iA₂}). It is harmless for the
certificate:
- the corner chart has A₂ ≤ r₀ = 0.05, so cos A₂ > 0;
- E, the Λ's, U₃ and U₃′ are all invariant under t₂ → t₂ − π.

## 3. §6.2 Lemma C — OK

Script `b3_lemmaC.py` → `b3_lemmaC.out`.

- **Symbolic (sympy), with *arbitrary* smooth β, Λ̃₁, Λ̃₂:**
  - τ·∇Ĝ = 0;
  - −Hess Ĝ(τ, τ) − J ≡ 0;
  - J/(U₁U₂)² − I ≡ 0, with cᵢ = 1/Uᵢ and aᵢ = Uᵢ′/Uᵢ².
- **Vertex-1 reduction, re-derived by hand.** With X = (s, −t₂): X₁ + X₂ = t₁, Λ̃₁ = Λ₃, Λ̃₂(X₂) = Λ₂(−X₂). Then
  Ũ₁ = U₃, Ũ₂ = −U₂, Ũ₂′ = +U₂′(t₂), and U₃′ = −U₃(cot s + cot(s + A₃)). Dividing by U₂² gives the stated J₁/U₂².
- **Why the sign of J is the intrinsic curvature sign.**
  - Ĝ₃ = (K − E) sin A₃/(2U₃) and Ĝ₁ = (K − E) sin A₁/(2U₁). Both have the zero set ∂M′ with ∇ ≠ 0 there, and they
    differ by positive factors.
  - On the zero set, Hess(gG)(τ, τ) = g Hess G(τ, τ).
  - The coordinate change is linear, which preserves local convexity.
- **Quantitative check (not only signs).** At 400 boundary points for each α, found by my own bisection on E itself:

  | Check | α = ½ (K = 4.5) | α = 10 (K = 23.05) |
  |---|---|---|
  | κ = I·(U₁U₂)²/\|∇Ĝ\|³, κ from E's derivatives (vertex 3) | to 1e-39 | to 3e-38 |
  | J₁/\|∇Ĝ₁\|³ = curvature in X-coordinates (vertex 1) | to 5e-38 | to 2e-36 |
  | Sign disagreements | 0 | 0 |
  | Negative points | none | one point where I, J₁ and κ are all negative |

  The α = 10 negative point is consistent with MIXED.md §8.
- **Mirror.** I(shape, t) = I(mirrored shape, swapped t) exactly (difference 0.0). So certifying the upper arc for all
  (A₁₂, λ) covers the lower arc.

## 4. §6.3 certificate — OK, one ISSUE (minor, fix verified)

**(a) `ia.py` soundness.**
- I read it line by line:
  - outward rounding: 2 ulps for IEEE operations, 4 ulps for library transcendentals;
  - sin/cos include interior extrema, with a 1e-12 guard band for the float π/2 and 2π;
  - sinc, asinc and asinhc use their monotonicity in |t| (true facts), with a 5-ulp widening of the float ratio;
  - div by an interval containing 0 returns (−∞, ∞), and recip raises;
  - sqrt, arccos and arcsin clip, which is sound because the true arguments lie in the domain;
  - NaN never certifies (`ok &= isfinite`).
- Unit tests against 50-digit true ranges (`b4b_ia_tests.py`): 1,500 random intervals per operation (narrow, medium
  and wide; containing 0 or extrema), with **0 failures**. The edge cases are correct.
- Measured libm error on this machine (numpy 2.5.3): at most **1.26 ulp** over 11 functions × 20,000 arguments, against
  an assumed 2 ulp and a 4-ulp widening.

**(b) Chart formulas.**
- All three were re-derived by hand:
  - main: the quadratic in rationalised form and β − 1 as a product of sines;
  - band: ψ_max = 2 arcsin(q sin(e/2)), the e-scalings of S, β′, cᵢ and Δ, the I/e⁴ expansion, and the mean-value
    split of (a₁ − a₂)/e with ∂a/∂L = (1 + cos A cosh L)/(cosh L + cos A)²;
  - corner: the vertex-1 quadratic.
- Containment test (`b4c_chart_containment.py`): their chart functions on random boxes (widths 1e-7 to 1e-1 of the
  domain, plus point boxes) against my independent references (boundary point by bisection, not by their quadratic).
  **0 failures in 4,000 reference points** across main, band and corner ±.

**(c) Coverage.**
- The three charts plus mirrors cover all shapes: band A₁₂ ≤ e₀; main A₁₂ ≥ e₀; corner A₁ ≥ π − r₀, and its mirror for
  A₂ ≥ π − r₀.
- A main-chart box is dropped only if it lies wholly in a corner region; the 1e-15 margin absorbs the float rounding of
  the product.
- For shapes with A₁, A₂ < π − r₀, no box containing the shape or its mirror is ever dropped, so both arcs are
  certified.
- The chart boundaries overlap: e₀ appears in both band and main, and π − r₀ in both main and corner.

**ISSUE (minor): the π-sliver.**
- `mixed_cert_all.branch_and_bound` clamps the main chart's top A₁₂ edge to `math.pi`, which is 1.22e-16 *below* π
  (`b4a_where_small.out`).
- So triangles with 0 < A₃ < 1.22e-16 and A₁, A₂ < π − r₀ are covered by no chart.
- **Fix:** use `ia.up(np.pi)`, as `mixed_cert.run` already does.
- Verified both ways:
  - with the fix, their code still certifies (127,403 boxes, 0 failed; `b4d_coverage_sliver.out`);
  - my rigorous evaluator certifies all 728 top-layer boxes extended to the float above π (`b4g_sliver_mpiv.out`).

**(d) Re-runs** (`b4_run_cert_0.2_0.05.out`, `b4_run_cert_other_params.out`). All three rows of MIXED.md's table are
reproduced exactly.

| e₀ | r₀ | Main boxes | Band | Corner ± | Failed | Time |
|---|---|---|---|---|---|---|
| 0.2 | 0.05 | 127,403 | 2,048 | 1,024 + 1,024 | 0 | 0.8 s |
| 0.1 | 0.1 | 2,027,388 | 2,048 | 1,024 + 1,024 | 0 | 19 s |
| 0.4 | 0.02 | 21,711 | 2,048 | 1,024 + 1,024 | 0 | 0.3 s |

The smallest certified lower bound (2e-10) sits at A₁₂ ≈ 0.20–0.32, where I = O(e⁴) and the bisection stops at the
first positive bound. It is not a near-failure.

**(e) Independent rigorous re-certification.** This replaces "spot-check boxes with mpmath".
- I wrote the three chart quantities from scratch in `mpmath.iv` (80-bit, correctly rounded).
- sinh, cosh, acos, sinc, asinc and asinhc are built from the rigorous primitives by monotonicity, and these helpers
  were tested with 0 failures (`b4e_helper_test.out`).
- Every box of the project's run was evaluated, with local bisection allowed.

| Chart | Boxes | Result | Smallest lower bound |
|---|---|---|---|
| Main (`b4e_mpiv_recert_main.out`) | 127,403 | all certified; 104 needed a local bisection; 0 failures | 2.04e-10 |
| Band (`b4f_mpiv_band_corner.out`) | 2,048 | all certified at first evaluation | 0.00934 |
| Corner + | 1,024 | all certified at first evaluation | 1.574 |
| Corner − | 1,024 | all certified at first evaluation | 1.548 |

- The smallest lower bounds agree with the project's to all printed digits, which independently confirms that both
  implementations compute the same quantities.
- The main chart's boxes tile [e₀, π] × [0, 1] × [−1, 1] minus the dropped corner boxes (volume deficit 5.28e-4).
- Claim M (Coulomb) therefore stands without the libm assumption, once the π-sliver fix is applied.

## 5. §7 assembly — OK (NOTE N3)

- **(a) Re-derived step by step.**
  1. Saddles in R₋ have D < 0.
  2. Every critical point with D ≥ c (a minimum or a degenerate point) lies in R₊ᵒᵘᵗ and maps injectively into
     M̄′ = {E ≤ K}. ∇E ≠ 0 on {E = K}, so this set is the closure of M′.
  3. Two such points are impossible: the open segment lies in the strictly convex M′, where Hess H ≺ 0. That needs M′
     connected (a disc) and negative definiteness at t° ∈ M′; both hold.
  4. Four saddles would persist (as a *finite* set) for a nearby Morse q′ with the same sign pattern. But q′ has at most
     1 extremum, hence at most 3 saddles.
  5. Total ≤ 4, degenerate critical points included.
- **Supporting numerics.**
  - Hess V_q = 2ακ L∘DΨ for mixed signs (κ < 0 in R₋): 120 critical points, to 9e-29; L has rank 2 (`b5_section7.out`).
  - Sard then gives dense Morse charges.
- **(b) Collinear** (`b5b_collinear_fixed.out`, 150 configurations).
  - On-axis equilibria number 0 or 2 (34 and 116 cases), with a far-field scan to 1e7.
  - After moving one charge 1e-3 off the axis, every non-degenerate on-axis equilibrium has a nearby critical point.
  - The perturbed configurations have exactly 2 critical points, with the index count satisfied.
  - (My first version of this test, in `b5_section7.out` part 3, searched blindly in 3D and missed points. It is
    superseded by the planar search in `b5b`.)
- **NOTE N3.** Theorem 2(b) is stated in ℝ³, where off-axis equilibria form circles. In ℝ² they would be isolated
  mirror pairs, so (b) must not be read as covering n = 2, unlike (a). This affects wording only.

## 6. Independent counterexample search — none found

- **Direct search** (`b6a_direct_search.out`; no reduction used to count).
  - 1,200 configurations, each with a minimum forced through my geometric inverse map:
    - 600 random shapes;
    - 600 near-degenerate shapes, in five families: A₁₂ down to 1e-3 rad, A₁ or A₂ within 1e-3 of π, A₃ small, A₁
      small.
  - All critical points in the plane were found from about 1,100 seeds to 80× the triangle size.
  - The 98 configurations that looked suspicious in float64 were re-verified by 50-digit Newton refinement and Hessian
    classification.
  - **Result: (1 extremum, 3 saddles) in all 1,200, in every family.**
  - The 11 "x₀ not recovered" entries are tolerance artifacts (`b6e_x0_check.out`). The unique extremum satisfies
    Q(x_min) ∝ q to 1e-48, so it is the forced minimum; it lies 400–1,900 triangle-sizes away.
- **The one float "2 extrema" case** from my first (crashed) run, reconstructed and verified at 50 digits, has
  (1 extremum, 3 saddles) (`b6d_candidate_check.out`).
- **All charges at once** (`b6b_winding_scan_step3.out`).
  - The number of minima equals the winding of ∇H(∂M′), which is orientation-preserving because Hess H ≺ 0.
  - 941 shapes: a 3° grid with A₁ ≤ A₂, plus 32 near-degenerate shapes (A₁₂ down to 1e-3, A₂ → π, A₃ → 0).
  - **max |winding| = 1** everywhere.
  - Float64 winding 0 appeared at two A₁₂ = 0.001 shapes. At 40 digits the winding is exactly 1 at all 225 interior
    images (`b6c_flagged_mp.out`).

## Files (all in `audit/`)

| Script | Output | Checks |
|---|---|---|
| `b0_common.py` | — | shared mpmath helpers (auditor's own) |
| `b1_sections1to5.py`, `b1b_far_saddle.py` | `.out` files | §§1–5, index count |
| `b2_section61.py` | `.out` | §6.1 |
| `b3_lemmaC.py` | `.out` | Lemma C |
| `b4a_where_small.py` | `.out`, `b4a_main_certified_boxes.npz` | replay of the main chart, keeping its boxes |
| `b4b_ia_tests.py` | `.out` | `ia.py` unit tests and libm ulps |
| `b4c_chart_containment.py` | `.out` | chart formulas against independent references |
| `b4d_coverage_sliver.py`, `b4g_sliver_mpiv.py` | `.out` files | the π-sliver |
| `b4e_mpiv_recert.py`, `b4e_helper_test.py`, `b4f_mpiv_band_corner.py` | `.out` files | rigorous re-certification |
| (their driver, re-run) | `b4_run_cert_*.out` | certificate re-runs |
| `b5_section7.py`, `b5b_collinear_fixed.py` | `.out` files | §7 |
| `b6a_direct_search.py`, `b6b_winding_scan.py`, `b6c_flagged_mp.py`, `b6d_candidate_check.py`, `b6e_x0_check.py` | `.out` files | counterexample search |
