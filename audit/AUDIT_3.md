# AUDIT 3: Theorem 3, every exponent α ∈ [1/4, 2] (the α-sweep certificate)

**Scope.**
- Theorem 3, i.e. Theorem 1.4 and Proposition 12.1 of the merged paper `paper/main.tex` (file time 2026-10-24 15:01).
  The same statements appear as Theorem 1.3 and Proposition 6.1 in `paper_mixed/main.tex` (14:46), and as MIXED.md §8
  "Theorem 3" (14:46).
- Everything Theorem 3 reuses, read for α-generality: MIXED.md §§1–5 and 7, and Part II of the merged paper
  (Sections 7–13), including its references back into Part I (Lemmas 3.2(b) and 4.1, Step 3 of Section 5).
- The certificate: `code/mixed_certF.py`, `code/mixed_alpha_sweep.py`, `code/mixed_certF_sliver.py`,
  `code/sweep_coverage.py`, `code/check_certF_containment.py`, `code/convexity_threshold.py`; `cert/runs.txt`,
  `cert/alpha_sweep*.out` and `.log` (all seven runs finished, the last at 14:44), `cert/sweep_coverage.out`,
  `cert/mixed_certF_sliver.out` (recorded at 15:02, during the audit), `cert/check_certF_containment.out`,
  `cert/convexity_threshold.out`.
- AUDIT_2's notes N1–N3 and its coverage fix.
- For item 8: the abstract, Section 1 (theorem statements) and Part II of `paper/main.tex`.
- Not in scope: Part I and Theorem 1 (AUDIT_1), and the Coulomb certificate (AUDIT_2), except where Theorem 3 reuses
  them.

**Rules followed.**
- All verification code is the auditor's own: `audit/c*.py`, with outputs in `audit/c*.out`.
- The reference `audit/c0_ref.py` shares nothing with `code/` or with AUDIT_2's `audit/b0_common.py`.
  - It evaluates the curvature numerator N = E₁₁E₂² − 2E₁₂E₁E₂ + E₂₂E₁² of {E < K} from exact derivatives of E.
  - It finds boundary points by root finding on E itself, along the line s = const (vertex 3) or t₁ = const (vertex 1).
  - The certificate's normalisations were derived independently: I = g₃³N/(U₁U₂)² with g₃ = sin s sin(s + A₃)/2, and
    J₁/U₂² = g₁³N/U₂² with g₁ = sin t₁ sin(A₁ − t₁)/2. Agreement with the certificate therefore tests the closed forms
    and Lemma C's algebra together.
- Project modules are imported only as objects under test: chart functions, `ia`, and the F-box function.
- The project scripts that were re-run (`mixed_certF.py`, `mixed_certF_sliver.py`, `sweep_coverage.py`) write no files.
  Their output was captured in `audit/`.
- Nothing outside `audit/` was modified, and no process started by others was touched. At most one heavy computation ran
  at a time (each under about 2.5 GB, the largest about 2 minutes).

## Verdict

**No GAP was found in the logic or in the coverage of Theorem 3.**
- **Reduction.** Every step of §§1–5 and 7 holds verbatim for every α > 0, up to the wording point N10. The only
  exponent-dependent facts used are β > 2, b > 1, K > 4, F > 1 and γ = 2α + 2 > 1. These were checked by hand, and
  numerically at 40 digits at α = 1/4, 0.7, 1, 3/2, 2.
- **Formulas.** For general F, the certificate computes the right quantity.
  - The auditor's reference agrees with the paper's closed forms for I and J₁/U₂² to 2e-49 relative at
    F = 13/12, 5/4, 9/5, 13/5.
  - The interval-F chart functions contain the reference value at all **7,164** test points (**0 failures**), including
    near-degenerate shapes.
- **Bookkeeping.**
  - There are 163 certified slices. Each has a log entry with the four charts at "FAILED 0".
  - The union of their F-boxes is [1.083333333333333, 2.600000000000001], with no gap. It contains [13/12, 13/5] in
    exact rational arithmetic.
  - Eight slices were re-run (at least one from every output, including the last, [1.998046875, 2]). All certify with
    box counts identical to the logs, which also pins down (e₀, r₀) = (0.4, 0.05).
- **Sliver.** The strip run certifies when re-run (59,901 boxes, 0 failed), identical to `cert/mixed_certF_sliver.out`,
  which was recorded during the audit. The band and corner charts need no patch.
- **Adversarial checks.** At α = 1/4, 1 and 2, 1,800 configurations were counted directly with no part of the reduction
  used: none has two extrema or more than four critical points. A boundary-curvature search found no non-positive value.
  The large-exponent numbers of Remark 14.1 are reproduced.

**Findings: no GAP, ten NOTEs.**
- One sentence of the abstract overclaims (N9; fix before publication).
- Two concern records: a pointer to the strip run's output, which was itself recorded during the audit (N1), and
  (e₀, r₀) in the logs (N2).
- Six are wording: N3 and N4 (literally false as written, harmless), and N5, N7, N8, N10.
- One is a latent code issue (N6).

**Caveat.**
- The α-range computation (1.3·10⁹ boxes) was not re-certified in rigorous arithmetic. It still rests on the library
  accuracy assumption of `ia.py`, as the paper and MIXED.md say.
- As a robustness check, the first slice and the hardest slice ([1.998046875, 2]) still certify with 4-ulp rounding per
  operation and 16-ulp rounding per library function, with identical box counts (`c11`, §3).
- This is an internal audit, not a referee report.

| Item | Verdict |
|---|---|
| 1. α-generality of the reduction (planarity in ℝⁿ, signed identities, location of minima, gradient structure, injectivity, Hess H ≺ 0 on M′, product representation, M′ = {cosh Δ < β(s)}, Lemma C, index count, degenerate points, part (b)) | **OK** (NOTES N4, N5, N10: wording) |
| 2. `mixed_certF.py` against the formulas and the audited Coulomb code | **OK** (NOTE N6: NaN masking in `ia._ratio`, latent) |
| 3. Sweep and coverage logic, `runs.txt`, logs, re-runs | **OK** (NOTE N2: (e₀, r₀) not logged; NOTE N3: F-box wording) |
| 4. Sliver patch; band and corner charts | **OK** (NOTE N1: pointer to the strip run's record) |
| 5. Own enclosure test with interval F | **OK**: 0 failures in 7,164 points |
| 6. AUDIT_2 follow-ups (N1–N3, π-sliver fix) | **OK** |
| 7. Adversarial search at α = 1/4, 1, 2 | **none found** |
| 8. Paper (`paper/main.tex`: abstract, Section 1, Part II) | **OK, except the overclaim N9**; NOTES N3, N4, N5, N7, N10 (and N8 for MIXED.md) |

---

## 1. α-generality of the reduction — OK (NOTES N4, N5, N10)

**By hand.** Each step was re-derived with α, b = 2α + 1, β = b + 1, c = b/(b+1)², K = 1/c and F = K/2 − 1 symbolic.

| Step | What depends on α | Holds for |
|---|---|---|
| Planarity (Lemma 7.1(1)) | nothing: ∇V = −2α(Σwᵢ(y − pᵢ) + zΣwᵢ) with z ⊥ Π; works for any normal space, so in every ℝⁿ | all α > 0 |
| w = κλ, regions R± (Lemma 7.1(2), Cor. 7.3) | nothing | all α |
| Hess V = 2ακ(βS − I); extremum ⇔ D > c; no extremum with both eigenvalues of βS − I negative | β > 2 | α > 0 |
| Signed identities (P), (D), (E); no maxima; minima in R₊ᵒᵘᵗ | c > 0 only | all α |
| Separated charges, ∇H = ψ_q | b appears as an exponent; algebra is exponent-free | all α |
| Injectivity of x ↦ t | geometric | all α |
| E = 1/D; M̄′ = {E ≤ K} = closure of M′ | K = 1/c | all α |
| e₂(d) = (b+1)² − bE, so Hess H is definite on M′ | polarisation, any b | all α |
| t° ∈ M′ (E(t°) = 4 < K) | K − 4 = (b − 1)²/b > 0 | α > 0 |
| Hess H ≺ 0 at t° | trace = (1 − b)(cot(A₁/2) + cot(A₂/2) − 2cot(A₁₂/2)) < 0 (cot convex and decreasing) | b > 1 |
| Index at infinity, Q_tot = 0 (Lemma 9.2) | degree of (1 − γ/2) − (γ/2)e^{2iφ} is 2 iff γ/2 > \|1 − γ/2\| | γ > 1 |
| Product representation; E < K ⇔ cosh Δ < β(s) | β(s) = [(F−1)cos A₃ − (F+1)cos(2s + A₃)]/2 for any K = 2(F+1) | all F |
| ψ_max well defined, ψ_max ∈ (0, π), ψ_max < A₁₂; boundary smooth (Lemma 10.2) | (2 − (F−1)cos A₃)/(F+1) ∈ (−1, 1) needs F > 1; ψ_max < A₁₂ ⇔ cos A₃ > −1 | F > 1 |
| Vertex-1 form, ψ₁,max < A₁, monotonicity of Λ₃ − Λ₂ in t₂ | ψ₁,max < A₁ ⇔ cos A₁ < 1 | F > 1 |
| Lemma C (Lemma 11.1), mirror symmetry I(1↔2, S → −S) = I | formula for arbitrary β, Λ̃; mirror leaves s, ψ, β, β′, P unchanged | all F |
| Degenerate points / Sard (Hess V_q = 2ακ ℒ∘DΨ) | the map Ψ depends smoothly on α | all α (NOTE N10) |
| (b): Hess on the axis = −2α(Σwᵢ)(I − βe₁e₁ᵀ), non-degenerate iff Σwᵢ ≠ 0; perturbation into a plane | none | all α, n ≥ 3 |

F(α) = (b+1)²/(2b) − 1 is increasing, with F′ = (b² − 1)/b² > 0. F(1/4) = 13/12 and F(2) = 13/5, both exact.

**Numerical check** (`c5_alpha_general.py` → `c5_alpha_general.out`; 40 digits, own geometry):
- 59 critical points were built per α (by q = κλᵢrᵢ^β, at points of R₊ and R₋, a quarter of them near the excentre),
  with 3–13 extrema per α.
- Results:
  - ∇V = 0 to 4e-40.
  - Hess V (numerical second derivatives of V itself) = 2ακ(βS − I) to 3e-39.
  - Type ⇔ D > c: 0 failures. No extremum in R₋. Every extremum is in R₊ᵒᵘᵗ with t ∈ Q.
  - ∇H(t(x)) = ψ_q with b = 2α + 1 to 8e-40.
  - e₂(d) (d by numerical differentiation of φᵢ) = (b+1)² − bE to 1e-37.
  - det Hess H = e₂(d) to 2e-38.
  - Hess H ≺ 0 and E(t°) = 4 at t°: 0 failures.
  - Winding numbers: 1 on every small circle; on the large circle 1 (Q_tot ≠ 0) or 2 (Q_tot = 0) in 98 of 98 cases.
- Planarity in ℝ⁴ (α = 1/4, 1, 2; 120 configurations): all 873 converged critical points lie in the plane, to 7e-28.

**Formulas for general F** (`c1_ref_selftest.py` → `c1_ref_selftest.out`). This covers 1,200 boundary points at
F ∈ {13/12, 5/4, 9/5, 13/5}, on generic, thin (A₁₂ down to 1e-5), flat (A₃ down to 1e-9), lopsided (A₂/A₁ down to 1e-7)
and corner shapes (r down to 1e-5).

| Check | Max relative error |
|---|---|
| Product representation E = 2(cosh Δ + cos A₃)/(sin s sin(s + A₃)), random t ∈ Q | 7e-54 |
| Region formula cosh Δ = β(s) at the boundary point found on E | 2e-54 |
| Vertex-1 region formula cosh(Λ₃ − Λ₂) = β₁(t₁) | 1e-54 |
| I from E's derivatives vs the closed formula for I | 2e-49 |
| J₁/U₂² from E's derivatives vs the closed formula for J₁/U₂² | 5e-49 |
| ψ_max < A₁₂, ψ₁,max < A₁; arc signs | 0 failures |

**NOTE N4 (wording: Theorem 1.4 via Theorem 1.2(b), `paper/main.tex` lines 90 and 105–107; proof line 871; MIXED.md §8).**
- Theorem 1.4 extends part (b), "the critical points off the line form circles about it", to ℝⁿ for n ≥ 3.
- For n ≥ 4, the off-axis critical points form (n − 2)-spheres about the line, not circles, so the sentence is literally
  false there.
- The count is unaffected: these sets are still not isolated and consist of degenerate points.
- Suggested wording: "form (n − 2)-spheres about it (circles when n = 3)".

**NOTE N5 (wording: Lemma 11.1, `paper/main.tex` line 740; MIXED.md line 167).**
- "{Ĝ > 0} is strictly convex near a boundary point iff J > 0 there": only "if" is true. At J = 0 the set can still be
  strictly convex, as {y > x⁴} at 0 shows.
- The proof uses only "if".
- Suggested wording: "the boundary has positive curvature at a point iff J > 0 there".

**NOTE N10 (wording: Section 13, `paper/main.tex` line 865, "Step 3 of Section 5 applies verbatim").**
- Step 3 (lines 426–445) is written for x ∈ T and positive charges: Ψ = (log(Q₁/Q₃), log(Q₂/Q₃)) on T, and
  κ = q₃/Q₃ > 0.
- For mixed signs:
  - Qᵢ/Q₃ < 0, so the logarithms need absolute values;
  - Ψ must be taken on R₊ ∪ R₋;
  - κ = q₃/Q₃ is negative in R₋.
- With these changes the argument goes through (`ℒ` is invertible because x₀ is off the side lines; AUDIT_2 checked
  Hess V_q = 2ακ ℒ∘DΨ for mixed signs numerically). It is not verbatim, though.
- Suggested wording: "Step 3 applies, with log|Qᵢ/Q₃| on R₊ ∪ R₋ and κ of either sign".

## 2. `mixed_certF.py` against the formulas and the Coulomb code — OK (NOTE N6)

**Code comparison** (`c10_code_diff.py` → `c10_code_diff.out`). The three interval-F chart functions differ from the
audited Coulomb functions only by replacing the constants F − 1 = 1/4, F + 1 = 9/4 and q = 1/3 with the intervals
Fm = F − 1, Fp = F + 1 and q = √(Fm/Fp), all formed with `ia` operations, plus cosmetic changes. There are 22, 17 and
9 differing lines, each listed in the output.

**Other code facts** (by reading):
- No Coulomb constant remains: every numeric literal is a structural 0.5, 1, 2 or 4.
- F enters only through `ia.sub(F, one)` and `ia.add(F, one)`, so every F-dependent quantity is rounded outward.
- The chart parametrisations are identical to the Coulomb ones:
  - main: (A₁₂, λ, ρ);
  - band: (e, λ, ρ), with ψ_max = 2 arcsin(q sin(e/2)) and the same e-scalings and mean-value step;
  - corner: (r, θ, ρ₁, ±), with A₁ = ia.PI − r, so it contains π.
- The ζ quadratic is in rationalised form; the arctan branch is valid because A₂ ≤ r₀ < π/2.
- The corner drop rule is the same. After all roundings, a dropped box lies strictly inside a corner region (worst-case
  margin 4e-17 beyond π − r₀), i.e. inside the corner chart's domain or its mirror.
- `bnb`:
  - the initial grid has exact outer edges;
  - neighbouring boxes share the same float edges, and bisection children share the float midpoint;
  - the F-scale floor of 0.05 only changes which side is split;
  - certification requires `isfinite(lo) & lo > 0`;
  - the cap and round limits leave boxes unevaluated only as FAILED.

**Coulomb consistency** (`c3_enclosure.py`, part 2). At the point F-box F = 5/4, the interval-F functions and the
audited Coulomb functions were compared on 1,600 boxes:
- all enclosures overlap;
- lower ends agree to 7e-14 relative.

**NOTE N6 (code, latent: `code/ia.py` lines 163–188).**
- `_ratio` returns 1.0 wherever `y > 0` is false, which includes y = NaN.
- So `sinc`, `asinc` and `asinhc` of a NaN interval return about [1, 1] instead of NaN (`c7_ia_nan.out`, part A). Every
  other `ia` function propagates NaN, or returns a valid global bound.
- This is not reachable in the certificate. The band chart (the only caller) builds every argument of these functions
  from finite box coordinates by finite operations.
- Instrumented band-chart runs on the slices [0.25, 0.3] and [1.998046875, 2] made 6.1 million `sinc` calls and
  0.47 million calls each of `asinc` and `asinhc`, with **0** non-finite inputs.
- Suggested fix: propagate NaN in `_ratio`.

## 3. Sweep, coverage, `runs.txt`, logs, re-runs — OK (NOTES N2, N3)

**Own coverage check** (`c2_coverage.py` → `c2_coverage.out`; it does not import `sweep_coverage.py`).

| Output | Certified slices | Not closed | Log problems | Final line |
|---|---|---|---|---|
| `alpha_sweep.out` | 34 | 8 | 0 | none (stopped at 1.55, as planned) |
| `alpha_sweep_mid.out` | 24 | 0 | 0 | SWEEP … ALL CERTIFIED |
| `alpha_sweep_hi.out` | 2 | 3 | 0 | none (stopped) |
| `alpha_sweep_hi2.out` | 23 | 2 | 0 | none (stopped) |
| `alpha_sweep_t1.out`, `_t2`, `_t3` | 27, 27, 26 | 0 | 0 | SWEEP … ALL CERTIFIED |

- **Logs.** Every counted slice has, in its .log, the identical header, four chart lines with `FAILED 0`, and
  `TOTAL failed: 0 -> CLAIM M CERTIFIED for alpha in [%g, %g]` with matching values.
- **Not-closed slices.** All 13 are re-covered by certified sub-slices.
- **Old 6-decimal outputs.**
  - These are matched uniquely by the auditor's own replay of the slicing.
  - The copy of `code/mixed_alpha_sweep.py` staged at 13:11 prints
    `%.6f` and otherwise equals the current driver; both pass `repr(lo)` to `mixed_certF.py`.
  - The coverage verdict is the same if the 6-decimal text is read literally instead.
- **F-boxes.**
  - The F-boxes (dn(F(lo)), up(F(hi))) were recomputed and agree bit for bit with `mixed_certF._F_of_alpha` and
    `ia.dn`/`ia.up`.
  - Their union is [1.083333333333333, 2.600000000000001], with no gaps. The smallest overlap of consecutive boxes is
    8.9e-16.
  - The union contains [13/12, 13/5] in exact rational arithmetic.
  - The strip run's single F-box [1.083333333333333, 2.600000000000001] also contains [13/12, 13/5].
- **Project script.** `sweep_coverage.py` itself, run on the finished `runs.txt`, prints YES
  (`c2b_project_sweep_coverage.out`, identical to `cert/sweep_coverage.out`).

**Re-runs** (`c4_rerun.py` → `c4_rerun.out`). `mixed_certF.py` was run with e₀ = 0.4, r₀ = 0.05 and exact float end
points; every slice certified, with per-chart counts identical to the log.

| Slice (output) | main | band | corner +1 | corner −1 | time |
|---|---|---|---|---|---|
| [0.25, 0.3] (alpha_sweep) | 192,633 | 4,096 | 2,048 | 2,048 | 1 s |
| [1.0, 1.05] (alpha_sweep) | 7,585,347 | 91,568 | 2,048 | 2,350 | 87 s |
| [1.55, 1.55625] (mid) | 6,260,950 | 120,658 | 2,048 | 43,264 | 69 s |
| [1.7, 1.70625] (hi) | 9,803,591 | 157,476 | 2,048 | 410,656 | 101 s |
| [1.834375, 1.8375] (hi2) | 7,146,327 | 113,404 | 2,048 | 284,672 | 73 s |
| [1.84375, 1.845703125] (t1) | 4,675,242 | 114,900 | 2,048 | 284,672 | 43 s |
| [1.947265625, 1.94921875] (t2) | 6,747,333 | 175,232 | 2,048 | 308,224 | 62 s |
| [1.998046875, 2] (t3) | 8,065,005 | 234,943 | 2,048 | 311,291 | 84 s |

- The computation is deterministic, and the main counts depend on e₀ and the corner counts on r₀. So the identical
  counts show that every output used (e₀, r₀) = (0.4, 0.05), as `runs.txt` states.
- The driver passes the same e₀ and r₀ to every slice of a run.

**Robustness to the libm assumption** (`c11_wide_rounding.py` → `c11_wide_rounding.out`).
- The four charts were re-run with every operation widened by 4 ulps and every library function by 16 ulps.
- Both [0.25, 0.3] and the hardest slice [1.998046875, 2] still certify. The box counts are identical to the standard
  run (192,633 and 8,065,005 main boxes): the wider rounding did not change the outcome of the bisection.
- So on these slices the certificate tolerates library errors of up to about 14 ulps, against the 1.26 ulp that
  AUDIT_2 measured. This is evidence, not a proof, for the other 161 slices.

**Totals** (`c9_log_totals.py` → `c9_log_totals.out`).
- The 163 counted slices hold 1.295·10⁹ certified boxes, with 3.5 h of logged chart time (3.9 h including the
  not-closed attempts).
- Widths range from 2⁻⁹ to 0.05.
- This matches "about 1.3·10⁹ boxes, about four hours".

**NOTE N2 (record).**
- The logs do not record e₀, r₀ or the command line. `mixed_certF.py`'s first output line, which shows the F-box, is
  filtered out by the driver.
- The values are established above by re-running, but the certificate should carry them.
- Suggested fix: log `sys.argv` and the F-box per slice.

**NOTE N3 (wording: `paper/main.tex` line 835 "for all F in a box containing the values of F on the slice"; code
comment `code/mixed_certF.py` line 205 "cover the exact F-range").**
- The 2-ulp widening of the binary64 value of F does not always contain the exact value, because F(α) in binary64 can
  be off by more than 2 ulps (up to about 3.6 ulps here).
- For 3 of the 163 slices, the box misses the exact F-range of its own slice by 1 to 1.6 ulps at one end:

  | Slice | End | Margin |
  |---|---|---|
  | [0.9500000000000003, 1.0000000000000002] | upper | −3.5e-16 |
  | [1.1500000000000004, 1.2000000000000004] | upper | −2.7e-16 |
  | [1.3000000000000005, 1.3250000000000006] | lower | −2.1e-16 |

- This is **harmless**. The coverage argument is made in F-space on the actual boxes (`sweep_coverage.py`, `c2`), and
  the missed values lie in the neighbouring slices' boxes.
- Suggested wording: "for all F in a floating-point box around the values of F on the slice; the union of these boxes is
  then checked to contain [13/12, 13/5]".

## 4. The sliver patch; band and corner charts — OK (NOTE N1)

- **Main chart.**
  - `mixed_certF.py` runs the main chart on A₁₂ ∈ [e₀, math.pi], with math.pi < π.
  - `mixed_certF_sliver.py` covers A₁₂ ∈ [math.pi − 1e-6, up(π)] for F ∈ [dn(F(1/4)), up(F(2))] ⊇ [13/12, 13/5], with
    the same r₀ and drop rule.
  - The two overlap by 1e-6, so together they cover [e₀, π), plus non-triangles above π.
  - The auditor's re-run certifies: 59,901 boxes, 0 failed, in 19 rounds (`c4b_sliver_rerun.out`), identical to
    `cert/mixed_certF_sliver.out`.
- **Band.** e ∈ [0, e₀] and main A₁₂ ∈ [e₀, …] share the same float e₀ = 0.4. There is no gap.
- **Corner.**
  - The corner chart is parametrised by r = π − A₁ ∈ [0, r₀], with exact float edges.
  - A₁ is formed as `ia.PI − r`, where `ia.PI` = [dn(math.pi), up(math.pi)] contains π. So A₁ ∈ [π − r₀, π] is covered
    exactly, and no strip is missing.
  - The mirror corner A₂ ≥ π − r₀ follows by the 1 ↔ 2 symmetry.
- **Shape space.** It is covered by band ∪ main ∪ strip ∪ corner ∪ mirror corner. Main-chart boxes are dropped only
  inside the corners (§2).
- **Lower bound near A₃ → 0.** At α = 1/4 the scale-free minimum of I/A₁₂⁴ sits at A₃ → 0 and tends to +0.00657 (§7).
  So the strip is a genuinely needed but comfortably positive region, consistent with the fast certification.

**NOTE N1 (record; mostly resolved during the audit).**
- When the audit started, the strip run had no output in `cert/`. MIXED.md §8 and Section 12 of the paper ("A separate
  run covers the strip …") cited it without a record.
- `cert/mixed_certF_sliver.out` (15:02) now holds it, and it is identical to the auditor's re-run.
- What remains: MIXED.md §8 (line 294) and `cert/runs.txt` still do not point to that file. The certificate for
  Theorem 3 is complete only together with this run, so the pointer should be added.

## 5. Own enclosure test with interval F — OK

`c3_enclosure.py` → `c3_enclosure.out`.
- **Setup.**
  - The chart functions of `mixed_certF.py` were evaluated on random boxes in all four coordinates, with
    F ∈ [13/12, 13/5].
  - Three families, 150 boxes each per chart: generic (widths 10^U(−7,−1) of each side), point boxes, and degenerate
    boxes.
  - Degenerate boxes push A₁₂ → π (A₃ down to 1e-12), λ → 0 or 1 (down to 1e-9), e → 0 or r → 0 (down to 1e-8),
    θ → 0 or 1, and ρ → ±1 (to 1e-9).
  - The reference (`c0_ref.py`, 50–80 digits) was evaluated at the centre and at 3 random points of each box.

| Chart | Reference points | Containment failures | Non-finite enclosures | Smallest relative slack |
|---|---|---|---|---|
| main | 1,792 | **0** | 2 | 3.6e-14 (point boxes) |
| band (I/e⁴) | 1,800 | **0** | 0 | 2.8e-9 (point boxes) |
| corner + | 1,800 | **0** | 0 | 9.2e-15 |
| corner − | 1,772 | **0** | 7 | 9.5e-15 |

- Non-finite enclosures can never certify, so they are harmless.
- The band chart's point-box enclosures are wide because of the mean-value step, as expected.
- With the project's own test (`cert/check_certF_containment.out`, 3,965 points against `b0_common`), two independent
  references now agree with the interval-F enclosures.

## 6. AUDIT_2 follow-ups — OK

| Item | MIXED.md | Paper (`paper/main.tex`) | Code |
|---|---|---|---|
| N1 (sheet argument; only E(t°) = 4 < K used) | §4 ✓ | proof of Lemma 9.1 ✓ | — |
| N2 (arctan branch; arg form; A₂ ≤ r₀ < π/2) | §6.2 ✓ | Section 12, corner chart ✓ | `J1n_corner` (A₂ = θr ≤ r₀) ✓ |
| N3 (part (b) needs n ≥ 3) | §7 ✓ | Remark 14.2 ✓ | — |
| π-sliver fix | §6.3 ✓ | "double just above π" ✓ | `mixed_cert_all.py` line 76 `float(ia.up(np.pi))` ✓; the interval-F code uses the separate strip run instead (§4) |

## 7. Adversarial search at α = 1/4, 1, 2 — none found

**Direct count** (`c6_direct_search.py` → `c6_direct_search.out`; no part of the reduction is used to count).
- **Configurations.** Each α had 300 forced-minimum configurations and 300 random-charge configurations, in four shape
  families: generic, thin (A₁₂ down to 1e-3), flat (A₃ down to 1e-3) and corner (A₁ up to π − 1e-3).
  - Forced minima: t was drawn in M′ from E's definition, half of them within 5% of the boundary in E. x(t) came from
    the auditor's inverse map (inversion at p₃), and qᵢ = λᵢrᵢ^β.
- **Search.**
  - Damped Newton from about 3,500 seeds per configuration, out to 10⁶ triangle sizes, plus the forced point.
  - Flagged configurations were re-solved at 50 digits; a point counts only if |∇V| < 1e-35·scale and the last Newton
    step is below 1e-30.

| α | Forced: (1 ext, 3 sad) | Random: (0 ext, 2 sad) | Two or more extrema, or more than 4 critical points |
|---|---|---|---|
| 1/4 | 299 of 300 (11 after the 50-digit step), plus 1 (see below) | 300 | **0** |
| 1 | 300 (6 after the 50-digit step) | 300 | **0** |
| 2 | 300 (6 after the 50-digit step) | 300 | **0** |

- **The one α = 1/4 case** (`c6b_flagged_detail.out`, `c6c_flagged_case.py` → `c6c_flagged_case.out`).
  - A thin triangle (A₁ = 1.2e-4, A₂ = 9.0e-4), with the forced minimum 0.4% from the fold (K − E small) and nearly
    cancelling charges (Q_tot/Σ|q| = 2e-8).
  - With binary64 charges the minimum and a nearby saddle annihilate.
  - With the construction redone at 50 digits there are exactly (1 minimum, 3 saddles), index sum −2, found from
    564 seeds.
- The float runs' spurious "many extrema" lists all disappear at 50 digits; they were non-converged points.

**Boundary curvature** (`c8_curvature_search.py` → `c8_curvature_search.out`; reference values, 40 digits).
- 1,000 random (shape, ρ) per α in four families (A₃ down to 1e-8, A₁₂ down to 1e-5, A₂/A₁ down to 1e-6), plus local
  refinement and 300 corner-chart points per α.

| α | Non-positive values | min I/A₁₂⁴ (where) | min J₁/U₂² (corner) |
|---|---|---|---|
| 1/4 | 0 | 0.00657 (A₃ → 0, A₁ ≈ A₂ ≈ 90°, ρ = −1) | 0.666 |
| 1 | 0 | 0.0411 (A₃ → 0) | 5.33 |
| 2 | 0 | 0.0727 (90.1°, 83.7°, 6.2°) | 12.8 |

These agree with MIXED.md's table (+0.0066 at α = 1/4, +0.073 at α = 2).

## 8. Paper (`paper/main.tex`, 15:01: abstract, Section 1, Part II) — OK except N9

- Part II (Sections 7–14) agrees with `paper_mixed/main.tex` apart from the renumbering and the internal references.
- Those references were checked:
  - Lemma 7.2's proof cites Lemmas 3.2(b) and 4.1, which are (P) and (D)/(E) for x ∈ T, where σᵢ = θᵢ;
  - Section 13 cites Step 3 of Section 5 (see N10);
  - Corollary 1.3 cites Theorem 1.1.
- The LaTeX log shows no undefined references.

**Checked against the evidence above.**
- Abstract, Theorem 1.4, Section 7 ("α ∈ [1/4, 2] for Theorem 1.4") and Proposition 12.1.
- The "Other exponents" paragraph:
  - 163 slices, widths 0.05 to 2⁻⁹, narrower for larger α ✓;
  - (e₀, r₀) = (0.4, 0.05) ✓;
  - about 1.3·10⁹ boxes and four hours ✓;
  - the strip run ✓ (record: N1);
  - exact-rational coverage ✓;
  - about 4,000 containment points ✓;
  - the libm disclosure ✓.
- Section 13 for general α ✓.

**Large-exponent remark, spot check** (`c8`, part 2; 50 digits, the auditor's reference). Along A₁ + A₂ = 10⁻³ with
A₂/A₁ = 10⁻⁵, the grid minima over ρ ∈ [0.80, 0.94] reproduce `cert/convexity_threshold.out` to all printed digits.

| α | Grid minimum | Refined minimum (ρ) |
|---|---|---|
| 5.6 | +0.010709 | +0.010645 (0.867) |
| 5.7 | +0.0052544 | +0.005218 (0.867) |
| 5.8 | −0.0003443 | −0.00036083 (0.866) |

At α = 10, the triangle (10.6°, 29.3°, 140.1°) has min I = −0.0434 over its two arcs, so negative curvature is confirmed.

**NOTE N9 (overclaim, fix before publication: abstract, `paper/main.tex` line 38).**
- The abstract says: "Together, these results prove Maxwell's bound (k − 1)² for k = 3 Coulomb charges". The
  introduction (line 51) defines Maxwell's bound for *isolated* equilibrium points.
- For collinear charges of different signs, Theorem 1.2(b) and Corollary 1.3 bound only the non-degenerate equilibria,
  and isolated degenerate ones do occur.
- Explicit example (`c12_collinear_degenerate.py` → `c12_collinear_degenerate.out`):
  - charges (−2, 24, −64) at positions 0, 1, 3 on a line in ℝ³;
  - the point x* = (−1, 0, 0) is critical with Hess V(x*) = 0;
  - the Jacobian of the off-axis equations (Σwᵢ = 0, Σwᵢtᵢ = 0) in (x₁, ρ²) at (−1, 0) has determinant −2.53, so no
    circle of critical points comes near x*;
  - on the axis, f′ vanishes to order exactly 2 at x*.
  So x* is an isolated equilibrium that the theorems do not count.
- The corollary itself is stated correctly ("at most four non-degenerate equilibrium points"). The abstract's sentence
  claims more than is proved for collinear charges.
- Suggested fix: either qualify the sentence ("… for non-collinear charges, and for non-degenerate equilibria of
  collinear ones"), or add a proof that the on-axis equilibria (zeros of f′) of collinear charges number at most four.

**NOTE N7 (wording: Remark 14.1, `paper/main.tex` line 887; MIXED.md line 283).**
- "The range of α … is limited by computing time rather than by the method" overstates the case at both ends.
  - The next paragraph says that below 1/4 the charts must be rescaled, which is a change of method.
  - Near α ≈ 5.75 the margin tends to 0, so no finite computation closes there.
- Suggested wording: "between 2 and about 5.7 the range is limited by computing time".

**NOTE N8 (wording: MIXED.md).**
- §6.1 line 152 (P ≥ 3.5, |Λᵢ|, |Δ| ≤ log 2) and §6.3 line 212 (q = 1/3) are Coulomb values. They are not used by the
  certificate and are false for other α; for example, at F = 13/5 one only has P ≥ 0.8.
- The title and status block of MIXED.md still describe the Coulomb result only.
- Suggested fix: mark these values as Coulomb-only, and mention Theorem 3 and this audit in the status block.

N3, N4, N5 and N10 (above) also concern the paper.

## Files (all in `audit/`)

| Script | Output | Checks |
|---|---|---|
| `c0_ref.py` | — | the auditor's reference (E-derivatives, root finding on E) |
| `c1_ref_selftest.py` | `.out` | product and region formulas, I and J₁/U₂² formulas for general F |
| `c2_coverage.py` | `.out` | slices, logs, replay, F-union in exact arithmetic, F-box wording |
| (project `sweep_coverage.py`, re-run) | `c2b_project_sweep_coverage.out` | the project's coverage verdict |
| `c3_enclosure.py` | `.out` | interval-F enclosures vs reference; Coulomb consistency |
| `c4_rerun.py` | `.out` | 8 slices re-run, counts vs logs |
| (project `mixed_certF_sliver.py`, re-run) | `c4b_sliver_rerun.out` | the strip near A₁₂ = π |
| `c5_alpha_general.py` | `.out` | exponent-dependent steps at five α |
| `c6_direct_search.py`, `c6c_flagged_case.py` | `c6_direct_search.out`, `c6b_flagged_detail.out`, `c6c_flagged_case.out` | direct counterexample search |
| `c7_ia_nan.py` | `.out` | NaN handling in `ia.py` |
| `c8_curvature_search.py` | `.out` | curvature search; large-exponent spot check |
| `c9_log_totals.py` | `.out` | box and time totals |
| `c10_code_diff.py` | `.out` | interval-F code vs Coulomb code |
| `c11_wide_rounding.py` | `.out` | two slices with 4/16-ulp rounding |
| `c12_collinear_degenerate.py` | `.out` | an isolated degenerate collinear equilibrium (N9) |
