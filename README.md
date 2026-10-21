# Recurrence-Free Survival in Breast Cancer: A Cox Model with External Validation

**[Read the rendered report](https://suparnaperumal11.github.io/Survival-analysis-breast-cancer/)** — the full analysis with figures and diagnostics, no R required.

A survival analysis of recurrence-free survival (RFS) in breast cancer, built on two published cohorts, to answer one question:

> **What can a clinician honestly tell a patient about their recurrence risk at diagnosis — and where does a simple, auditable model stop being trustworthy?**

The hazard ratios are not the deliverable. The boundary is.

---

## 1. The decision this supports

A patient asks, once surgery is done and pathology is back, what the next ten years look like. The honest answer has two halves: what the measured prognostic factors actually support, and the point past which the model answering the question is no longer reliable enough to quote. Both halves are the output here.

Age, tumour size, nodal involvement, grade and treatment are all known at that moment. This analysis asks how much of the subsequent risk of recurrence or death those factors explain in a straightforward, auditable Cox model — and then spends most of its length establishing where that model's honesty breaks down. The proportional hazards assumption is tested and reported rather than assumed; four of seven covariates violate it and are re-fit accordingly; discrimination, calibration and prediction error are measured on data the model was never fit to.

Three of those boundaries are firm enough to state up front. The model ranks patients only moderately well (concordance 0.674, optimism-corrected 0.6725), so it separates higher-risk from lower-risk profiles without resolving individuals. Four of its predictors have no single hazard ratio to quote, because their effect on the hazard measurably changes over follow-up. And both cohorts predate current treatment practice by three decades. What that supports is a conversation about relative risk — nodal burden dominates the profile, and by roughly this much. What it does not support is a number for one patient's ten-year recurrence probability.

This is a methods project on published cohorts, not a validated clinical tool. Section 9 states what would have to change for it to become one.

**My role.** An independent methods project: the question, the endpoint construction, the modelling, the diagnostics and this write-up are all mine.

## 2. Who would use this estimate, and at what point in care

A model like this — hazard ratios for baseline prognostic factors, checked for proportional hazards violations, internally and externally validated — is the kind of analysis a hospital's clinical data or biostatistics team would build to support **prognostic counselling at the point of diagnosis and treatment planning**: after surgery, once pathology (tumour size, grade, nodal status) is known, before adjuvant treatment decisions are finalised. It is not a real-time risk calculator and it is not this analysis's job to be one — the model here is fit on data from the 1970s–90s and is not a calibrated instrument ready for a 2026 Singaporean patient. Section 9 says explicitly what stands between this and that.

## 3. Cohorts and methods

**Development cohort**: `rotterdam` (R `survival` package), n = 2982, breast cancer surgeries between 1978 and 1993 (Foekens et al., 2000).
**External validation cohort**: `gbsg` (R `survival` package), n = 686, drawn from the GBSG-2 trial (Schumacher et al., 1994) — which enrolled **node-positive patients only**, a case-mix difference from `rotterdam` that matters for how the validation result should be read (Section 4).

**Endpoint**: recurrence-free survival — time to whichever of recurrence or death comes first, in years. `rotterdam` records recurrence and death as two separate time/event pairs; naively taking `pmin()` of the two time columns silently mis-times the event for a small subset of patients where the last disease-free assessment predates the actual death date. The endpoint was constructed to route on which event was observed, not on which time value was smaller — see `analysis.Rmd` for the full check (43 of 2982 patients affected by this if handled naively).

**Covariates**: age, menopausal status, tumour size (categorical, ≤20/20–50/>50mm), nodes (log-transformed after a martingale-residual check showed the raw count's effect on log-hazard saturates rather than staying linear), grade, hormonal therapy, chemotherapy.

**Model**: Cox proportional hazards (Cox, 1972), adjusted for all baseline covariates simultaneously. Proportional hazards was checked with `cox.zph()` — printed test and residual plots both — and violated for four of the seven covariates; those four were re-fit with time-varying coefficients (`tt()`) rather than reported as if the assumption held. R's standard validation tools (`concordance()`, `survfit()`, `pec()`) cannot process a `coxph` object with `tt()` terms, so a separate flat (time-fixed) model is used for calibration, Brier score, and bootstrap validation, while the time-varying model's own discrimination is assessed at fixed follow-up horizons using its own fitted coefficients. External validation uses a version of the flat model without `chemo` (not recorded in `gbsg`).

Full derivations, diagnostic checks, and code are in `analysis.Rmd` / `analysis.html`.

## 4. Key results

- **Unadjusted**, patients who received hormonal therapy had *worse* recurrence-free survival (median 4.78 vs. 7.15 years, log-rank p = 3×10⁻⁴). **Adjusted** for age, size, nodes, grade, and chemotherapy, the direction reverses: hazard ratio 0.68 (95% CI 0.58–0.79, p = 1.1×10⁻⁶) — a 32% lower hazard. This is confounding by indication, not a treatment paradox: hormonal therapy was preferentially given to higher-risk patients, and the unadjusted comparison reflects that risk imbalance more than the treatment itself. Even the adjusted estimate only controls for the confounders measured here — see Section 8.
- **Strongest predictor**: nodal burden (log(nodes+1)), hazard ratio 1.76 per unit (95% CI 1.66–1.86, p < 10⁻⁸⁰) in the time-fixed model — though its effect is not actually time-fixed (Section 6).
- **Tumour size** and **grade 3** (vs. grade 2) both independently elevated hazard (size >50mm: HR 1.61, 95% CI 1.37–1.90; grade 3: HR 1.40, 95% CI 1.25–1.58).
- **Age and menopausal status** showed no independent association in the time-fixed model (both p > 0.8) — but age's effect is not actually time-fixed either (Section 6).
- **RMST at 10 years**, hormonal therapy vs. none (naive, unadjusted): −0.77 years (95% CI −1.20 to −0.34, p < 0.001) — hormonal therapy patients had fewer recurrence-free years within the horizon. This carries the identical confounding-by-indication caveat as the unadjusted KM comparison above; it is not an adjusted "life-years gained from treatment" estimate, and shouldn't be read as one.
- **Cause-specific sensitivity check**: decomposing the composite RFS endpoint into recurrence-specific and death-specific cause-specific hazards, hormonal therapy's association is in the same direction for both (HR 0.69 for recurrence, HR 0.75 for death) — the composite result isn't an artefact of one component alone.

## 5. Model validation

**Internal validation.** Apparent concordance for the flat model is 0.674 (se 0.007). A bootstrap optimism-correction procedure (B = 200 resamples) puts the corrected concordance at 0.6725 — optimism of only 0.0013, consistent with a model this size (7 predictors, n = 2982) having little room to overfit its own training data.

**Discrimination over time.** Time-dependent AUC for the flat model is fairly stable across follow-up: 75.2% at 1 year, 71.8% at 5 years, 69.8% at 10 years. The time-varying model's own discrimination — assessed with a landmark risk score built from its fitted coefficients at each horizon, since standard tools can't evaluate a `tt()` model directly — is **higher at 1 year (76.4%), similar at 5 years (71.9%), and lower at 10 years (63.0%)** than the flat model's. This decline at the longest horizon is reported as found, not explained away; one plausible cause is that the linear-in-time coefficient specification extrapolates less reliably the further it's pushed past where most of the data sit.

**Calibration.** The flat model's calibration slope on its own development data is 1.00 — expected by construction of the Cox partial-likelihood fit, and not read as evidence of good calibration on its own. The external calibration slope, calculated on `gbsg` (data the model never saw), is 1.13 (se 0.12) — consistent with reasonable calibration, though the uncertainty is wide enough that this shouldn't be oversold as precise. No calibration or Brier/IBS figure exists for the time-varying model, for the same tooling reason as above.

**Prediction error.** Brier score and integrated Brier score (IBS) both favour the flat Cox model over a covariate-free null model at every follow-up horizon checked (e.g., at 5 years: 20.7% vs. 24.5%; IBS 0.179 vs. 0.205 overall).

**External validation.** External concordance on `gbsg`, using the flat model's coefficients frozen from `rotterdam` (not refit), is 0.664 (se 0.016) — a small drop (~0.01) from internal concordance. External time-dependent AUC is 73.3% at 1 year and 71.2% at 5 years. `gbsg`'s maximum follow-up is 7.3 years, so none of these external figures extend to 10 years — that horizon simply isn't estimable from this cohort, and no external validation for the time-varying model was attempted (Section 6 explains why). Given the node-positive-only case mix in `gbsg`, all of this should be read as "holds up on a higher-risk subgroup," not "holds up on an independent general population."

## 6. What I checked and what failed

`cox.zph()` was run on the full adjusted model and rejected proportional hazards globally (p = 1.1×10⁻¹²). Read against the residual plots (not p-values alone — with n≈3000, `cox.zph` flags deviations too small to matter as often as it flags ones that do):

- **Real, practically meaningful violations**: `size`, `log(nodes+1)`, and `grade` all show their effect on hazard visibly declining over follow-up — strong early on, fading over years. `age` looked trivial on the plot (a near-flat line) but wasn't: quantifying it, the per-year age hazard ratio drifts from ≈0.98 early in follow-up to ≈1.05 by year 14, which compounds to an 18%-lower-hazard-to-68%-higher-hazard reversal for two patients 10 years apart in age — a real, clinically plausible pattern (younger patients: more aggressive early recurrence risk; older patients: rising competing-cause mortality later), not noise.
- **No meaningful violation**: `meno`, `hormon`, `chemo` — flat residual trends, consistent with the PH assumption.

The fix used was time-varying coefficients (`tt()` in `coxph`) for the four violators, not stratification — stratification would have dropped the hazard ratios for `size` and `grade`, which this analysis needs reported. A nested likelihood-ratio test confirmed the added time-interaction terms are warranted (χ² = 62.3, df = 5, p = 4.0×10⁻¹²), not overfitting. The cost: four of the model's eight predictors no longer have a single hazard ratio to quote — they're reported at 1/5/10-year horizons instead — and one of the four (`grade`'s time-interaction term) didn't reach significance on its own even though the joint test and the visual diagnostic supported including it.

A further, harder limitation surfaced once validation began: `concordance()`, `survfit()`, and `pec()` all refuse a `coxph` object containing `tt()` terms — a documented limitation of the R `survival` package, confirmed directly rather than assumed. That's why the time-varying model's discrimination is assessed separately (landmark scores at fixed horizons, Section 5) rather than with the same bootstrap/calibration/Brier machinery used for the flat model — there is currently no equivalent optimism-corrected, calibration, or Brier figure for the time-varying model in this analysis, and that gap is reported rather than papered over.

## 7. What I chose not to build, and why

**Time-varying coefficients instead of stratification.** Four of the seven covariates — `age`, `size`, `log(nodes+1)`, `grade` — violate proportional hazards. Two standard fixes exist: stratification, which fits a separate baseline hazard per level of the offending variable, and time-varying coefficients, which keep a hazard ratio but let it be a function of follow-up time. Stratification is the more conservative choice and it was rejected. `size` and `grade` are exactly the two variables a prognostic conversation needs a reported hazard ratio for, and stratifying a variable discards its hazard ratio entirely. `age` and `log(nodes+1)` are continuous, so stratifying them would have meant binning them first — an arbitrary cut made to satisfy an assumption, not because the data suggested one.

What the trade bought: a documented, interpretable effect for all four, and one finding stratification would have buried. `age`'s per-year hazard ratio drifts from ≈0.98 early in follow-up to ≈1.05 by year 14 — younger patients carrying more aggressive early recurrence risk, older patients rising competing-cause mortality later. Stratification would have absorbed that into the baseline hazard instead of reporting it.

What it cost, which is the larger half: `concordance()`, `survfit()` and `pec()` all refuse a `coxph` object containing `tt()` terms. So the model this analysis considers correct is the one that cannot be put through standard validation, and a separate time-fixed model carries the calibration, Brier and bootstrap results instead. Four predictors have no single hazard ratio to quote — they are reported at 1/5/10-year horizons — and `grade`'s time-interaction term did not reach significance on its own, even though the joint likelihood-ratio test and the residual plot both supported keeping it.

**No Fine-Gray / competing-risks model.** Death without recurrence is a competing event, and the composite RFS endpoint absorbs it rather than modelling it. The sensitivity check in Section 4 uses cause-specific hazards instead. The reason is that a Fine-Gray model answers a different question — a covariate's effect on the absolute, population-level cumulative incidence of one event type when a competing event can preclude it. That is the right tool for absolute-risk prediction and the wrong one for a prognostic-factor study, which this has been throughout; cause-specific hazards are the direct extension of the Cox models used everywhere else here.

The consequence is a question this analysis cannot answer: the absolute probability of recurrence as distinct from death. The cause-specific check establishes that hormonal therapy's association points the same way for both components of the endpoint (HR 0.69 recurrence, HR 0.75 death), so the composite result is not an artefact of one component. It does not quantify how much of the composite hazard ratio comes from either, and shouldn't be read as doing so. A patient asking "how likely is it to come back, specifically" needs the Fine-Gray model that was not fit.

**No causal claim about the treatment effect.** The most striking number in the analysis is a reversal: hormonal therapy looks harmful unadjusted (median RFS 4.78 vs. 7.15 years, log-rank p = 3×10⁻⁴) and protective adjusted (HR 0.68, 95% CI 0.58–0.79). That is confounding by indication — hormonal therapy was preferentially given to higher-risk patients — recovered by adjustment. It is an adjusted association, not a treatment-effect estimate, and no attempt was made to turn it into one.

The reason is that nothing in the design supports it. Treatment assignment in `rotterdam` is observational and non-randomised, and a Cox model adjusts for measured confounders only. Receptor status and comorbidity both plausibly drove who was treated, and neither is in the model. A propensity-score or marginal-structural approach would have produced a causal-looking estimate resting on exactly the same unverifiable assumption of no unmeasured confounding, in heavier machinery — the inferential problem here is the data, not the estimator.

The consequence is that every treatment number in this repository is associational. The RMST figure needs the most care, because "recurrence-free life-years" invites a causal reading: it is unadjusted, it carries the same confounding-by-indication caveat as the raw Kaplan-Meier comparison, and it is reported that way rather than quietly promoted to a benefit estimate.

**No receptor status — a scope decision, not a data limitation.** `rotterdam` and `gbsg` both carry `er` and `pgr` as continuous assay values (fmol/l), complete for all 2,982 development patients. The model does not use them. The covariate set was fixed to the baseline clinical and treatment factors shared across both cohorts, and receptor concentrations were left out of it. Of everything on this list, that is the choice I would revisit first.

What it rules out: hormonal therapy's benefit is well established to depend on ER/PR status, so the 0.68 hazard ratio is averaging over patients for whom the treatment mechanism does not apply at all. No receptor-stratified or receptor-interacted estimate is available, and the unmeasured-confounding caveat on the treatment association is partly self-inflicted — one of the named confounders is sitting in the dataset, unused. HER2 is genuinely absent from both cohorts, so a modern three-marker subtype model is out of reach from this data regardless of scope.

## 8. Limitations

- **Non-informative censoring is assumed and cannot be tested from these data.** If patients were lost to follow-up for reasons related to their unobserved prognosis (e.g., moved because they were too sick to continue treatment locally), the survival estimates here would be biased in a way this analysis cannot detect or correct.
- **Both cohorts predate current breast cancer treatment practice by three decades.** `rotterdam` spans surgeries from 1978–1993; treatment protocols, staging technology, and hormonal/chemotherapy regimens have changed substantially since. The hazard ratios reported here describe risk relationships under 1980s–90s-era treatment, not current standard of care.
- **No Singapore cohort.** Both `rotterdam` (Netherlands) and `gbsg` (Germany) are Western European cohorts. Nothing here speaks to whether these risk relationships hold in a Singaporean population, whose age distribution, tumour biology patterns, and treatment access may all differ.
- **No cost data exists in either cohort.** RMST here is presented purely as an effectiveness measure (recurrence-free life-years within a horizon). This is not a cost-effectiveness or economic evaluation, and nothing in this analysis should be read as implying one.
- **The external validation is not a clean transportability test.** `gbsg` enrolls node-positive patients only, with follow-up capped at 7.3 years; `rotterdam` does not have either restriction. A small drop in concordance from internal to external (0.674 → 0.664) is a reasonable sign, but it's evidence the model discriminates within a higher-risk subgroup over a shorter horizon, not evidence it generalizes to an unselected population at 10 years.
- **The time-varying model is validated less completely than the flat model, by necessity, not choice.** R's standard validation tools cannot process it; its own discrimination is checked at fixed horizons using its own coefficients, but no bootstrap-corrected, calibrated, or externally validated figure exists for it.
- **The competing-risk sensitivity check uses cause-specific, not Fine-Gray, hazards.** It answers whether `hormon` and `log(nodes+1)`'s associations look consistent across the recurrence and death components of RFS — it does not estimate either component's absolute, population-level cumulative incidence, which a Fine-Gray model would be needed for and was not fit here.
- **Adjustment is only as good as what was measured.** The hormonal-therapy hazard ratio (0.68) accounts for the confounders available in `rotterdam` — it cannot rule out unmeasured confounders (e.g., hormone receptor status, comorbidities) that also influenced who received treatment.

## 9. What would be needed to use this for real

- **A contemporary cohort**, ideally from the target population (Singapore), with modern staging (including receptor status — ER/PR/HER2 — which neither cohort's model here uses despite both containing `pgr`/`er` fields) and current treatment protocols.
- **Calibration and discrimination assessed on a cohort that isn't restricted to a single risk subgroup and to a 7.3-year horizon** — the calibration and time-dependent AUC work here (Section 5) is a start, not a substitute for testing on an unselected, longer-follow-up population.
- **Recalibration or refitting for current practice**, given the 30-year gap between these cohorts and current treatment.
- **A causally-motivated comparison for the hormonal-therapy effect** (e.g., propensity-score methods, marginal structural models) if the goal were an actual "treatment effect" estimate rather than an adjusted association — the Cox-adjusted HR and the RMST here are both associational, not causal, despite the RMST being framed in "life-years" language that invites a causal reading.
- **Receptor-status stratification.** Hormonal therapy's benefit is well established to depend heavily on ER/PR status; a model that doesn't stratify or interact on receptor status is averaging over patients for whom the treatment mechanism doesn't even apply.
- **Cost data**, if the eventual use case is a resource-allocation or cost-effectiveness decision rather than a purely clinical-risk one — nothing here currently supports that.

## 10. Selected methodological references

- Cox, D.R. (1972) — the proportional hazards model underlying every Cox fit in this analysis.
- Grambsch & Therneau (1994) — the scaled Schoenfeld residual test (`cox.zph()`) used for the proportional-hazards check (Section 6).
- Harrell, Lee & Mark (1996) and Harrell (2015) — the bootstrap optimism-correction procedure used for internal validation (Section 5).
- Uno et al. (2014) — the concordance (C-index) statistic used throughout for discrimination.
- Blanche, Dartigues & Jacqmin-Gadda (2013) — the time-dependent AUC methodology (`timeROC`) used for both the flat and time-varying models (Section 5).
- Graf et al. (1999) — the Brier score / integrated Brier score for censored survival data (Section 5).
- Andersen et al. (2012) — restricted mean survival time (RMST) as an effectiveness summary (Section 4).
- Royston & Altman (2013) — principles for external validation of a Cox prognostic model, which the `rotterdam`/`gbsg` external validation here follows (Section 5).
- Foekens et al. (2000) — the original cohort study underlying the `rotterdam` development data (Section 3).
- Schumacher et al. (1994) — the GBSG-2 trial underlying the `gbsg` external validation cohort, including its node-positive-only enrollment (Section 3).

## 11. Files in this project

- **`analysis.Rmd`** — the analysis itself: source code, methodological reasoning, and results, structured to be read top to bottom.
- **`analysis.html`** — the rendered report (knit from `analysis.Rmd`); open this to read the analysis with figures and tables inline without running R.
- **`index.html`** — redirect to `analysis.html`, so the GitHub Pages root URL opens the report.
- **`sessioninfo.txt`** — R and package versions from the session that produced `analysis.html`.
- **`README.md`** — this file.

## 12. Reproducing this analysis

Requires R with the `survival`, `survminer`, `dplyr`, `tidyr`, `ggplot2`, `gtsummary`, `broom`, `survRM2`, `riskRegression`, `pec`, and `timeROC` packages installed, plus Pandoc (bundled with RStudio, or installed separately) for HTML rendering. Exact package versions are recorded in [`sessioninfo.txt`](sessioninfo.txt) — `sessionInfo()` output from the session that knitted the committed `analysis.html`, so it describes what actually produced the numbers above rather than a later machine. The analysis was run under R 4.5.3 with `survival` 3.8-6, `riskRegression` 2026.03.11, `pec` 2025.06.24 and `timeROC` 0.4.1. Dependencies are not pinned with `renv`: the final chunk of `analysis.Rmd` calls `sessionInfo()`, so every rendered report carries its own provenance, and a lockfile built now would describe the current library rather than that run. Running substantially newer versions of `riskRegression`, `pec` or `timeROC` could change output formatting, though the underlying methods are stable.

From an R session in this directory:

```r
rmarkdown::render("analysis.Rmd")
```

This regenerates `analysis.html` from scratch, including all model fits, diagnostics, and figures. Both `rotterdam` and `gbsg` are loaded directly from the `survival` package — no external data files are needed.

---

See `analysis.Rmd` (source) / `analysis.html` (rendered) for the full derivation, code, and diagnostic plots underlying every number above.
