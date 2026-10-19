# CHANGELOG — lowvol_v4, revision 2 (relative to the delivered v4 of 7 October 2026)

Revision 2 implements an internal review of the delivered v4 written from the lecturer's perspective (open
inference questions, missing robustness checks, unclear or overstated statements, implementation details). It was
made in this folder; the delivered v4 is preserved unchanged in `../lowvol_v4.zip` (md5
71134d361c0f4ba63a899c019e95f1ce). `lowvol_v2` and `lowvol_v3` were not modified (all 60 files identical,
`logs/md5_v2v3_r2.txt` = `logs/md5_v2v3_at_delivery.txt`). The database was opened read-only and not copied. No
GitHub actions. Status words as in the v4 part below.

## R1. New computations (`build.R`, `R/functions.R`) — implemented and executed in the delivered run

| # | Review point | Change | Result and check |
|---|---|---|---|
| R1 | The valuation spread is persistent, so the OLS predictive slope is biased (Stambaugh 1999); a market-state confound was not controlled. | Reduced-bias estimator of Amihud and Hurvich (2004) on the non-overlapping December observations (`pred_ah`); predictive regression with the trailing 12-month market excess return as control. | 1968–2012 (45 years): reduced-bias t 1.47 (L–H) and 2.95 (BN); innovation correlation −0.44 / −0.21; the market control lowers the Newey-West t from 2.02 and 3.30 to 1.79 and 3.09 (spread–market correlation 0.33). 2013–2024 (12 years): 0.45 and −0.97. Claim checks. |
| R2 | GFKSW's delisting rule (−30% NYSE/AMEX, −55% NASDAQ) was described but not tested. | Variant `gfksw` of the holding returns: stocks dropped by their exchange (CIZ delisting action GDR on the active security record) get the GFKSW return after the last trading day. | 25 held stock-months 1929–2024, 6 in 1968–2012, one of them in an extreme beta quintile; mean L–H 1968–2012 unchanged within 0.01 bp, state returns within 1 pp (Appendix Table 10). Claim checks. |
| R3 | Sensitivity of the state classification (middle states). | Two variants added to `state_sens`: value-weighted B/P averages within the quintiles and calendar-year book equity. | 50% and 12% of the months change state; one-year L–H returns of states 2–4 move by at most 2.4 pp across all five variants; no variant reproduces GFKSW's middle states (Appendix Table 12). Claim checks. |
| R4 | Concentration of the long-only low-beta quintile. | Quintile with weights capped at 5% (`cap_weights`), compatible with the UCITS diversification rule. | Effective N 44 (Q1: 29), net Sharpe ratio 0.52 (Q1: 0.50). |
| R5 | Better covariance estimation with daily data (course material on nonlinear shrinkage). | Long-only global minimum variance from the daily returns of the last 12 months with the analytical nonlinear shrinkage of Ledoit and Wolf (2020) (`nl_shrink`, incl. p > n; `minvar_weights_nl`, 3% cap). | Solved in all 697 months (on average 550 eligible stocks; unit tests on the server). Volatility 11.5% (lowest), turnover 432% p.a., Sharpe ratio 0.53 gross and 0.49 net, break-even 23 bp: lowest realized volatility among the portfolios considered (frequency, window, estimator and stock selection differ at the same time, so this does not isolate the estimator), but after costs no better than the simple low-beta quintile (0.50). |
| R6 | The 25% beta-neutral portfolio has a realised beta of 0.27. | BN with betas shrunk towards one (0.6β + 0.4, Frazzini and Pedersen 2014), short ratio reset monthly. | Average short ratio 48%, realised beta −0.06; net Sharpe ratio 0.26 (1968–2024); 2013–2024 1.7% p.a. (t 0.61) (Appendix Table 23). |
| R7 | Multiple testing. | Holm adjustment of the long-only comparisons now over 8 candidates (v4: 6). | Only the Holm p-values of existing rows change. |

Effect on existing results: 55 of 75 result components are identical to the delivered v4; in the other 20, every
pre-existing row and column is identical (new rows, columns or variants only), except the Holm p-values of
`LO_post`, `LO_post_net` and `LO_pub_net` (larger test family).

## R2. Text, tables and references (`paper.qmd`, `refs.bib`)

* Introduction: explicit research question and three contributions; no references to lecture slides in the text.
* Literature: Novy-Marx and Velikov (2022), Asness et al. (2017), Ledoit and Wolf (2020) added.
* Section 3: GFKSW's delisting rule described and tested; one-month (ours) versus three-month T-bills (GFKSW);
  BN ex ante and shrunk; the logic of the circular block bootstrap explained; predictive regression with the
  Stambaugh bias correction.
* Section 4.1: the BN excess return (0.38% per month, t 2.87) is named as the number comparable with L–H;
  realised betas of all long-short constructions.
* Section 4.2: the extremes replicate, the middle states do not (cut-offs, next-month signs, classification
  variants); the contrast is borderline and rests on the technology bust and the 2009 rebound; predictive
  regression with reduced-bias t and market control.
* Table 3: our alphas with bootstrap t next to GFKSW's published alphas and t; our GFKSW-style t in the new
  Appendix Table 16.
* Section 5: post-2012 losses decomposed into beta and alpha; the real-time timing table moved to the main text.
* Sections 6–8 renamed Implementation I–III. Table 6 now has nine portfolios (Q1, buffer, quarterly, capped,
  daily-beta quintile, minimum variance monthly and daily, 1/N, cap-weighted) with information ratios;
  t-statistics and standard errors in separate rows so the table needs less scaling. New Appendix Table 22
  (diversification and liquidity); Appendix Table 23 adds shrunk-beta BN, t-statistics in separate rows.
* Corrections after the first render of revision 2: unrendered minus sign in the state cut-offs; "moving-block" →
  "circular block" bootstrap (Table 2 note, Table 8); delisting counts refer to 1968–2012; minimum variance
  (daily): "lowest risk but a lower net Sharpe ratio than the low-beta quintile" instead of "loses its advantage"
  (its net Sharpe ratio stays above the benchmark's); "65% of the zero-cost portfolio's CAPM alpha" instead of
  "of the premium"; timing rule "did not significantly improve"; shrunk-beta BN "nearly" beta-neutral; Table 7
  row label "Min. variance (monthly)"; Table 22 no longer exceeds the text width.
* Corrections after an independent check of the rendered PDF (a separate reviewer that had not written the text):
  Table 5 (timing) and Table 7 exceeded the right margin → t-statistics in separate rows and scaling (Table 5),
  scaling (Table 7); Table 10 headers on two lines; eligible stocks reported for both minimum-variance
  portfolios (517 monthly, 550 daily); shrunk-beta BN beta given for 1968–2012 and 1968–2024; "December 2012,
  when GFKSW's sample ends"; Figure 4 windows start in June 1932; Table 21 note states that the valuation effect
  uses 1968–2024; the market control is compared with the Newey-West t (2.02/3.30 → 1.79/3.09), not with the
  reduced-bias t; "robust predictive power" replaced by "predicts ... significantly in sample"; abstract,
  introduction and conclusion: the contrasts of both portfolios are borderline, the replication of the extremes
  refers to the one-year L–H returns; over-rejection of the baseline bootstrap for the state-5 means stated and
  the volatility-preserving t of the L–H state-5 mean reported (not significant); post-2012 predictive slopes
  reported (L–H similar size, BN negative) instead of "no power"; December-book-value restriction affects the
  middle quintiles; cheapest states "mainly" in the listed years; benchmarking explanation "consistent with"
  instead of "show why"; daily betas better forecasts "from long enough windows"; daily minimum variance "no
  higher" net Sharpe ratio than Q1 (0.49 vs 0.50); highest net Sharpe ratios given to three decimals; variance
  ratios to three decimals with the reason they differ from squared volatility ratios. New claim checks for each.
* Appendix C: AI statement describing the actual use; `repo_url` placeholder for the group's GitHub repository.
* `refs.bib`: politis1992, stambaugh1999, amihud2004, ledoit2020, novymarx2022, asness2017 (43(5), 72–87),
  ucits2009 added; blitz2014 as working paper; URL field for the French data library.
* Author list: on 7 October 2026 (21:31 local time) `paper.qmd` on the server was edited outside this revision
  and Gjergji Emini removed from the authors; revision 2 kept the author list as found; restored in revision 3 (R5).
* New claim checks for every new statement (rendering stops otherwise).

## R3. Delivered run of revision 2

| Step | Status | Evidence |
|---|---|---|
| Fresh build: `cache/` and `results/res.rds` deleted, then `Rscript build.R` (`logs/final_run.sh`) | executed, exit code 0 | `logs/build.log`; `logs/build_resources.txt`: started 2026-10-07 20:16:52 UTC, 12 min 54 s, peak resident memory 7.4 GB in the largest R process; stock panel built in the same run. A first attempt stopped in `add_bp` (column masking); the following development run (`logs/dev_run_r2.log`) preceded two further fixes (GFKSW variant on the active security records, market control) and is not the delivered run |
| Internal checks in `build.R` | executed and passed | the build stops on any failure; it completed |
| Render (`quarto render paper.qmd --to pdf`) | executed 2026-10-09 (final render after the independent check), all claim checks passed | `logs/render.log`; 0 unresolved references in the PDF |
| Page budget and layout | validated (text extraction of all pages; word positions: no text beyond the right margin, max. x = 529.4 pt; visual check of Tables 5, 6, 7, 10, 22, 23) | 27 pages: main text pp. 1–15 (the conclusion ends in the upper half of p. 15, incl. Figures 1–3 and Tables 1–7), references pp. 15–18, appendix pp. 19–27 (state after revision 3: see R5) |
| Independent check | executed | a separate reviewer checked about 150 numbers in the text against the tables and all quoted GFKSW values against the original article, and reported two misattributed numbers, several overstated statements and layout issues; all were corrected as listed in R2 |
| `lowvol_v2`, `lowvol_v3` unchanged | validated | `logs/md5_v2v3_r2.txt` identical to `logs/md5_v2v3_at_delivery.txt` (60 files) |

## R4. Open points

* The group must enter its GitHub URL (`repo_url`), confirm the AI statement, and check the university's rules on
  AI use and on WRDS/CRSP data before submission.
* A turnover penalty for the daily minimum-variance portfolio was not implemented.
* The discrepancy in GFKSW's middle states remains unexplained (CRSP vintage, unpublished details).

## R5. Revision 3 (9 October 2026): five targeted corrections

| # | Correction | Change | Validation |
|---|---|---|---|
| 1 | Significance wording | Abstract, introduction, Section 4.2 and conclusion: "The cheap-minus-expensive contrasts are not robustly significant at the 5% level." (replaces "borderline"); the bias-corrected BN predictive regression (t = 2.95) is kept and described as a test of a different hypothesis (linear relation between spread and next-year return, not the difference between the extreme states). | Sentence present in all four places; "borderline" no longer in the PDF. |
| 2 | Table 7: monthly market excess return | `build.R`, section 7: `prod(1 + mkt) - 1` replaced by `prod(1 + mkt + rf) - prod(1 + rf)`; stop if the market's beta on itself is not 1. Affected: monthly betas and standard errors of Table 7 (`res$freq_beta`). | `update_r3.R` first reproduces the stored monthly betas and standard errors with the previous formula (max. abs. difference 8.9e-16), then recomputes them: market beta on itself 1.0000000 (previously 1.0037); other monthly betas change by at most 0.006 (e.g. Q1 0.6218 → 0.6195, Q5 1.4420 → 1.4365); rounded to two decimals, no entry of Table 7 changes (the market's beta of 1.0037 was already shown as 1.00). |
| 3 | Appendix Table 15: calendar-time HAC | `R/functions.R`: `ct_reg` now uses `wls_hac`, which keeps months with zero weight in the calendar with a zero score, so that 12 lags are 12 calendar months; estimate and regression definition unchanged; GFKSW replication method unchanged. Recomputed: `res$T34_ct` (Table 15) and `res$T2_rt_ct` (real-time variant, not shown). | Stored tables reproduced with the previous HAC before the change; point estimates identical after it (tolerance 1e-12); without gaps the new HAC equals sandwich::NeweyWest (difference < 1e-12). t-statistics change by at most about 0.06; all statements and claim checks unchanged. |
| 4 | Covariance-estimation claim | Section 7, introduction, conclusion and abstract: "The daily nonlinear-shrinkage implementation has the lowest realized volatility among the portfolios considered."; Section 7 adds that frequency, window, estimator and stock selection change at the same time, so this is no isolated evidence for the estimator. No new analysis. | No "better covariance estimate" or "better estimation" left in the PDF. |
| 5 | Author list | Gjergji Emini re-inserted after Benjamin Halavac; no other names changed. | Title page shows five authors. |

Only `res$freq_beta`, `res$T34_ct` and `res$T2_rt_ct` were changed (by `update_r3.R`, log in `logs/update_r3.log`;
`res$meta$update_r3` records time and checksums); the bootstraps and size simulations were not rerun. `build.R` and
`R/functions.R` contain the same corrections for future fresh builds; a fresh build was not run in revision 3.
Render: all claim checks passed; 27 pages, main text pp. 1–15 (the conclusion ends at about 60% of p. 15, Figures 1–3
in the main text), references pp. 15–18, appendix pp. 19–27; no text beyond the right margin (max. x = 529.4 pt); pages
1, 14, 15 and 22 checked visually. `lowvol_v2`, `lowvol_v3` and `../lowvol_v4.zip` unchanged.

---

# CHANGELOG — lowvol_v4 (relative to lowvol_v3)

Status words: **implemented** = code or text changed; **executed** = ran in the delivered fresh run on the
university server; **validated** = an explicit check confirms the result (the check is named). Numbers quoted
here come from the delivered run (`results/res.rds`, `res$meta$built`); the paper never reads this file.

`lowvol_v3` and `lowvol_v2` were not modified (md5 of all their files identical before and after the delivered
run, see section 7). The v3 results are kept only as `reference/res_v3_reference_only.rds` for comparisons; neither
`build.R` nor `paper.qmd` reads it. The database was opened read-only and not copied.

Classification used below: **error** = confirmed mistake in code or text; **open** = question the data cannot
settle; **limitation** = legitimate restriction of the design, documented rather than fixed.

## 1. Daily returns and delistings

| # | Finding (classification) | Change | Check |
|---|---|---|---|
| 1 | **Documentation gap, not a numerical error.** v3 rebuilt daily paths from `crsp_daily`, which omits CRSP's delisting-day records, and booked the residual (1 + r_month)/(1 + r_daily) on the last trading day without showing which event it represents. | Daily paths now come from `dsf_v2` (`daily_panel_dsf`, R/functions.R) with the delisting-day records on their own dates. Every stock-month is classified: exact match (≤ 1 bp), delisting record dated in the next month, conflict, unexplained. | Cap-weighted universe portfolio 1968–2024: 378,343 of 378,535 stock-months match; 188 are delisting payments that `dsf_v2` dates in the first days of the next month while `msf_v2` credits them to the held month (booked on the held month's last trading day, proceeds in cash until rebalancing); 4 conflicts; **0 unexplained** (stop in paper.qmd otherwise). Every baseline stock-month compounds to its monthly return (< 1e-10, stop in build.R). **Validated.** |
| 2 | **Open (CRSP inconsistency).** In 4 held stock-months `msf_v2` excludes the delisting-day return that `dsf_v2` records: permno 22357 (Apr 1968, Q5, weight 0.52%), 14867 (Mar 1980), 49365 (Jan 1984: monthly +3.0%, daily −100%), 15923 (Sep 1985). Only 22357 is in a beta quintile; all four are in the cap-weighted (≤ 0.08% weight) and equal-weighted (≤ 0.25%) universe portfolios; the minimum-variance portfolio holds them at numerically zero weight. A fifth case (43546, Dec 1973) is not held by any analysed portfolio. | Baseline follows the monthly file (consistent with the monthly engine). Variants: (i) daily-file return in the conflict months, (ii) v3 reconstruction from `crsp_daily`, (iii) no delisting returns at all. Monthly impact of (i) shown in Appendix Table 9. | Largest change of any portfolio's daily volatility 0.01 pp, Sharpe ratio 0.002, excess kurtosis 0.01, maximum drawdown 0.04 pp, beta 0.0001 (Appendix Table 10). Monthly: mean L–H 0.1562% → 0.1560% per month. **Executed; impact quantified.** |
| 3 | **Error in v3 text.** "Higher annualised daily volatility indicates negative autocorrelation." | Lo–MacKinlay variance ratio of 21-day returns with robust z added (`vr_test`). Text: consistent with mild negative serial correlation, not established. | VR(21) = 0.90 (Q1, z = −0.98) and 0.90 (benchmark, z = −1.02). Claim check. |
| 4 | **Explained difference.** Monthly columns of Table 6 (frequency) vs Table 5 (long-only). | `res$t7t8`: identical monthly returns (max diff 4e-15); Table 6 subtracts the T-bill return compounded from daily rates (differs from the monthly series by up to 11 bp in a month); with the monthly series the Sharpe ratios of Table 5 are reproduced exactly; for the market, Table 6 compounds the daily FF factor (≤ 18 bp per month from the monthly factor). | Claim check (< 1e-8). **Validated.** |
| 5 | Financing and cost basis of the daily paths. | Unchanged and documented: L–H = capital in T-bills + long Q1 − short Q5 (value G1 − G5 + Gf), BN = G1 − 0.25 G5 + 0.25 Gf, legs reset at month ends; Table 6 is gross of costs. | Compounded daily L–H equals the monthly engine (< 1e-10). |
| 6 | Market betas now with Newey-West standard errors (Table 1 text; Section 8 text). | BN realised beta 0.27 (s.e. 0.04), daily 0.32 (0.02): not market-neutral. | Claim check (t > 3). |

## 2. Bootstrap and inference

| # | Finding (classification) | Change | Check |
|---|---|---|---|
| 7 | **Error in v3 interpretation.** The moving-block bootstrap with the state calendar held fixed approximates the distribution of the state averages **under the null** that the states carry no information; v3 reported "normal 95% intervals" from it as if they were confidence intervals. | Estimand and null hypothesis stated for every analysis (Section 3.3): state averages of window returns or window alphas (H0: zero), contrast state 5 − state 1 (H0: no information), calendar-time WLS (different estimand, H0: zero), monthly means (ordinary t). The interval column of Table 2 replaced by the t-statistic of a volatility-preserving bootstrap; contrast reported with p-value from centred draws. | — |
| 8 | Block length. Politis–White (with Patton–Politis–White correction): 2.6 months for the mean of L–H returns, 37.9 for squared returns (BN: 0.3 / 26.6). | Blocks of 12/24/36/60 months reported (Appendix Table 16). | Contrast t (L–H): 1.96 / 1.99 / 2.03 / 2.14; BN 2.04 / 2.02 / 2.20 / 2.18. The 24-month indices are identical to the main bootstrap (stop in build.R). |
| 9 | **Confirmed problem (assumption violated).** Volatility depends on the state: realised L–H volatility over the year after a state-5 month 28.0% p.a., after other states 16.9–18.5%. | Volatility-preserving bootstrap (standardised returns, rescaled with the local volatility of the original calendar position); also for the window alphas (`state_window_stats(vp = TRUE)`, same indices). | Contrast t 1.74 (L–H), 1.77 (BN). One-year means significant only for BN in states 4–5; state-5 CAPM alphas significant only at 48/60 months (L–H) and 24/48/60 months (BN); no four-factor alpha (max |t| 1.74). Point estimates identical to the baseline (stop in build.R). |
| 10 | v3 validated the bootstrap only under iid returns. | Size simulation `sim_tests` (1,000 samples × 299 draws, actual state calendar; Appendix Table 17) under iid normal returns, AR(1)-GARCH(1,1)-t fitted to each strategy, and the fitted AR(1) with the observed volatility path. | iid/GARCH: 24-month bootstrap 3.6–6.1% at nominal 5%; GFKSW-style NW 12.6–17.8%; formation-month ("pairs") bootstrap 10.5–14.7%. Observed volatility: 24-month bootstrap rejects the L–H contrast null in 9.3% (size-adjusted critical value 2.29), BN 5.6% (2.03); volatility-preserving 7.0% (2.10) and 3.8% (1.91). Monte Carlo s.e. of rejection rates ≤ 1.2 pp. **Executed.** |
| 11 | Monte Carlo uncertainty of borderline results. | Delta-method MC s.e. of the bootstrap s.e. reported (Appendix Table 16). | L–H contrast: s.e. 8.40 pp with MC s.e. 0.19 pp → t = 1.99 ± 0.05; p = 0.045 (MC s.e. ≈ 0.007). |
| 12 | **Central result re-classified.** v3: one-year contrast "at the margin of significance (t = 1.99)". | v4: "not statistically established" (abstract, introduction, Section 4.2, conclusion). | Claim checks: t within 1.9–2.2 across blocks; < 1.96 with the volatility-preserving bootstrap; size-adjusted critical value ≥ observed t. |
| 13 | Exploratory test family for beta forecasts incl. D36 vs M36. | Holm adjustment over all ten DM comparisons added (`p_holm_all`). | All five results reported in Section 6 stay significant (largest adjusted p = 0.021). Claim check. |
| 14 | **Error in v3 table note.** Table 4 stated 6 Newey-West lags for Panel A; the rule floor(4(T/100)^(2/9)) gives 5 for T = 540 (6 only for 684). | Lags now computed in the note (5 / 4 / 4). | Same function as `nw_lag()` in R/functions.R. |
| 15 | Calendar gaps, common samples. | GFKSW-style statistic kept only as replication of their (presumed) procedure and labelled as treating episode boundaries as neighbours; Sharpe tests on common samples (unchanged). | — |

## 3. Dates and return definitions

* Re-checked, unchanged (assertions in build.R): fiscal-year rule at May/June/December/January boundaries;
  `mid` = 12·year + month − 1, return month = formation month + 1; 540 formation months Dec 1967–Nov 2012 with
  108 per state; no replication window uses a return after Dec 2012; post-sample 2013–2024 (144 months) and
  post-publication 2016–2024 (108 months); delisting-day records dated in the month after a held month are
  mapped to that held month.
* All components of `results/res.rds` that v4 did not intend to change are identical to v3 (tolerance 1e-10);
  the components that differ are only those listed in sections 1–2 (`T1*` gains a standard-error column,
  `conv_sens` a collateralised column, `data_diag` the conflict list and the daily-file variant, `beta_dm` the
  ten-test Holm column, `recon`, `freq`, `freq_beta`, `daily_check`, `cfg`, `meta`).
* Total, excess and active returns, L–H/BN financing and chained multi-year returns unchanged; the collateralised
  L–H total-return variant is added to `conv_sens` for completeness.

## 4. Claims, figures and text

| v3 statement | v4 |
|---|---|
| One-year contrast "at the margin of significance" | Not statistically established (range of t, volatility-preserving t, size-adjusted critical value) |
| "95% intervals" of state means from the null-model bootstrap | Removed; bootstrap described as test under the null |
| State-5 CAPM alphas significant at every horizon | Only with the baseline bootstrap; volatility-preserving: only at long horizons |
| Higher daily volatility "indicates negative autocorrelation" | Consistent with, not evidence of (variance-ratio test) |
| Low-beta portfolios "not less liquid" | Similar dollar-volume ranks describe the stocks, not tradability or capacity |
| Q1 return "largely spanned" by the four factors | Four-factor alpha insignificant, 95% interval −1.3% to 2.2% p.a.; co-movement, not full explanation |
| Daily betas better forecasts (general) | For these portfolios, this target and this sample |
| BN beta-neutral | Realised beta 0.27 (s.e. 0.04), daily 0.32 (0.02) |
| Break-even costs without context | Values ≥ 50 bp extrapolate beyond the cost scenarios; refer to insignificant differences |

* Figures: three main-text figures — Figure 1 cumulative growth (Section 4.1; Q1, Q5, market; L–H, BN, market,
  T-bills), Figure 2 B/P spread with the state breakpoints, 1967–2024 (Section 4.2), Figure 3 rolling 36/60-month
  Sharpe-ratio differences, CAPM alphas and volatilities (Section 8). Appendix: factor correlations, B/P levels
  since 1929, cost scenarios. The v3 appendix figure of cumulative returns and the duplicated spread panel were
  removed. One colour per series across all figures (Q1 blue, Q5 orange, L–H green, BN purple, market/benchmark
  grey dashed, minimum variance magenta); captions state period, units, portfolios, return and cost basis.
* Tables moved to the appendix: GFKSW Table 2 Panel A (now Appendix Table 8) and the beta-estimator table
  (Appendix Table 19). New appendix tables: daily delisting sensitivity, contrast inference, simulated size.
* New claim checks for every new statement (rendering stops otherwise).

## 5. Requirements, sources and documentation

The v3 CHANGELOG stated that the "Projects and Examination" slides were not available. That statement is
superseded: the slides and the syllabus were provided for v4 and were used for the table below.

| Requirement | Source | Type | Where |
|---|---|---|---|
| Replicate a current empirical-finance paper or test own ideas; focus on the investment strategies | Projects and Examination, slide 2; Syllabus §4 | course requirement | Sections 4–8 |
| Topic: "Elaborate on the low-volatility cycles (Garcia-Feijóo et al. 2015)" | Projects and Examination, slide 6, topic 8 | course topic list | whole paper |
| Written report on methodology and empirical results, ca. 10–15 pages | slide 2 | course requirement | body pp. 1–15 (incl. title page); references and appendix separate |
| Groups of at most 3/4 | slide 2 | course rule | five authors: exception approved for the group (as stated by the group; no document seen here) |
| Project description by 25.09.2026; report by Monday 19.10.2026, 24:00 | slide 2 | deadlines | for the group |
| Presentation ca. 30 min + 5–10 min discussion on 20.10.2026; slides due 19.10.2026 | slide 3 | course requirement | not part of this delivery |
| Submit project description, paper, code and presentation to the lecturer | slide 3 | course requirement | for the group |
| Code in R; code shared with the report | Syllabus §4; slide 4 | course requirement | build.R, R/functions.R, paper.qmd |
| All code managed via GitHub with version control | slide 4 | course requirement | for the group; no GitHub actions here |
| Data: CRSP via Tidy Finance, Fama–French library | slide 4 | course guidance | Section 3.1 |
| Daily vs monthly returns; different rolling windows | De Nard e-mail (as relayed) | agreed extensions | Sections 6 and 8 |
| Realistic implementation constraints; well-diversified long-only comparison; transaction costs | De Nard e-mail (as relayed) | recommendation | Section 7 |
| Evaluation windows/frequency; multiple testing; long-only constraints; min-variance with shrinkage | Lecture 3 slide 7; Lecture 5 slides 14–15, 54, 58; Lecture 4 slides 14, 19–20, 30–33 | course material | Sections 3, 7, 8 |
| Beta-measurement comparison, real-time timing rule, inference simulations | — | own extensions | Sections 3.3, 5, 6 |

* Lecture reference corrected: the minimum-variance slides of Lecture 4 carry the printed numbers 30–33 (PDF
  pages 31–34); v3 cited 31–34. The other cited slide numbers were checked against the PDFs.
* References added: Künsch (1989), Politis and White (2004), Patton, Politis and White (2009), Bollerslev (1986),
  Lo and MacKinlay (1988). De Nard and Zhao is cited as published (IREF 80, 2022, 654–676; the lecture slide
  dates it 2023).
* AI statement: describes the actual use; states that the syllabus and the examination slides contain no rules on
  AI use. **The group must check the university's general rules and confirm the statement before submission.**
* README.md updated for v4.

## 6. Open points and limitations (effect on conclusions)

* Few episodes per valuation state; the size of every test depends on how the volatility of these episodes is
  modelled. Effect: the conditional evidence cannot be made conclusive with these data.
* Four CRSP delisting months with inconsistent monthly and daily files (#2). Effect: immaterial for all reported
  statistics.
* CRSP vintage and delisting conventions differ from GFKSW; effect on the middle states unknown.
* Cost model proportional with scenario parameters; no market-impact or borrowing data. Effect: break-even costs
  are sensitivities.
* Forecast target for betas (12 monthly returns) is noisy; a daily-return target was not implemented.
* The local-volatility estimate of the volatility-preserving bootstrap uses a centred window (future months); it is
  a device for inference under the null, not a trading signal.
* `logs/build.log` contains data.table messages about shallow copies (memory only, no effect on results).

## 7. Delivered run

| Step | Status | Evidence |
|---|---|---|
| Fresh build: `cache/` and `results/res.rds` deleted, then `Rscript build.R` (script `logs/final_run.sh`) | executed, exit code 0 | `logs/build.log`; `logs/build_resources.txt`: started 2026-10-07 06:40:56 UTC, 10 min 45 s, 32 worker cores, peak resident memory 11.1 GB in the largest R process; `res$meta$panel_built_this_run = TRUE` |
| Internal checks in `build.R` (data reconciliation, dates, 540/108 states, no post-2012 return in the replication, daily = monthly engine for every stock-month, same bootstrap draws across checks) | executed and passed | the build stops on any failure; it completed |
| Reproducibility | validated | all 76 result components (everything except the run metadata) bit-identical to the preceding development run with the same code |
| Code state of the results | validated | md5 of `R/functions.R` equals the cache key stored in `res$meta`; `build.R` and `R/functions.R` unchanged since the development run |
| Render (`quarto render paper.qmd --to pdf`, LuaLaTeX) | executed, all claim checks passed | `logs/render.log`: the only warnings are Quarto's notes that the kableExtra tables use LaTeX labels (cross-references resolved by LaTeX; 0 unresolved references in the PDF) |
| Page budget | validated by inspecting the PDF | 26 pages: main text pp. 1–15 (title page with abstract to the end of the conclusion, incl. Figures 1–3 and Tables 1–6), references pp. 16–18, appendix pp. 19–26 |
| `lowvol_v2`, `lowvol_v3` unchanged | validated | md5 of all 60 files identical before the final run and at delivery (`logs/md5_v2v3_before_final.txt`, `logs/md5_v2v3_at_delivery.txt`) |
