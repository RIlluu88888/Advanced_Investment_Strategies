# Low-Volatility Cycles Revisited



## Contents

| Path | Purpose |
|---|---|
| `build.R` | Computes every number, table and series used in the paper; writes `results/res.rds` |
| `update_r3.R` | Revision 3: applies the corrected Table 7 market betas and calendar-time standard errors to `results/res.rds` without rerunning the bootstraps and simulations (log: `logs/update_r3.log`) |
| `R/functions.R` | Function library (data loading, betas, portfolios, turnover, linear and nonlinear shrinkage, daily paths, inference incl. reduced-bias predictive regression, simulations) |
| `paper.qmd` | Quarto source of the paper; reads only `results/res.rds` |
| `paper.pdf` | Rendered on the university server from `paper.qmd` and `results/res.rds` (see `logs/render.log`) |
| `refs.bib` | References |
| `data/DFF_BE_With_Nonindust.txt` | Davis–Fama–French historical book equity (Kenneth French data library) |
| `results/res.rds` | All results of the delivered `build.R` run (see `res$meta`) |
| `diag/` | Data diagnostics (delisting months, daily vs monthly files) with their scripts |
| `logs/` | `build.log`, `build_resources.txt`, `render.log` and `final_run.sh` (the script that ran the fresh build and logged its resources) belong to the delivered run of revision 2; `md5_v2v3_*.txt` document that v2 and v3 were not changed; `dev_run*.log`, `render_dev*.log` and `diag_*.log` are development and diagnostic runs |
| `reference/res_v3_reference_only.rds` | v3 results, kept **only** for comparison; never read by `build.R` or `paper.qmd` |
| `CHANGELOG.md` | Revisions 2 and 3 (top) and the findings, changes, checks, requirement match and open points of v4 |

## Prerequisites

* R ≥ 4.4 with `data.table`, `DBI`, `RSQLite`, `sandwich`, `quadprog` (build) and `knitr`,
  `kableExtra`, `ggplot2`, `scales`, `patchwork` (paper). Versions used: Appendix C of the paper or
  `res$meta$packages`.
* The university SQLite database `/home/shared/data/tidy_finance.sqlite` (tables `crsp_monthly`, `crsp_daily`,
  `msf_v2`, `dsf_v2`, `compustat_annual` and the Fama–French factor tables). It is opened
  **read-only** and not copied. Another path can be passed as the first argument.
* Quarto ≥ 1.4 with a LaTeX engine (rendered with Quarto 1.9.38 and LuaHBTeX, TeX Live 2022).
* Resources of the delivered run: see `logs/build_resources.txt` (wall time 12 min 54 s, peak memory 7.4 GB in the
  largest R process). The build uses up to 32 worker cores; results do not depend on the number of cores.

## How to run

```bash
cd lowvol_v4
rm -rf cache results/res.rds                      # force a fresh build
Rscript build.R > logs/build.log 2>&1             # about 13 minutes on the university server (32 cores)
quarto render paper.qmd --to pdf > logs/render.log 2>&1
```

The size simulations of the inference procedures run inside `build.R` (`res$sim_inference`,
`res$sim_contrast`) and take most of the run time.

`build.R` stores intermediate data in `cache/` together with a key (database path, size and modification time;
md5 of `R/functions.R`, of the data section of `build.R` and of the book-equity file). The cache is reused only
when the key matches; otherwise it is rebuilt. `res$meta$panel_built_this_run` records whether the delivered
results were built from scratch.

