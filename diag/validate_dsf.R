# Diagnostic 3: validate the recovered delisting-month returns against dsf_v2 for ALL recovered stock-months,
# describe the non-delisting months that exist only in msf_v2, and check the volume fields.
# Run from lowvol_v3/ after diag/missing_returns.R:  Rscript diag/validate_dsf.R
suppressPackageStartupMessages({library(data.table); library(DBI); library(RSQLite)})
source("R/functions.R")
DB <- "/home/shared/data/tidy_finance.sqlite"
sink("diag/validate_dsf_out.txt", split = TRUE); options(width = 200)
h <- readRDS("diag/missing_returns.rds")
D <- h[startsWith(cat, "D")]
con <- dbConnect(SQLite(), DB, flags = SQLITE_RO)
ids <- paste(unique(D$permno), collapse = ",")
dsf <- setDT(dbGetQuery(con, paste0("SELECT permno, dlycaldt, dlyret, dlydelflg, tradingstatusflg FROM dsf_v2 WHERE permno IN (", ids, ")")))
dsf[, rmid := mid_of(to_idate(dlycaldt))]
agg <- dsf[, .(r_dsf = prod(1 + dlyret[is.finite(dlyret)]) - 1, n_del = sum(dlydelflg == "Y"), n_days = .N), by = .(permno, rmid)]
D[agg, ":="(r_dsf = i.r_dsf, n_del = i.n_del, n_dsf = i.n_days), on = .(permno, rmid)]
cat("== recovered stock-months (only in msf_v2): msf_v2 mthret vs compounded dsf_v2 dlyret ==\n")
print(D[, .(n = .N, with_dsf = sum(is.finite(r_dsf)), max_abs_diff = max(abs(msf_ret - r_dsf), na.rm = TRUE),
            share_diff_gt_1bp = mean(abs(msf_ret - r_dsf) > 1e-4, na.rm = TRUE), with_delisting_day = sum(n_del > 0, na.rm = TRUE))])
cat("cases with |mthret - compounded dlyret| > 1 bp:\n"); print(D[abs(msf_ret - r_dsf) > 1e-4, .(rmid = format(date_of_mid(rmid), "%Y-%m"), q, permno, w = round(w, 4), msf_ret, r_dsf, rd, ndays, n_dsf, n_del, msf_del, msf_retflg)])
print(D[, .(n = .N, with_delisting_day = sum(n_del > 0, na.rm = TRUE)), keyby = .(msf_del, msf_retflg, msf_trd)])
cat("\n== non-delisting months only in msf_v2 (mthdelflg N) ==\n")
x <- D[msf_del == "N"]
m2 <- setDT(dbGetQuery(con, paste0("SELECT permno, mthcaldt, primaryexch, tradingstatusflg, sharetype, securitytype, conditionaltype, mthretflg FROM msf_v2 WHERE permno IN (",
                                   paste(unique(x$permno), collapse = ","), ")")))
m2[, rmid := mid_of(to_idate(mthcaldt))]
x[m2, ":="(pexch = i.primaryexch, trd = i.tradingstatusflg, shrt = i.sharetype, cond = i.conditionaltype), on = .(permno, rmid)]
print(x[, .N, keyby = .(msf_retflg, pexch, trd, cond)])
print(x[order(-w)][1:10, .(rmid = format(date_of_mid(rmid), "%Y-%m"), q, permno, w = round(w, 4), msf_ret, rd, ndays, msf_retflg, pexch, trd, cond)])
cat("\n== volume fields: mthprcvol vs mthvol * |mthprc| (sample of 2010-2012 rows) ==\n")
v <- setDT(dbGetQuery(con, "SELECT permno, mthcaldt, mthvol, mthprcvol, mthprc, shrout, primaryexch FROM msf_v2 WHERE mthcaldt BETWEEN 14610 AND 15705 AND mthvol > 0 LIMIT 200000"))
v[, ratio := mthprcvol / (mthvol * abs(mthprc))]
print(v[, as.list(quantile(ratio, c(.01, .1, .5, .9, .99), na.rm = TRUE))])
print(v[, .(med_turnover = median(mthvol / (shrout * 1000), na.rm = TRUE), med_turnover_raw = median(mthvol / shrout, na.rm = TRUE)), by = primaryexch])
cat("volume coverage by decade (share of msf_v2 rows with mthvol > 0):\n")
cv <- setDT(dbGetQuery(con, "SELECT CAST(yyyymm / 1000 AS INTEGER) AS dec3, AVG(CASE WHEN mthvol > 0 THEN 1.0 ELSE 0.0 END) AS cov FROM msf_v2 GROUP BY dec3"))
print(cv)
dbDisconnect(con)
sink()
