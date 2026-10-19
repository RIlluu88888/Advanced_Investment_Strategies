# Diagnostic: why do held large-cap stocks lack a next-month return in crsp_monthly?
# Run from lowvol_v3/:  Rscript diag/missing_returns.R
# Input (read-only): the university database and, for the holdings only, the v2 panel cache.
# Output: diag/missing_returns_out.txt, diag/missing_returns.rds. Changes no result file.
suppressPackageStartupMessages({library(data.table); library(DBI); library(RSQLite)})
source("R/functions.R")
DB <- "/home/shared/data/tidy_finance.sqlite"
sink("diag/missing_returns_out.txt", split = TRUE); options(width = 200)
cm <- readRDS("../lowvol_v2/cache/panel.rds")          # v2 panel, used here only to define holdings
h <- form_holdings(cm, "beta_M60")
h <- h[mid >= ym(1967L, 12L) & mid <= ym(2024L, 11L)]
h[, w := mktcap / sum(mktcap), by = .(mid, q)]
h[, wvw := mktcap / sum(mktcap), by = mid]               # weight in the cap-weighted universe
h[, rmid := mid + 1L]
con <- dbConnect(SQLite(), DB, flags = SQLITE_RO)
msf <- setDT(dbGetQuery(con, paste("SELECT permno, mthcaldt, mthret, mthretx, mthdelflg, mthretflg, primaryexch,",
                                   "tradingstatusflg, conditionaltype, sharetype, securitytype FROM msf_v2")))
cmr <- setDT(dbGetQuery(con, "SELECT permno, date, ret, exchange FROM crsp_monthly"))
sih <- setDT(dbGetQuery(con, paste("SELECT permno, secinfostartdt, secinfoenddt, securityenddt, delactiontype, delstatustype,",
                                   "delreasontype, delpaymenttype FROM stksecurityinfohist")))
dly <- setDT(dbGetQuery(con, "SELECT permno, date, ret FROM crsp_daily"))
dbDisconnect(con)
msf[, rmid := mid_of(to_idate(mthcaldt))]; cmr[, rmid := mid_of(to_idate(date))]
stopifnot(!anyDuplicated(msf, by = c("permno", "rmid")), !anyDuplicated(cmr, by = c("permno", "rmid")))
cat("== 0. sanity: crsp_monthly vs msf_v2 on common rows ==\n")
cc <- merge(cmr[, .(permno, rmid, ret)], msf[, .(permno, rmid, mthret)], by = c("permno", "rmid"))
cat("common rows:", nrow(cc), " max |ret - mthret|:", cc[is.finite(ret) & is.finite(mthret), max(abs(ret - mthret))],
    " rows only in crsp_monthly:", nrow(cmr) - nrow(cc), "\n")
h[cmr, ":="(cm_ret = i.ret, cm_exch = i.exchange, cm_row = TRUE), on = .(permno, rmid)]
h[msf, ":="(msf_ret = i.mthret, msf_del = i.mthdelflg, msf_retflg = i.mthretflg, msf_trd = i.tradingstatusflg,
            msf_row = TRUE), on = .(permno, rmid)]
h[is.na(cm_row), cm_row := FALSE][is.na(msf_row), msf_row := FALSE]
dd <- dly[, .(ndays = .N, rd = prod(1 + ret[is.finite(ret)]) - 1), by = .(permno, rmid = mid_of(to_idate(date)))]
h[dd, ":="(ndays = i.ndays, rd = i.rd), on = .(permno, rmid)]
h[is.na(ndays), ndays := 0L]
h[, cat := fcase(
  cm_row & cm_exch %in% c("NYSE", "AMEX", "NASDAQ") & is.finite(cm_ret), "A: crsp_monthly return (listed exchange)",
  cm_row & !(cm_exch %in% c("NYSE", "AMEX", "NASDAQ")) & is.finite(cm_ret), "B: crsp_monthly return, exchange 'Other' (dropped by v2 filter)",
  cm_row & !is.finite(cm_ret), "C: crsp_monthly row, missing return",
  !cm_row & msf_row & is.finite(msf_ret), "D: only msf_v2 has the month, return available",
  !cm_row & msf_row & !is.finite(msf_ret), "E: only msf_v2 row, return missing",
  !cm_row & !msf_row & ndays > 0, "F: no monthly record, daily records exist",
  default = "G: no monthly and no daily record")]
cat("\n== 1. held stock-months by availability of the next-month return (formation 1967-12 to 2024-11) ==\n")
print(h[, .(stock_months = .N, avg_q_weight_pct = 100 * sum(w) / uniqueN(paste(mid, q)),
            avg_vw_weight_pct = 100 * sum(wvw) / uniqueN(mid)), keyby = cat])
cat("\nby quintile (average weight per month, %):\n")
print(dcast(h[!startsWith(cat, "A")], cat ~ q, value.var = "w", fun.aggregate = function(x) 100 * sum(x) / uniqueN(h$mid)))
cat("\nD: delisting flag / return flag in msf_v2:\n"); print(h[startsWith(cat, "D"), .N, keyby = .(msf_del, msf_retflg)])
cat("\nD: distribution of the recovered returns:\n"); print(h[startsWith(cat, "D"), as.list(quantile(msf_ret, c(0, .01, .05, .25, .5, .75, .95, .99, 1)))])
cat("\nD: largest weighted effects in Q1 and Q5 (weight x return, bp):\n")
print(h[startsWith(cat, "D") & q %in% c(1, 5)][order(-abs(w * msf_ret))][1:15, .(mid = format(date_of_mid(rmid), "%Y-%m"), q, permno, w = round(w, 4),
      msf_ret = round(msf_ret, 3), daily_partial = round(rd, 3), ndays, msf_del)])
last_sih <- sih[order(permno, -secinfostartdt)][, .SD[1], by = permno]
h[last_sih, ":="(delact = i.delactiontype, delreason = i.delreasontype, delpay = i.delpaymenttype), on = "permno"]
cat("\nE/F/G and C: details\n")
print(h[substr(cat, 1, 1) %in% c("C", "E", "F", "G"), .(n = .N, w_q = round(100 * sum(w) / uniqueN(h$mid), 4)),
        keyby = .(cat, delact, delreason)])
print(h[substr(cat, 1, 1) %in% c("C", "E", "F", "G")][order(-w)][1:20, .(mid = format(date_of_mid(rmid), "%Y-%m"), q, permno, w = round(w, 4),
      cat = substr(cat, 1, 1), ndays, rd = round(rd, 3), delact, delreason, delpay)])
lastm <- msf[, .(last_rmid = max(rmid)), by = permno]
h[lastm, last_rmid := i.last_rmid, on = "permno"]
cat("\nG: stock has msf_v2 records after the missing month (gap) vs. none (exit):\n")
print(h[startsWith(cat, "G"), .(n = .N, gap = sum(last_rmid > rmid, na.rm = TRUE), exit_or_unknown = sum(!(last_rmid > rmid) | is.na(last_rmid)))])
saveRDS(h[!startsWith(cat, "A")], "diag/missing_returns.rds")
sink()
