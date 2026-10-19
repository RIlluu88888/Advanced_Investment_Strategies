# Diagnostic 2: (a) are the 'G' stock-months (no record after formation) stocks whose security ended on or
# before the formation date? (b) how do msf_v2 / dsf_v2 / stksecurityinfohist record delistings for sample cases?
# Run from lowvol_v3/ after diag/missing_returns.R:  Rscript diag/delist_cases.R
suppressPackageStartupMessages({library(data.table); library(DBI); library(RSQLite)})
source("R/functions.R")
DB <- "/home/shared/data/tidy_finance.sqlite"
sink("diag/delist_cases_out.txt", split = TRUE); options(width = 220)
h <- readRDS("diag/missing_returns.rds")
con <- dbConnect(SQLite(), DB, flags = SQLITE_RO)
G <- h[startsWith(cat, "G")]
ids <- paste(unique(G$permno), collapse = ",")
msf <- setDT(dbGetQuery(con, paste0("SELECT permno, mthcaldt, mthret, mthretx, mthdelflg, mthretflg, mthprcdt, tradingstatusflg FROM msf_v2 WHERE permno IN (", ids, ")")))
msf[, mid := mid_of(to_idate(mthcaldt))]
sih <- setDT(dbGetQuery(con, paste0("SELECT permno, secinfostartdt, secinfoenddt, securitybegdt, securityenddt, tradingstatusflg, primaryexch, ",
                                    "delactiontype, delstatustype, delreasontype, delpaymenttype FROM stksecurityinfohist WHERE permno IN (", ids, ")")))
for (v in c("secinfostartdt", "secinfoenddt", "securitybegdt", "securityenddt")) set(sih, j = v, value = to_idate(sih[[v]]))
# formation-month record of each G stock-month
G[msf, ":="(f_ret = i.mthret, f_del = i.mthdelflg, f_retflg = i.mthretflg, f_prcdt = to_idate(i.mthprcdt)), on = .(permno, mid)]
send <- sih[, .(securityenddt = max(securityenddt)), by = permno]
G[send, securityenddt := i.securityenddt, on = "permno"]
G[, formation_end := as.IDate(date_of_mid(mid + 1L)) - 1L]          # last calendar day of the formation month
cat("== (a) G stock-months: delisting flag of the FORMATION-month record and security end date ==\n")
print(G[, .N, keyby = .(f_del, f_retflg)])
cat("\nsecurity end date relative to the formation month:\n")
print(G[, .(n = .N, end_in_or_before_formation_month = sum(securityenddt <= formation_end),
            end_after = sum(securityenddt > formation_end), missing = sum(is.na(securityenddt)))])
cat("\nlast trading date (mthprcdt) of the formation-month record vs. security end date:\n")
print(G[, .(same_day = sum(f_prcdt == securityenddt, na.rm = TRUE), n = .N)])
# delisting reasons from the record that carries them (any record of the permno with a non-N/A action)
dr <- sih[!is.na(delactiontype) & delactiontype != "N/A", .SD[.N], by = permno]
G[dr, ":="(delact = i.delactiontype, delreason = i.delreasontype, delpay = i.delpaymenttype), on = "permno"]
print(G[, .N, keyby = .(delact, delreason)])
cat("\nexamples:\n")
print(G[order(-w)][1:8, .(mid = format(date_of_mid(mid), "%Y-%m"), q, permno, w = round(w, 4), f_ret, f_del, f_retflg, f_prcdt, securityenddt, delact, delreason)])
cat("\nsecurity-info records of the largest G case:\n")
print(sih[permno == G[order(-w)][1, permno]][order(secinfostartdt)])

cat("\n== (b) daily records around delistings (dsf_v2 vs crsp_daily) ==\n")
cases <- data.table(permno = c(47079L, 80599L, 92257L, 30971L, 55685L), y = c(1998L, 2008L, 2023L, 1975L, 1981L), m = c(10L, 9L, 11L, 6L, 4L))
idc <- paste(cases$permno, collapse = ",")
dsf <- setDT(dbGetQuery(con, paste0("SELECT * FROM dsf_v2 WHERE permno IN (", idc, ")")))
dly <- setDT(dbGetQuery(con, paste0("SELECT permno, date, ret FROM crsp_daily WHERE permno IN (", idc, ")")))
dbDisconnect(con)
cat("dsf_v2 columns:", paste(names(dsf), collapse = ", "), "\n")
dsf[, date := to_idate(dlycaldt)]; dly[, date := to_idate(date)]
for (k in seq_len(nrow(cases))) {
  m <- ym(cases$y[k], cases$m[k])
  x <- dsf[permno == cases$permno[k] & mid_of(date) %in% (m + 0:1), intersect(c("permno", "date", "dlyret", "dlyretx", "dlydelflg", "dlyretmissflg", "dlyprc", "dlyprcflg", "tradingstatusflg"), names(dsf)), with = FALSE]
  x[dly[permno == cases$permno[k]], crsp_daily_ret := i.ret, on = .(permno, date)]
  cat("\npermno", cases$permno[k], format(date_of_mid(m), "%Y-%m"), ": last 4 dsf_v2 rows of the month and any later rows\n")
  print(tail(x[order(date)], 6))
  cat("compounded dsf_v2 dlyret in month:", round(x[mid_of(date) == m, prod(1 + dlyret[is.finite(dlyret)]) - 1], 4),
      "| compounded crsp_daily:", round(x[mid_of(date) == m, prod(1 + crsp_daily_ret[is.finite(crsp_daily_ret)]) - 1], 4), "\n")
}
sink()
