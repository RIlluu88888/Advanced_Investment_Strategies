# Diagnostic (v4, exploration): data flows between msf_v2, dsf_v2 and crsp_daily for the stocks held by the
# beta-sorted portfolios, the four delisting months in which CRSP's monthly and daily files disagree, and the
# timing of delisting-day records. Read-only. Uses the holdings of the v3 panel cache (read-only) because the
# v4 cache does not exist yet; the numbers used in the paper come from build.R.
# Run from lowvol_v4/:  Rscript diag/daily_delist_explore.R > logs/diag_daily_delist_explore.log 2>&1
suppressPackageStartupMessages({library(data.table); library(DBI); library(RSQLite)})
source("R/functions.R")
options(width = 220, datatable.print.nrows = 200, datatable.print.class = FALSE)
DB <- "/home/shared/data/tidy_finance.sqlite"
con <- dbConnect(SQLite(), DB, flags = SQLITE_RO)
cm <- readRDS("../lowvol_v3/cache/panel.rds")
held <- cm[large == TRUE & mid >= ym(1967L, 12L) & mid <= ym(2024L, 11L),
           .(permno, fmid = mid, rmid = mid + 1L, mktcap, ret_next, src_next)]
cat("held stock-months (formation 1967-12 .. 2024-11):", nrow(held), "\n")
hq <- form_holdings(cm, "beta_M60"); hq[, w := mktcap / sum(mktcap), by = .(mid, q)]

# ---- 1. the four conflict cases ------------------------------------------------------------------
cases <- data.table(permno = c(22357L, 14867L, 49365L, 15923L), rmid = c(ym(1968L, 4L), ym(1980L, 3L), ym(1984L, 1L), ym(1985L, 9L)))
ids <- paste(cases$permno, collapse = ",")
m <- setDT(dbGetQuery(con, paste0("SELECT permno, yyyymm, mthcaldt, mthprc, mthprcdt, mthdelflg, mthret, mthretx, mthretflg,
  mthprevprc, mthprevdt, tradingstatusflg, mthcompflg FROM msf_v2 WHERE permno IN (", ids, ")")))
m[, rmid := mid_of(to_idate(mthcaldt))]
d <- setDT(dbGetQuery(con, paste0("SELECT permno, dlycaldt, dlydelflg, dlyprc, dlyprcflg, dlyprevprc, dlyprevdt, dlyret, dlyretx,
  dlyretmissflg, dlyretdurflg, dlydistretflg, tradingstatusflg, dlyorddivamt, dlynonorddivamt, dlyvol FROM dsf_v2 WHERE permno IN (", ids, ")")))
d[, ":="(date = to_idate(dlycaldt), rmid = mid_of(to_idate(dlycaldt)))]
cd <- setDT(dbGetQuery(con, paste0("SELECT permno, date, ret FROM crsp_daily WHERE permno IN (", ids, ")")))
cd[, ":="(date = to_idate(date), rmid = mid_of(to_idate(date)))]
si <- setDT(dbGetQuery(con, paste0("SELECT permno, secinfostartdt, secinfoenddt, securitybegdt, securityenddt, tradingstatusflg,
  securityactiveflg, delactiontype, delstatustype, delreasontype, delpaymenttype, securitynm FROM stksecurityinfohist WHERE permno IN (", ids, ")")))
for (v in c("secinfostartdt", "secinfoenddt", "securitybegdt", "securityenddt")) set(si, j = v, value = as.character(to_idate(si[[v]])))
for (i in seq_len(nrow(cases))) {
  p <- cases$permno[i]; r <- cases$rmid[i]
  cat("\n================ permno", p, "| return month", format(date_of_mid(r), "%Y-%m"), "================\n")
  h <- hq[permno == p & mid == r - 1L]
  if (nrow(h)) cat("held at formation", format(date_of_mid(h$mid), "%Y-%m"), "| quintile", h$q, "| weight in quintile", round(h$w, 5),
                   "| src_next:", h$src_next, "| ret_next (used by the monthly engine):", h$ret_next, "\n")
  x <- m[permno == p & rmid %in% (r - 1L):(r + 1L)]
  x[, ":="(mthcaldt = as.character(to_idate(mthcaldt)), mthprcdt = as.character(to_idate(mthprcdt)), mthprevdt = as.character(to_idate(mthprevdt)))]
  cat("-- msf_v2\n"); print(x[, !"rmid"])
  y <- d[permno == p & rmid == r]; setorder(y, date)
  y[, ":="(dlycaldt = as.character(date), dlyprevdt = as.character(to_idate(dlyprevdt)))]
  cat("-- dsf_v2 (return month)\n"); print(y[, !c("date", "rmid")])
  z <- cd[permno == p & rmid == r]; setorder(z, date)
  cat("-- crsp_daily (return month): n =", nrow(z), "| last date", as.character(max(z$date)), "| compounded", prod(1 + z$ret) - 1, "\n")
  cat("-- compounded dsf_v2 all records:", y[is.finite(dlyret), prod(1 + dlyret) - 1],
      "| without delisting-day record:", y[is.finite(dlyret) & dlydelflg != "Y", prod(1 + dlyret) - 1], "\n")
  cat("-- stksecurityinfohist (last rows)\n"); print(tail(si[permno == p][order(secinfostartdt)], 3))
}

# ---- 2. all held stock-months: compounded dsf_v2 vs msf_v2 vs crsp_daily -------------------------
pl <- unique(held$permno)
cat("\nloading dsf_v2 and crsp_daily for", length(pl), "permnos\n")
tmp <- data.table(permno = pl); dbWriteTable(con, "tmp_ids", tmp, temporary = TRUE)
ds <- setDT(dbGetQuery(con, "SELECT d.permno, d.dlycaldt, d.dlyret, d.dlydelflg, d.tradingstatusflg FROM dsf_v2 d JOIN tmp_ids t ON d.permno = t.permno"))
ds[, ":="(date = to_idate(dlycaldt), rmid = mid_of(to_idate(dlycaldt)))][, dlycaldt := NULL]
ds <- ds[held[, .(permno, rmid)], on = .(permno, rmid), nomatch = 0L]
cr <- setDT(dbGetQuery(con, "SELECT c.permno, c.date, c.ret FROM crsp_daily c JOIN tmp_ids t ON c.permno = t.permno"))
cr[, ":="(date = to_idate(date), rmid = mid_of(to_idate(date)))]
cr <- cr[held[, .(permno, rmid)], on = .(permno, rmid), nomatch = 0L]
a1 <- ds[, .(r_dsf = prod(1 + dlyret[is.finite(dlyret)]) - 1, n_dsf = .N, n_dsf_ret = sum(is.finite(dlyret)),
             n_del = sum(dlydelflg == "Y"), del_date = if (any(dlydelflg == "Y")) max(date[dlydelflg == "Y"]) else as.IDate(NA),
             last_trade = if (any(dlydelflg != "Y")) max(date[dlydelflg != "Y"]) else as.IDate(NA),
             r_del = if (any(dlydelflg == "Y")) prod(1 + dlyret[dlydelflg == "Y" & is.finite(dlyret)]) - 1 else NA_real_,
             stat = paste(sort(unique(tradingstatusflg)), collapse = "")), by = .(permno, rmid)]
a2 <- cr[, .(r_cd = prod(1 + ret[is.finite(ret)]) - 1, n_cd = .N, last_cd = max(date)), by = .(permno, rmid)]
H <- a2[a1[held, on = .(permno, rmid)], on = .(permno, rmid)]
H[, ":="(dev_dsf = abs((1 + ret_next) / (1 + r_dsf) - 1), dev_cd = abs((1 + ret_next) / (1 + r_cd) - 1))]
cat("\n== held stock-months by source of the monthly return ==\n")
print(H[, .(n = .N, with_dsf = sum(is.finite(r_dsf)), with_cd = sum(is.finite(r_cd)),
            dsf_gt1bp = sum(dev_dsf > 1e-4, na.rm = TRUE), dsf_gt10bp = sum(dev_dsf > 1e-3, na.rm = TRUE),
            cd_gt1bp = sum(dev_cd > 1e-4, na.rm = TRUE), cd_gt10bp = sum(dev_cd > 1e-3, na.rm = TRUE),
            with_del_record = sum(n_del > 0, na.rm = TRUE)), keyby = src_next])
cat("\n== delisting-day records: date relative to the last trading record (calendar days) ==\n")
dd <- H[n_del > 0]
dd[, gap := as.integer(del_date - last_trade)]
print(dd[, .(n = .N, gap_min = min(gap, na.rm = TRUE), gap_med = median(gap, na.rm = TRUE), gap_max = max(gap, na.rm = TRUE),
             same_month = mean(mid_of(del_date) == rmid), del_ret_na = sum(!is.finite(r_del)), del_minus100 = sum(r_del <= -0.999, na.rm = TRUE))])
print(dd[, .N, keyby = .(gap = pmin(gap, 10L))])
cat("\n== stock-months where compounded dsf_v2 differs from msf_v2 by more than 10 bp ==\n")
print(H[dev_dsf > 1e-3, .(rmid = format(date_of_mid(rmid), "%Y-%m"), permno, src_next, ret_next, r_dsf, r_cd, n_dsf, n_dsf_ret, n_cd, n_del, stat)][order(rmid)])
cat("\n== stock-months without any dsf_v2 return record ==\n")
print(H[!is.finite(r_dsf) | n_dsf_ret == 0, .N, keyby = src_next])
cat("\n== stock-months where crsp_daily (used by v3) differs from msf_v2 by more than 10 bp: by source ==\n")
print(H[dev_cd > 1e-3, .N, keyby = .(src_next, has_del = n_del > 0)])
dbDisconnect(con)
saveRDS(H, "diag/daily_delist_explore.rds")
