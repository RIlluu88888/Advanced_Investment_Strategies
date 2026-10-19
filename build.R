# =============================================================================
# Low-Volatility Cycles (Garcia-Feijoo, Kochard, Sullivan & Wang, FAJ 2015)
# Replication + extensions, version 4 (revision 2). Computes every number used in paper.qmd.
# Usage (from this folder):  Rscript build.R [db_path] [cores]
# Output: results/res.rds (all tables and series), cache/*.rds (intermediate), logs via stdout.
# The database is opened read-only. The panel cache is reused only if its key (database path,
# size and modification time; md5 of R/functions.R, of this file's data section and of the
# book-equity file) matches; otherwise it is rebuilt. A fresh run: rm -rf cache && Rscript build.R
# Random numbers are used only for the block bootstraps and the size simulations (seeds derived from cfg$seed,
# L'Ecuyer-CMRG; bootstrap indices are drawn in the main process and every simulated sample has its own seed,
# so results do not depend on the number of cores).
# =============================================================================
args  <- commandArgs(trailingOnly = TRUE)
DB    <- if (length(args) >= 1) args[1] else "/home/shared/data/tidy_finance.sqlite"
CORES <- if (length(args) >= 2) as.integer(args[2]) else max(1L, min(32L, parallel::detectCores() - 2L))
source("R/functions.R")
set.seed(20261006)
dir.create("results", showWarnings = FALSE); dir.create("cache", showWarnings = FALSE)
t_start <- Sys.time()
logmsg("database:", DB, "| cores:", CORES)
stopifnot(file.exists(DB))
DFF <- "data/DFF_BE_With_Nonindust.txt"
db_info <- file.info(DB)
self <- readLines("build.R")
data_section <- self[seq(grep("^# ---- 1. data", self)[1], grep("^logmsg\\(\"panel ready\"\\)", self)[1])]
tf <- tempfile(); writeLines(data_section, tf)
cache_key <- list(db = normalizePath(DB), db_size = db_info$size, db_mtime = format(db_info$mtime, "%Y-%m-%d %H:%M:%S"),
                  functions_md5 = unname(tools::md5sum("R/functions.R")), dff_md5 = unname(tools::md5sum(DFF)),
                  data_code_md5 = unname(tools::md5sum(tf)))
unlink(tf)

cfg <- list(
  rep  = c(ym(1968L, 1L), ym(2012L, 12L)),   # GFKSW main sample (return months)
  long = c(ym(1930L, 1L), ym(2012L, 12L)),   # GFKSW longer evidence
  post = c(ym(2013L, 1L), ym(2024L, 12L)),   # after GFKSW's sample end (post-sample extension)
  pub  = c(ym(2016L, 1L), ym(2024L, 12L)),   # after publication (FAJ, May/June 2015)
  full = c(ym(1968L, 1L), ym(2024L, 12L)),   # full modern sample
  horizons = c(1L, 12L, 24L, 36L, 48L, 60L),
  costs = c(0, 5, 10, 20, 50) / 1e4,           # one-way cost per dollar traded
  base_cost = 10 / 1e4, borrow = 0.0025,       # 25 bp p.a. lending fee on short notional
  bn_short = 0.25, minvar_window = 60L, minvar_wmax = 0.03, cap = 0.05,
  buffer = c(0.2, 0.3), bp_trim = 0.005, ct_lag = 12L, rt_min = 120L, timing_states = 4:5,
  boot_B = 999L, boot_block = 24L, seed = 20261006L, sim_n = 1000L, sim_B = 299L,
  specs = c("beta_M60", "beta_M36", "beta_M120", "beta_D3", "beta_D12", "beta_D36", "beta_DD12"))
in_period <- function(m, p) m >= p[1] & m <= p[2]

# ---- 1. data ---------------------------------------------------------------------------
old_key <- if (file.exists("cache/key.rds")) readRDS("cache/key.rds") else NULL
panel_built <- !(file.exists("cache/panel.rds") && identical(old_key, cache_key))
if (!panel_built) logmsg("cache key matches: reusing cache/panel.rds built", format(file.info("cache/panel.rds")$mtime))
if (panel_built) {
  if (dir.exists("cache")) unlink(list.files("cache", full.names = TRUE))
  con <- dbConnect(SQLite(), DB, flags = SQLITE_RO)
  fac <- load_factors(con); fd <- load_daily_factors(con)
  logmsg("monthly CRSP")
  cm <- load_monthly(con)
  msf <- load_msf(con)
  # crsp_monthly and msf_v2 must agree wherever both have a return (no double counting of delisting returns)
  chk_cm <- msf[cm[, .(permno, mid, ret)], on = .(permno, mid), nomatch = 0L]
  data_check <- list(common_rows = nrow(chk_cm), cm_rows = nrow(cm),
                     max_abs_diff = chk_cm[is.finite(ret) & is.finite(mthret), max(abs(ret - mthret))],
                     cm_rows_without_msf = nrow(cm) - nrow(chk_cm))
  stopifnot(data_check$max_abs_diff < 1e-12, data_check$cm_rows_without_msf == 0L)
  cm <- add_monthly_betas(cm, fac)
  be <- load_book_equity(con, DFF)
  cm <- add_bp(cm, be)
  cm <- add_universe(cm, msf)
  logmsg("monthly panel:", nrow(cm), "rows; large-cap stock-months:", cm[large == TRUE, .N],
         "; excluded (delisted in formation month):", cm[large_all & delisted_t, .N], "; with B/P:", cm[univ == TRUE, .N])
  rets_m <- msf[is.finite(mthret), .(permno, mid, ret = mthret)]; setkey(rets_m, permno, mid)
  rm(msf, chk_cm); invisible(gc())
  logmsg("daily CRSP (one pass)")
  d <- setDT(dbGetQuery(con, "SELECT permno, date, ret, ret_excess FROM crsp_daily"))
  d[, date := to_idate(date)]
  d[, mid := year(date) * 12L + month(date) - 1L]
  held <- cm[large_all == TRUE, .(permno, mid = mid + 1L)]
  dret <- d[held, on = .(permno, mid), nomatch = 0L, .(permno, date, mid, ret)]
  setkey(dret, permno, date)
  saveRDS(dret, "cache/daily_held_cd.rds"); rm(dret)
  # daily returns for the covariance matrix of the nonlinear-shrinkage minimum-variance portfolio: every
  # stock-month within the 12 months up to a formation month (from December 1967) in which the stock is
  # in the large-cap universe
  uk <- cm[large == TRUE & mid >= ym(1967L, 12L), .(permno, t = mid)]
  uk <- unique(rbindlist(lapply(0:11, function(k) uk[, .(permno, mid = t - k)])))
  du <- d[uk, on = .(permno, mid), nomatch = 0L, .(permno, date, mid, ret)]
  setkey(du, mid)
  saveRDS(du, "cache/daily_univ.rds"); rm(uk, du); invisible(gc())
  # CRSP's daily file (dsf_v2) for the held stocks, including the delisting-day records (dlydelflg 'Y') that
  # crsp_daily omits. In-month records feed the daily portfolio returns (section 7); delisting-day records dated
  # in the month after a held month are kept separately (CRSP's monthly file credits them to the held month).
  logmsg("daily CRSP file dsf_v2 for the held stocks")
  dbWriteTable(con, "tmp_held_permno", data.table(permno = unique(held$permno)), temporary = TRUE)
  ds <- setDT(dbGetQuery(con, paste("SELECT d.permno, d.dlycaldt, d.dlyret, d.dlydelflg FROM dsf_v2 d",
                                    "JOIN tmp_held_permno t ON d.permno = t.permno")))
  # CIZ delisting information: delisting action and reason, carried by the security's active records (the record
  # that starts after the delisting is inactive and has none); used only for the variant with GFKSW's imputed
  # returns for performance-related delistings
  dli <- setDT(dbGetQuery(con, paste("SELECT permno, secinfoenddt, securityenddt, securityactiveflg, delactiontype, delreasontype",
                                     "FROM stksecurityinfohist")))
  dli <- dli[securityactiveflg == "Y"][order(permno, secinfoenddt)][, .SD[.N], by = permno]
  saveRDS(dli, "cache/delist_info.rds"); rm(dli)
  dbDisconnect(con)
  ds[, ":="(date = to_idate(dlycaldt), del = dlydelflg %in% "Y")]
  ds[, ":="(mid = year(date) * 12L + month(date) - 1L, dlycaldt = NULL, dlydelflg = NULL)]
  dsh <- ds[held, on = .(permno, mid), nomatch = 0L, .(permno, date, mid, ret = dlyret, del)]
  dnx <- ds[del == TRUE, .(permno, mid = mid - 1L, ret_nd = dlyret, date_nd = date)][held, on = .(permno, mid), nomatch = 0L]
  setkey(dsh, permno, date)
  saveRDS(dsh, "cache/daily_held.rds"); saveRDS(dnx, "cache/daily_next_del.rds"); rm(ds, dsh, dnx, held); invisible(gc())
  logmsg("daily moments")
  d <- fd[, .(date, mx = mkt, ml = mkt_lag)][d[, .(permno, date, mid, ret_excess)], on = "date"]
  d <- d[is.finite(ret_excess) & is.finite(mx) & is.finite(ml)]
  d[, ":="(xy = mx * ret_excess, xx = mx * mx, ll = ml * ml, xl = mx * ml, ly = ml * ret_excess)]
  agg <- d[, .(n = .N, sx = sum(mx), sy = sum(ret_excess), sl = sum(ml), sxx = sum(xx),
               sxy = sum(xy), sll = sum(ll), sxl = sum(xl), sly = sum(ly)), by = .(permno, mid)]
  rm(d); invisible(gc())
  cm <- add_daily_betas(cm, agg); rm(agg); invisible(gc())
  saveRDS(cm, "cache/panel.rds"); saveRDS(fac, "cache/fac.rds"); saveRDS(fd, "cache/fd.rds")
  saveRDS(rets_m, "cache/rets_m.rds"); saveRDS(data_check, "cache/data_check.rds"); saveRDS(cache_key, "cache/key.rds")
} else {
  cm <- readRDS("cache/panel.rds"); fac <- readRDS("cache/fac.rds"); fd <- readRDS("cache/fd.rds")
  rets_m <- readRDS("cache/rets_m.rds"); data_check <- readRDS("cache/data_check.rds")
}
logmsg("panel ready")

# ---- 1b. integrity checks: dates, timing of information, completeness of holding returns ----------
# book equity: December FYE of year y-1 is used from the June formation of year y (>= 6 months lag)
fy_rule <- function(m) (m %/% 12L) - fifelse(m %% 12L >= 5L, 1L, 2L)
stopifnot(fy_rule(ym(1990L, 5L)) == 1988L, fy_rule(ym(1990L, 6L)) == 1989L, fy_rule(ym(1990L, 12L)) == 1989L,
          fy_rule(ym(1991L, 1L)) == 1989L, all(cm$fy == fy_rule(cm$mid)))
stopifnot(mid_of(as.IDate("2012-12-31")) == ym(2012L, 12L), ym(2013L, 1L) == ym(2012L, 12L) + 1L,
          format(eom_of_mid(ym(2012L, 12L))) == "2012-12-31", format(date_of_mid(ym(2013L, 1L))) == "2013-01-01")
# every eligible large-cap holding has a next-month return once delisting months come from msf_v2
# (the last formation month of the database has no next month and is never used as a holding month)
last_m <- max(cm$mid)
hold_check <- cm[large == TRUE & mid < last_m, .(holdings = .N, no_return = sum(!is.finite(ret_next)),
                                                  from_msf_only = sum(src_next == "msf_only"), no_trade = sum(src_next == "no_trade"),
                                                  no_trade_from_1967_12 = sum(src_next == "no_trade" & mid >= ym(1967L, 12L)),
                                                  no_trade_last = if (any(src_next == "no_trade")) max(mid[src_next == "no_trade"]) else NA_integer_)]
excl_check <- cm[large_all == TRUE & delisted_t == TRUE & mid < last_m,
                 .(excluded = .N, with_next_month_record = sum(is.finite(ret_next)))]
logmsg("holdings:", hold_check$holdings, "| without next-month return:", hold_check$no_return,
       "| return only in msf_v2:", hold_check$from_msf_only, "| no-trade months (zero return):", hold_check$no_trade,
       "| excluded at formation:", excl_check$excluded)
hold_check[, no_return_in_samples := cm[large == TRUE & mid < last_m & mid + 1L >= cfg$long[1], sum(!is.finite(ret_next))]]
# holdings without any return can only occur before the first return month analysed (January 1930)
stopifnot(hold_check$no_return_in_samples == 0L, excl_check$with_next_month_record == 0L)

# B/P outliers: GFKSW delete "the top and bottom 0.5% of the B/P observations" but keep all 540
# months (108 per state). Baseline: stock-level B/P trimmed cross-sectionally every month within
# the large-cap universe (uses no future information). Sensitivity: trimming over the pooled
# 1968-2012 panel (uses future information) and the first draft's month-level trimming.
cm[, bp_t := NA_real_]
cm[large == TRUE & is.finite(bp), bp_t := trim_cs(bp, cfg$bp_trim), by = mid]
cm[, bp_alt_t := NA_real_]                                         # sensitivity: book equity of calendar year y-1
cm[large == TRUE & is.finite(bp_alt), bp_alt_t := trim_cs(bp_alt, cfg$bp_trim), by = mid]
pl <- cm[large == TRUE & is.finite(bp) & in_period(mid + 1L, cfg$rep),
         quantile(bp, c(cfg$bp_trim, 1 - cfg$bp_trim), names = FALSE)]
cm[, bp_pool := fifelse(large == TRUE & is.finite(bp) & bp >= pl[1] & bp <= pl[2], bp, NA_real_)]
cm[large == TRUE, szp := frank(mktcap, ties.method = "average") / .N, by = mid]   # size percentile
# dollar-volume percentile within the universe (formation month; NA where CRSP reports no volume)
cm[large == TRUE, dvp := { r <- frank(dvol, ties.method = "average", na.last = "keep"); r / sum(is.finite(dvol)) }, by = mid]
BT <- cm[large == TRUE, c("permno", "mid", cfg$specs), with = FALSE]; setkey(BT, permno, mid)

# ---- 2. beta-sorted portfolios for every beta specification -------------------------------
quintile_series <- function(spec, ucol = "large") {
  h <- form_holdings(cm, spec, base_col = if (spec == "beta_M60") NULL else "beta_M60", ucol = ucol)
  stopifnot(h[mid < last_m & mid + 1L >= cfg$long[1], all(is.finite(ret_next))] || ucol == "large_all")
  h[, w := mktcap / sum(mktcap), by = .(mid, q)]
  r <- vw_returns(h)
  to <- rbindlist(lapply(1:5, function(k) turnover(h[q == k])[, q := k]))
  r[to, turnover := i.turnover, on = .(mid, q)]
  bp <- h[, .(bp = mean(bp_t, na.rm = TRUE), bp_raw = mean(bp, na.rm = TRUE), bp_pool = mean(bp_pool, na.rm = TRUE),
              bp_alt = mean(bp_alt_t, na.rm = TRUE),
              bp_vw = { ok <- is.finite(bp_t); sum(mktcap[ok] * bp_t[ok]) / sum(mktcap[ok]) },
              bp_cov = mean(is.finite(bp)), nform = .N), by = .(mid, q)]
  r[bp, ":="(bp = i.bp, bp_raw = i.bp_raw, bp_pool = i.bp_pool, bp_alt = i.bp_alt, bp_vw = i.bp_vw, bp_cov = i.bp_cov,
             nform = i.nform), on = .(mid, q)]
  # value-weighted forecasts of every beta estimator for this sort's quintiles (formation month)
  hb <- BT[h[, .(permno, mid, q, mktcap)], on = .(permno, mid)]
  bx <- hb[, lapply(.SD, function(b) { ok <- is.finite(b); if (!any(ok)) NA_real_ else sum(mktcap[ok] * b[ok]) / sum(mktcap[ok]) }),
           by = .(mid, q), .SDcols = cfg$specs]
  w <- dcast(r, mid ~ q, value.var = c("ret", "beta", "n", "turnover", "bp", "bp_raw", "bp_pool", "bp_alt", "bp_vw",
                                       "bp_cov", "nform"))
  setnames(w, gsub("_([1-5])$", "\\1", names(w)))
  w[, ":="(spread = bp1 - bp5, spread_raw = bp_raw1 - bp_raw5, spread_pool = bp_pool1 - bp_pool5,
           spread_alt = bp_alt1 - bp_alt5, spread_vw = bp_vw1 - bp_vw5)]
  w[, fmid := mid][, mid := mid + 1L]              # return month = formation + 1
  w <- fac[w, on = "mid"]
  w[, s_dyn := pmin(pmax(beta1 / beta5, 0), 1)]
  # ex-ante short ratio with betas shrunk towards one as in Frazzini and Pedersen (2014): 0.6 beta + 0.4
  w[, s_shr := pmin(pmax((0.6 * beta1 + 0.4) / (0.6 * beta5 + 0.4), 0), 1)]
  # excess returns: L-H is self-financing; BN = (r1 - rf) - 0.25 (r5 - rf) = r1 - 0.25 r5 - 0.75 rf
  w[, ":="(LH = ret1 - ret5, BN = ret1 - cfg$bn_short * ret5 - (1 - cfg$bn_short) * rf,
           BN_raw = ret1 - cfg$bn_short * ret5, BNd = ret1 - s_dyn * ret5 - (1 - s_dyn) * rf,
           BNs = ret1 - s_shr * ret5 - (1 - s_shr) * rf)]
  for (v in paste0("turnover", 1:5)) set(w, which(is.na(w[[v]])), v, 0)
  w[, ds := abs(s_dyn - shift(s_dyn, fill = s_dyn[1]))]
  w[, ds_shr := abs(s_shr - shift(s_shr, fill = s_shr[1]))]
  c <- cfg$base_cost; f <- cfg$borrow / 12
  w[, ":="(LH_net = LH - c * (turnover1 + turnover5) - f,
           BN_net = BN - c * (turnover1 + cfg$bn_short * turnover5) - cfg$bn_short * f,
           BNd_net = BNd - c * (turnover1 + s_dyn * turnover5 + ds) - s_dyn * f,
           BNs_net = BNs - c * (turnover1 + s_shr * turnover5 + ds_shr) - s_shr * f)]
  list(series = w, holdings = h, bx = bx)
}
Q <- list(); BX <- list(); HOLD <- NULL
for (s in cfg$specs) {
  logmsg("portfolios for", s)
  qs <- quintile_series(s)
  Q[[s]] <- qs$series; BX[[s]] <- qs$bx
  if (s == "beta_M60") HOLD <- qs$holdings
}
S <- Q[["beta_M60"]]
saveRDS(Q, "cache/quintiles.rds")

# ---- 3. replication: Table 1 --------------------------------------------------------------
# t: ordinary t-statistic, as in GFKSW's Table 1 (replication only)
table1 <- function(S, p) {
  x <- S[in_period(mid, p)]
  one <- function(name, raw, ex, b_ante) {
    rb <- suppressWarnings(nw_reg(ex, cbind(mkt = x$mkt)))[term == "mkt"]   # realised beta, Newey-West standard error
    data.table(strategy = name, arith = mean(raw), geo = prod(1 + raw)^(1 / length(raw)) - 1,
               beta_ante = b_ante, beta_real = rb$est, beta_real_se = rb$se,
               sharpe = sqrt(12) * mean(ex) / sd(ex), t = tstat(raw), n = length(raw))
  }
  rbind(
    rbindlist(lapply(1:5, function(k) one(paste0("Q", k), x[[paste0("ret", k)]],
                                          x[[paste0("ret", k)]] - x$rf, mean(x[[paste0("beta", k)]])))),
    one("Low - high", x$LH, x$LH, mean(x$beta1 - x$beta5)),
    one("Beta neutral", x$BN_raw, x$BN, mean(x$beta1 - cfg$bn_short * x$beta5)),
    one("Beta neutral (ex ante)", x$ret1 - x$s_dyn * x$ret5, x$BNd, mean(x$beta1 - x$s_dyn * x$beta5)),
    one("Market", x$mkt + x$rf, x$mkt, 1))
}

# ---- valuation states (quintiles of formation months on the B/P spread) --------------------
make_states <- function(S, p, breaks = NULL, col = "spread", trim_months = FALSE) {
  x <- S[in_period(mid, p) & is.finite(get(col)), .(mid, fmid, spread = get(col))]
  if (is.null(breaks)) {
    if (trim_months) { lim <- quantile(x$spread, c(0.005, 0.995)); x <- x[spread >= lim[1] & spread <= lim[2]] }
    breaks <- quantile(x$spread, seq(0, 1, 0.2), names = FALSE)
    setorder(x, spread, mid)
    x[, state := as.integer(ceiling(5 * seq_len(.N) / .N))]
  } else {
    x[, state := findInterval(spread, breaks[2:5]) + 1L]
  }
  setorder(x, mid)
  list(states = x, breaks = breaks)
}

# Forward compounding of monthly excess returns over return months mid .. mid+h-1 ("chained
# excess return"), annualised geometrically for h > 1. Only complete windows inside period p are
# used, so no return after p[2] enters (asserted below).
forward_returns <- function(S, st, p, cols = c("LH", "BN", "BNd"), hs = cfg$horizons) {
  x <- S[in_period(mid, p)]
  setkey(x, mid)
  rbindlist(lapply(hs, function(h) {
    rbindlist(lapply(cols, function(cl) {
      v <- x[[cl]]; m <- x$mid
      lg <- log1p(v); cs <- c(0, cumsum(lg))
      i <- match(st$mid, m); j <- i + h - 1L
      ok <- !is.na(i) & j <= length(m)
      ok[ok] <- (m[j[ok]] - m[i[ok]]) == (h - 1L)
      cum <- rep(NA_real_, nrow(st)); last <- rep(NA_integer_, nrow(st))
      cum[ok] <- exp(cs[j[ok] + 1] - cs[i[ok]]) - 1; last[ok] <- m[j[ok]]
      ann <- if (h == 1L) cum else (1 + cum)^(12 / h) - 1
      data.table(mid = st$mid, state = st$state, strategy = cl, h = h, value = ann, last = last)
    }))
  }))[is.finite(value)]
}

panelA <- function(S, st) {
  x <- merge(st$states[, .(mid, state, sspread = spread)], S, by = "mid")
  x[, .(from = min(sspread), to = max(sspread), months = .N, nstocks = mean((n1 + n5) / 2),
        Q1 = mean(ret1), Q2 = mean(ret2), Q3 = mean(ret3), Q4 = mean(ret4), Q5 = mean(ret5),
        LH = mean(LH), t_LH = tstat(LH), BN = mean(BN), t_BN = tstat(BN)), keyby = state]
}

# Tables 3/4 of GFKSW: regressions over windows t+1..t+h for every formation month, averaged by
# state. t_gfk: Newey-West with h lags on the state-filtered sequence of window estimates, our reading
# of GFKSW's procedure (replication only). That sequence treats the last month of one episode and the
# first month of the next as adjacent; the block bootstrap below gives standard errors for the same
# estimand, and ct_table() an alternative calendar-time estimator.
window_regressions <- function(S, st, p, cols = c("LH", "BN"), hs = c(12L, 24L, 36L, 48L, 60L)) {
  x <- S[in_period(mid, p)]; setkey(x, mid)
  X4 <- cbind(1, as.matrix(x[, .(mkt, smb, hml, mom)])); X1 <- X4[, 1:2]
  jobs <- CJ(i = seq_len(nrow(st)), h = hs)
  f <- function(r) {
    i0 <- match(st$mid[jobs$i[r]], x$mid); h <- jobs$h[r]
    if (is.na(i0) || i0 + h - 1L > nrow(x)) return(NULL)
    idx <- i0:(i0 + h - 1L)
    if (x$mid[i0 + h - 1L] - x$mid[i0] != h - 1L) return(NULL)
    rbindlist(lapply(cols, function(cl) {
      y <- x[[cl]][idx]
      b1 <- .lm.fit(X1[idx, , drop = FALSE], y)$coefficients
      b4 <- .lm.fit(X4[idx, , drop = FALSE], y)$coefficients
      data.table(mid = st$mid[jobs$i[r]], state = st$state[jobs$i[r]], h = h, strategy = cl,
                 a_capm = b1[1], a_4f = b4[1], b_hml = b4[4], b_mom = b4[5])
    }))
  }
  res <- parallel::mclapply(seq_len(nrow(jobs)), f, mc.cores = CORES)
  rbindlist(res)
}
summ_windows <- function(W) {
  W[, .(a_capm = mean(a_capm), t_gfk_capm = nw_t(a_capm, max(h)), a_4f = mean(a_4f),
        t_gfk_4f = nw_t(a_4f, max(h)), b_hml = mean(b_hml), b_mom = mean(b_mom), n = .N),
    keyby = .(strategy, state, h)]
}
# calendar-time WLS alternative (different estimand; one observation per calendar month, Newey-West with 12 lags)
ct_table <- function(S, st, p, cols = c("LH", "BN"), hs = c(12L, 36L, 60L)) {
  x <- S[in_period(mid, p)]; setkey(x, mid)
  rbindlist(lapply(cols, function(cl) rbindlist(lapply(hs, function(h) rbindlist(lapply(1:5, function(s) {
    starts <- st[state == s, mid]
    a <- ct_reg(x, starts, h, cl, character(0), cfg$ct_lag)
    if (is.null(a)) return(NULL)
    b <- ct_reg(x, starts, h, cl, "mkt", cfg$ct_lag)
    d <- ct_reg(x, starts, h, cl, c("mkt", "smb", "hml", "mom"), cfg$ct_lag)
    data.table(strategy = cl, h = h, state = s, mean = a$est, se_mean = a$se, t_mean = a$t,
               a_capm = b[term == "alpha", est], t_capm = b[term == "alpha", t],
               a_4f = d[term == "alpha", est], t_4f = d[term == "alpha", t],
               months = a$months, windows = a$windows)
  }))))))
}

logmsg("replication tables")
res <- list(cfg = cfg)
res$T1 <- table1(S, cfg$rep)
# robustness: only large caps with a December-FYE book-to-price
res$T1_bpuniv <- table1(quintile_series("beta_M60", ucol = "univ")$series, cfg$rep)
st_rep <- make_states(S, cfg$rep)
res$T2A <- panelA(S, st_rep)
res$T2_breaks <- st_rep$breaks
stopifnot(nrow(st_rep$states) == 540L, all(st_rep$states[, .N, by = state]$N == 108L),
          min(st_rep$states$fmid) == ym(1967L, 12L), max(st_rep$states$fmid) == ym(2012L, 11L))
fw_rep <- forward_returns(S, st_rep$states, cfg$rep)
stopifnot(max(fw_rep$last) <= cfg$rep[2])                         # no post-2012 return in replication
res$T2BC <- fw_rep[, .(value = mean(value), n = .N, last = max(last)), keyby = .(strategy, state, h)]
W_rep <- window_regressions(S, st_rep$states, cfg$rep)
res$T34 <- summ_windows(W_rep)
# block-bootstrap standard errors for exactly these state averages (state calendar held fixed)
logmsg("block bootstrap of the state averages: B =", cfg$boot_B, "| block =", cfg$boot_block, "months")
xb <- S[in_period(mid, cfg$rep)]; setkey(xb, mid)
stopifnot(nrow(xb) == 540L, all(diff(xb$mid) == 1L))
pos_state <- st_rep$states[match(xb$mid, mid), state]
RNGkind("L'Ecuyer-CMRG"); set.seed(cfg$seed)
BI <- replicate(cfg$boot_B, mbb_index(nrow(xb), cfg$boot_block))       # drawn once: reproducible for any core count
BTS <- state_window_stats(xb, pos_state, hs = c(12L, 24L, 36L, 48L, 60L), cols = c("LH", "BN"), idx = BI, cores = CORES)
chk_pt <- merge(BTS[, .(strategy, h, state, ret, a_capm, a_4f)], res$T2BC[, .(strategy, h, state, value)], by = c("strategy", "h", "state"))
chk_pt <- merge(chk_pt, res$T34[, .(strategy, h, state, a1 = a_capm, a4 = a_4f)], by = c("strategy", "h", "state"))
stopifnot(chk_pt[, max(abs(ret - value), abs(a_capm - a1), abs(a_4f - a4))] < 1e-12)   # same estimator as the tables
res$T2BC[BTS, ":="(se_boot = i.se_ret, t_boot = value / i.se_ret), on = .(strategy, h, state)]
res$T34[BTS, ":="(se_boot_capm = i.se_capm, t_boot_capm = a_capm / i.se_capm, se_boot_4f = i.se_4f,
                 t_boot_4f = a_4f / i.se_4f), on = .(strategy, h, state)]
res$boot_info <- list(B = cfg$boot_B, block = cfg$boot_block, seed = cfg$seed)
res$boot_d51 <- copy(attr(BTS, "d51"))[, ":="(t_d_ret = d_ret / se_d_ret, t_d_capm = d_capm / se_d_capm, t_d_4f = d_4f / se_d_4f)]

# ---- checks of the inference for the state averages of returns and the central contrast (v4) ----------
# (i) block lengths 12-60 months for the null-model bootstrap (24 = the indices used above), (ii) a
# volatility-preserving null model (blocks of standardised returns rescaled with the local volatility of the
# original calendar position), (iii) the observed-relation bootstrap (blocks of formation months with their
# states), (iv) Monte Carlo error of the bootstrap standard errors, (v) serial dependence of the strategy
# returns and volatility by state, (vi) size of the tests when returns are dependent (simulation).
logmsg("inference checks: block lengths, bootstrap variants, dependence")
bstats <- list()
for (cl in c("LH", "BN")) for (h in c(12L, 60L)) {
  y <- xb[[cl]]; a <- win_ret(y, h); n <- length(a); s <- pos_state[seq_len(n)]; th <- state_means(a, s)
  for (b in c(12L, 24L, 36L, 60L)) {
    I <- if (b == cfg$boot_block) BI else { set.seed(cfg$seed + b); replicate(cfg$boot_B, mbb_index(nrow(xb), b)) }
    bstats[[length(bstats) + 1L]] <- boot_summary(th, boot_null(y, pos_state, h, I))[, ":="(method = "null", block = b, strategy = cl, h = h)]
  }
  if (h == 12L) {
    bstats[[length(bstats) + 1L]] <- boot_summary(th, boot_vp(y, pos_state, h, BI))[, ":="(method = "vp", block = cfg$boot_block, strategy = cl, h = h)]
    set.seed(cfg$seed + 1000L); Ip <- replicate(cfg$boot_B, mbb_index(n, cfg$boot_block))
    bstats[[length(bstats) + 1L]] <- boot_summary(th, boot_pairs(a, s, Ip))[, ":="(method = "pairs", block = cfg$boot_block, strategy = cl, h = h)]
  }
}
res$boot_checks <- rbindlist(bstats)
chk_b <- merge(res$boot_checks[method == "null" & block == cfg$boot_block & stat == "d51", .(strategy, h, se)], res$boot_d51[, .(strategy, h, se_d_ret)], by = c("strategy", "h"))
stopifnot(chk_b[, max(abs(se - se_d_ret))] < 1e-12)            # same draws as the main bootstrap
# (vii) the volatility-preserving null model for the window alphas of Table 3 as well (same block indices BI;
# strategy and factor returns standardised with their own local volatility)
BTSv <- state_window_stats(xb, pos_state, hs = c(12L, 24L, 36L, 48L, 60L), cols = c("LH", "BN"), idx = BI, cores = CORES, vp = TRUE)
stopifnot(BTSv[chk_pt, on = .(strategy, h, state), max(abs(ret - i.ret), abs(a_capm - i.a_capm), abs(a_4f - i.a_4f))] < 1e-12)
chk_v <- merge(BTSv[h == 12L & state %in% c(1L, 5L), .(strategy, stat = paste0("state", state), se_ret)],
               res$boot_checks[method == "vp" & h == 12L, .(strategy, stat, se)], by = c("strategy", "stat"))
stopifnot(nrow(chk_v) == 4L, chk_v[, max(abs(se_ret - se))] < 1e-12)   # same estimator as boot_vp()
res$boot_vp_windows <- BTSv[, .(strategy, h, state, ret, se_ret, a_capm, se_capm, a_4f, se_4f,
                                t_ret = ret / se_ret, t_capm = a_capm / se_capm, t_4f = a_4f / se_4f)]
res$boot_vp_d51 <- copy(attr(BTSv, "d51"))[, ":="(t_d_ret = d_ret / se_d_ret, t_d_capm = d_capm / se_d_capm, t_d_4f = d_4f / se_d_4f)]
dep <- function(y, name) {
  ac <- acf(y, lag.max = 12, plot = FALSE)$acf[-1]; ac2 <- acf(y^2, lag.max = 12, plot = FALSE)$acf[-1]
  data.table(strategy = name, acf1 = ac[1], acf_maxabs = max(abs(ac)), acf_nsig = sum(abs(ac) > 2 / sqrt(length(y))),
             vr12 = vr_simple(y, 12), vr36 = vr_simple(y, 36), vr60 = vr_simple(y, 60), pw_block = pw_block(y),
             pw_block_sq = pw_block(y^2), acf2_1 = ac2[1], acf2_mean12 = mean(ac2))
}
res$dependence <- rbind(dep(xb$LH, "LH"), dep(xb$BN, "BN"))
gfit <- list(LH = garch_fit(xb$LH), BN = garch_fit(xb$BN))
res$garch <- rbindlist(lapply(names(gfit), function(cl) with(gfit[[cl]], data.table(strategy = cl, mu = mu, phi = phi, omega = omega,
                                                                                    alpha = alpha, beta = beta, nu = nu, conv = conv))))
res$state_vol <- rbindlist(lapply(c("LH", "BN"), function(cl) {
  y <- xb[[cl]]; n <- length(y) - 11L; s <- pos_state[seq_len(n)]; sg <- gfit[[cl]]$sigma
  rv <- sqrt(12) * vapply(seq_len(n), function(i) sd(y[i:(i + 11L)]), 0)
  gv <- sqrt(12) * vapply(seq_len(n), function(i) sqrt(mean(sg[i:(i + 11L)]^2)), 0)
  data.table(strategy = cl, state = 1:5, realised = as.numeric(tapply(rv, s, function(z) sqrt(mean(z^2)))),
             garch = as.numeric(tapply(gv, s, mean)), all_realised = sqrt(mean(rv^2)), all_garch = mean(gv))
}))
logmsg("simulation: size of the tests for the central contrast under dependent returns")
res$sim_contrast <- rbindlist(lapply(c("LH", "BN"), function(cl) {
  y <- xb[[cl]]; g <- gfit[[cl]]; TT <- length(y)
  gens <- list(iid = function() rnorm(TT, mean(y), sd(y)),
               garch = function() garch_sim(g, TT),
               observed_vol = function() { z <- rt_std(TT, g$nu) * g$sigma; out <- numeric(TT); out[1] <- g$mu + z[1]
                                           for (t in 2:TT) out[t] <- g$mu + g$phi * (out[t - 1] - g$mu) + z[t]; out })
  rbindlist(lapply(names(gens), function(nm) { set.seed(cfg$seed + 2L)
    sim_tests(pos_state, gens[[nm]], 12L, cfg$sim_n, cfg$sim_B, CORES)[, ":="(dgp = nm, strategy = cl)] }))
}))
# size of the three candidate tests under iid returns and the actual state calendar (simulation)
logmsg("simulation: size of the tests for state averages")
set.seed(cfg$seed + 1L)
res$sim_inference <- sim_state_inference(st_rep$states[order(mid), state], cores = CORES)
res$T34_ct <- ct_table(S, st_rep$states, cfg$rep)
res$states_rep <- st_rep$states[, .(mid, fmid, state)]
# does one episode dominate the central result? State assignments stay fixed; one decade of formation
# months is dropped at a time (pre-specified check of the 12-month state-5 minus state-1 difference
# and of the state-5 average window CAPM alpha; not used to choose any specification)
dec_of <- function(m) 10L * (((m - 1L) %/% 12L) %/% 10L)          # decade of the formation month (m = return month)
res$lodo <- rbindlist(lapply(c(NA_integer_, sort(unique(dec_of(fw_rep$mid)))), function(dd) {
  f <- if (is.na(dd)) fw_rep else fw_rep[dec_of(mid) != dd]
  w <- if (is.na(dd)) W_rep else W_rep[dec_of(mid) != dd]
  d51 <- function(stg) f[strategy == stg & h == 12, mean(value[state == 5]) - mean(value[state == 1])]
  data.table(dropped = dd, lh12 = d51("LH"), bn12 = d51("BN"),
             a5_lh = w[strategy == "LH" & h == 12 & state == 5, mean(a_capm)],
             a5_bn = w[strategy == "BN" & h == 12 & state == 5, mean(a_capm)],
             n5 = f[strategy == "LH" & h == 12, sum(state == 5)], n1 = f[strategy == "LH" & h == 12, sum(state == 1)])
}))

# sensitivity of the state classification and of the compounding convention
sens_def <- list(
  list(key = "base", col = "spread", trim = FALSE),
  list(key = "pool", col = "spread_pool", trim = FALSE),
  list(key = "month", col = "spread_raw", trim = TRUE),
  list(key = "none", col = "spread_raw", trim = FALSE),
  list(key = "vw", col = "spread_vw", trim = FALSE),          # value-weighted average B/P within the quintiles
  list(key = "fy", col = "spread_alt", trim = FALSE))         # book equity of calendar year y-1 in all months of y
res$state_sens <- rbindlist(lapply(sens_def, function(v) {
  st <- make_states(S, cfg$rep, col = v$col, trim_months = v$trim)
  agree <- merge(st$states[, .(mid, s1 = state)], st_rep$states[, .(mid, s0 = state)], by = "mid")[, mean(s1 == s0)]
  fw <- forward_returns(S, st$states, cfg$rep, cols = c("LH", "BN"), hs = c(12L, 60L))
  out <- fw[, .(value = mean(value)), keyby = .(strategy, h, state)]
  out[st$states[, .(months = .N), by = state], months := i.months, on = "state"]
  out[, ":="(variant = v$key, agree = agree)]
}))
conv <- function(hh) {
  x <- S[in_period(mid, cfg$rep)]; setkey(x, mid); st <- st_rep$states
  m <- x$mid; c1 <- c(0, cumsum(log1p(x$ret1))); c5 <- c(0, cumsum(log1p(x$ret5))); cl <- c(0, cumsum(log1p(x$LH)))
  cf <- c(0, cumsum(log1p(x$rf))); ct <- c(0, cumsum(log1p(x$rf + x$LH)))      # T-bills; collateralised L-H (total return)
  i <- match(st$mid, m); j <- i + hh - 1L; ok <- !is.na(i) & j <= length(m)
  W1 <- exp(c1[j[ok] + 1] - c1[i[ok]]); W5 <- exp(c5[j[ok] + 1] - c5[i[ok]]); WL <- exp(cl[j[ok] + 1] - cl[i[ok]])
  WF <- exp(cf[j[ok] + 1] - cf[i[ok]]); WT <- exp(ct[j[ok] + 1] - ct[i[ok]])
  data.table(state = st$state[ok], h = hh, chained = WL^(12 / hh) - 1, collateral = WT^(12 / hh) - WF^(12 / hh),
             relative = (W1 / W5)^(12 / hh) - 1, wealth_diff = W1 - W5)[, lapply(.SD, mean), keyby = .(state, h)]
}
res$conv_sens <- rbind(conv(12L), conv(60L))

# ---- data diagnostics: delisting months, formation exclusions, impact of the v3 return treatment ----
# (i) holdings whose next-month return exists only in msf_v2 (mostly delisting months) and
# (ii) stocks that delisted in the formation month (excluded at formation), as weights of the quintiles
ha <- form_holdings(cm, "beta_M60", ucol = "large_all")           # v2 universe: includes (ii)
ha[cm, delisted_t := i.delisted_t, on = .(permno, mid)]
ha[, w := mktcap / sum(mktcap), by = .(mid, q)]
dq <- function(p) ha[in_period(mid + 1L, p), .(w_msf = sum(w[src_next == "msf_only"]), n_msf = sum(src_next == "msf_only"),
                                                 w_excl = sum(w[delisted_t]), n_excl = sum(delisted_t)), by = .(mid, q)][
  , .(w_msf = mean(w_msf), n_msf = sum(n_msf), w_excl = mean(w_excl), n_excl = sum(n_excl)), keyby = q]
flags <- cm[HOLD[src_next == "msf_only", .(permno, mid)], on = .(permno, mid), .(delflg = delflg_next, retflg = retflg_next)]
# alternative return treatments of the same beta sort (series built like S, states of the baseline)
alt_series <- function(h) {
  r <- dcast(vw_returns(h), mid ~ q, value.var = "ret")
  setnames(r, as.character(1:5), paste0("ret", 1:5))
  r[, ":="(fmid = mid, mid = mid + 1L)]
  r <- fac[r, on = "mid"]
  r[, ":="(LH = ret1 - ret5, BN = ret1 - cfg$bn_short * ret5 - (1 - cfg$bn_short) * rf)]
  setkey(r, mid); r
}
h_v2 <- copy(ha)[src_next %in% c("msf_only", "no_trade"), ret_next := NA_real_]   # v2: dropped, weights renormalised
h_cash <- copy(ha)[delisted_t == TRUE, ret_next := 0]              # (ii) kept, merger proceeds held at a zero return
# (iii) delisting months in which CRSP's monthly file excludes the delisting-day return that its daily file records:
# the daily file's compounded return (delisting day included) replaces the monthly return
dsh <- readRDS("cache/daily_held.rds"); dnx <- readRDS("cache/daily_next_del.rds")
smd <- dsh[, .(r_dsf = prod(1 + ret[is.finite(ret)]) - 1, has_d = any(del)), by = .(permno, mid)]
smd[cm[, .(permno, mid = mid + 1L, rm = ret_next)], rm := i.rm, on = .(permno, mid)]
conflicts <- smd[has_d == TRUE & is.finite(rm) & abs((1 + rm) / (1 + r_dsf) - 1) > 1e-4, .(permno, mid, rm, r_dsf)]
h_dsf <- copy(HOLD)[conflicts[, .(permno, mid = mid - 1L, r_dsf)], ret_next := i.r_dsf, on = .(permno, mid)]   # baseline holdings
# (iv) GFKSW's rule: "If the delisting was for performance-related reasons ... we set the delisting return equal to
# -55% if trading on NASDAQ or to -30% if trading on the NYSE or AMEX" (their p. 48). Performance-related is
# approximated by CIZ delisting action 'GDR' (dropped by the exchange; the legacy 500-series codes). The held
# month's return becomes (1 + return before the delisting day, from the daily file) x (1 + imputed return) - 1.
dli <- readRDS("cache/delist_info.rds")
gd <- HOLD[, .(permno, mid, q, mktcap, ret_next)][cm[, .(permno, mid, delflg_next, exchange)], on = .(permno, mid), nomatch = 0L][
  !is.na(delflg_next) & delflg_next != "N"]
gd[dli, ":="(act = i.delactiontype, reason = i.delreasontype), on = "permno"]
gd_all <- copy(gd)
gd <- gd[act %chin% "GDR"]
pre <- dsh[del == FALSE & is.finite(ret), .(r_pre = prod(1 + ret) - 1), by = .(permno, mid)]
gd[pre[, .(permno, mid = mid - 1L, r_pre)], r_pre := i.r_pre, on = .(permno, mid)]
gd[is.na(r_pre), r_pre := 0]
gd[, r_gfksw := (1 + r_pre) * (1 + fifelse(exchange == "NASDAQ", -0.55, -0.30)) - 1]
h_gfk <- copy(HOLD)[gd[, .(permno, mid, r_gfksw)], ret_next := i.r_gfksw, on = .(permno, mid)]
variants <- list(v3 = S, v2 = alt_series(h_v2), cash = alt_series(h_cash), dsf = alt_series(h_dsf), gfksw = alt_series(h_gfk))
imp_row <- function(x, p, nm) {
  x <- x[in_period(mid, p)]
  data.table(variant = nm, q1 = mean(x$ret1), q5 = mean(x$ret5), lh = mean(x$LH), bn = mean(x$BN),
             vol1 = sd(x$ret1), vol5 = sd(x$ret5), mdd1 = max_drawdown(x$ret1), mdd5 = max_drawdown(x$ret5),
             mdd_lh = max_drawdown(x$LH + x$rf), months = nrow(x))
}
imp_states <- function(x, nm) {
  f <- forward_returns(x, st_rep$states, cfg$rep, cols = c("LH", "BN"), hs = c(12L, 60L))
  f[, .(value = mean(value)), keyby = .(strategy, h, state)][, variant := nm]
}
mk <- cm[is.finite(ret_next) & is.finite(mktcap) & mktcap > 0, .(vwret = sum(mktcap * ret_next) / sum(mktcap)), by = mid][, mid := mid + 1L][fac, on = "mid", nomatch = 0L][in_period(mid, cfg$rep)]
res$data_diag <- list(
  check_cm_vs_msf = data_check, hold_check = hold_check, excl_check = excl_check,
  q_rep = dq(cfg$rep), q_full = dq(cfg$full),
  flags = flags[, .N, keyby = .(delflg, retflg)],
  impact = rbindlist(lapply(c("rep", "full"), function(k) rbindlist(lapply(names(variants), function(v)
    imp_row(variants[[v]], cfg[[k]], v)))[, period := k])),
  impact_states = rbindlist(lapply(names(variants), function(v) imp_states(variants[[v]], v))),
  mkt_panel = mean(mk$vwret), mkt_ff = mean(mk$mkt + mk$rf), mkt_cor = cor(mk$vwret, mk$mkt + mk$rf),
  bp_trim_pooled = pl,
  gfksw_cases = merge(gd[, .(permno, mid, q, exchange, reason, ret_next, r_pre, r_gfksw)],
                      HOLD[, .(permno, w = mktcap / sum(mktcap)), by = .(mid, q)], by = c("permno", "mid", "q")),
  delist_actions = gd_all[, .N, keyby = act],
  conflicts = merge(conflicts, copy(HOLD)[, w := mktcap / sum(mktcap), by = .(mid, q)][, .(permno, mid = mid + 1L, q, w)],
                    by = c("permno", "mid"), all.x = TRUE))

# longer evidence 1930-2012 (B/P from Davis-Fama-French book equity before Compustat)
st_long <- make_states(S, cfg$long)
fw_long <- forward_returns(S, st_long$states, cfg$long, cols = "LH")
stopifnot(max(fw_long$last) <= cfg$long[2])
res$T5 <- fw_long[, .(value = mean(value), n = .N), keyby = .(state, h)]
res$T5_breaks <- st_long$breaks
res$T5A <- panelA(S, st_long)

# time series for figures
# B/P levels are measured at the end of the formation month (fmid); returns belong to month mid = fmid + 1
bp_all <- cm[is.finite(bp), .(bp_market = mean(bp)), by = .(fmid = mid)]
bp_large <- cm[large == TRUE & is.finite(bp_t), .(bp_large = mean(bp_t)), by = .(fmid = mid)]
res$series <- S[, .(mid, fmid, ret1, ret2, ret3, ret4, ret5, LH, BN, BN_raw, BNd, BNs, mkt, rf, smb, hml,
                    mom, bp1, bp5, spread, spread_raw, spread_vw, spread_alt, beta1, beta5, s_dyn, s_shr, n1, n5, bp_cov1, bp_cov5)][
  bp_all, bp_market := i.bp_market, on = "fmid"][bp_large, bp_large := i.bp_large, on = "fmid"]
stopifnot(all(res$series$mid == res$series$fmid + 1L))

# rolling 48-month correlations of the low-minus-high portfolio with HML and MOM
rc <- S[is.finite(LH) & is.finite(mom), .(mid, LH, hml, mom)]
rc[, ":="(rho_hml = frollapply(seq_len(.N), 48, function(i) cor(LH[i], hml[i])),
          rho_mom = frollapply(seq_len(.N), 48, function(i) cor(LH[i], mom[i])))]
res$fig4 <- rc[, .(mid, rho_hml, rho_mom)]

# ---- 4. after the original sample, after publication, and in real time ----------------------
logmsg("post-sample and post-publication")
res$T1_post <- table1(S, cfg$post); res$T1_pub <- table1(S, cfg$pub); res$T1_full <- table1(S, cfg$full)
perf_block <- function(p) {
  x <- S[in_period(mid, p)]
  v <- list(Q1 = x$ret1 - x$rf, Q5 = x$ret5 - x$rf, LH = x$LH, BN = x$BN, BNd = x$BNd, BNs = x$BNs, MKT = x$mkt)
  rbindlist(lapply(names(v), function(nm) {
    e <- v[[nm]]; ps <- perf_stats(e + x$rf, x$rf, x)
    data.table(portfolio = nm, mean_ex = ps$mean_ex, t_ex = ps$t_ex, sharpe = ps$sharpe, beta = ps$beta,
               a_capm = ps$alpha_capm, t_capm = ps$t_capm, a_4f = ps$alpha_4f, t_4f = ps$t_4f, months = length(e))
  }))
}
res$perf_periods <- rbindlist(lapply(c("rep", "post", "pub", "full"), function(k) perf_block(cfg[[k]])[, period := k]))
res$sr_q1_mkt <- rbindlist(lapply(c("rep", "post", "pub", "full"), function(k) {
  x <- S[in_period(mid, cfg[[k]])]; z <- sharpe_diff_test(x$ret1 - x$rf, x$mkt)
  data.table(period = k, diff = z[["diff"]], se = z[["se"]], p = z[["p"]])
}))
st_post <- make_states(S, cfg$post, breaks = st_rep$breaks)   # breakpoints from 1968-2012 only
st_pub  <- make_states(S, cfg$pub,  breaks = st_rep$breaks)
res$state_counts <- rbind(st_post$states[, .(period = "post", months = .N), keyby = state],
                          st_pub$states[,  .(period = "pub",  months = .N), keyby = state])
res$T2_post <- forward_returns(S, st_post$states, cfg$post, hs = c(1L, 12L, 36L))[
  , .(value = mean(value), n = .N), keyby = .(strategy, state, h)]
res$T2_pub <- forward_returns(S, st_pub$states, cfg$pub, hs = c(1L, 12L, 36L))[
  , .(value = mean(value), n = .N), keyby = .(strategy, state, h)]
st_full <- make_states(S, cfg$full)
res$T2_full <- forward_returns(S, st_full$states, cfg$full)[, .(value = mean(value), n = .N), keyby = .(strategy, state, h)]
res$T2A_full <- panelA(S, st_full)

# real-time states: expanding-window quintile breakpoints from earlier formation months only.
# The spread itself uses cross-sectionally trimmed B/P, so no step uses future information.
rt <- S[in_period(mid, cfg$full) & is.finite(spread), .(mid, fmid, spread)]
rt[, state := NA_integer_]
for (i in seq_len(nrow(rt))) if (i > cfg$rt_min) {
  br <- quantile(rt$spread[1:(i - 1)], seq(0.2, 0.8, 0.2), names = FALSE)
  rt$state[i] <- findInterval(rt$spread[i], br) + 1L
}
rt <- rt[!is.na(state)]
res$rt_start <- min(rt$mid)
res$T2_rt <- forward_returns(S, rt, cfg$full, hs = c(1L, 12L, 36L))[, .(value = mean(value), n = .N), keyby = .(strategy, state, h)]
res$T2_rt_ct <- ct_table(S, rt, cfg$full, hs = 12L)
res$rt_counts <- rbind(rt[in_period(mid, cfg$post), .(period = "post", months = .N), keyby = state],
                       rt[in_period(mid, cfg$pub), .(period = "pub", months = .N), keyby = state])

# predictive regressions of future 12-month LH/BN returns on the standardised B/P spread:
# overlapping monthly observations with Newey-West (12 lags), and non-overlapping December
# formations (one per year, OLS) as a check on the overlap correction.
# market excess return over the 12 months up to the formation month (known at formation), for return month mid
M12 <- S[order(mid), .(mid, mkt12 = exp(frollsum(log1p(shift(mkt)), 12L)) - 1)]
pred <- function(p, cl) {
  st <- S[in_period(mid, p) & is.finite(spread), .(mid, fmid, state = 0L, spread)]
  fw <- forward_returns(S, st, p, cols = cl, hs = 12L)[st, on = "mid", nomatch = 0L]
  fw[M12, mkt12 := i.mkt12, on = "mid"]
  z <- (fw$spread - mean(fw$spread)) / sd(fw$spread)
  r <- nw_reg(fw$value, cbind(spread_z = z), lag = 12L)
  a <- fw[fmid %% 12L == 11L]; za <- (a$spread - mean(fw$spread)) / sd(fw$spread)
  ra <- summary(lm(a$value ~ za))$coefficients
  # Stambaugh bias: the spread is persistent and its innovations are correlated with the returns; the
  # reduced-bias estimator of Amihud and Hurvich (2004) on the non-overlapping December observations
  x1 <- (S[match(a$fmid + 12L, fmid), spread] - mean(fw$spread)) / sd(fw$spread)
  ok <- is.finite(x1)
  ah <- pred_ah(a$value[ok], za[ok], x1[ok])
  # valuation signal or market cycle? the same regression controlling for the market's trailing 12-month return
  zm <- (fw$mkt12 - mean(fw$mkt12)) / sd(fw$mkt12)
  rc <- nw_reg(fw$value, cbind(spread_z = z, mkt12_z = zm), lag = 12L)
  am <- fw$fmid %% 12L == 11L
  rca <- summary(lm(fw$value[am] ~ z[am] + zm[am]))$coefficients
  data.table(strategy = cl, slope = r[term == "spread_z", est], t = r[term == "spread_z", t],
             ctrl_slope = rc[term == "spread_z", est], ctrl_t = rc[term == "spread_z", t],
             mkt12_slope = rc[term == "mkt12_z", est], mkt12_t = rc[term == "mkt12_z", t],
             ctrl_t_annual = rca[2, 3], mkt12_t_annual = rca[3, 3], cor_spread_mkt12 = cor(z, zm),
             r2 = summary(lm(fw$value ~ z))$r.squared, n = nrow(fw),
             slope_annual = ra[2, 1], t_annual = ra[2, 3], n_annual = nrow(a),
             ah_n = ah$n, ah_b_ols = ah$b_ols, ah_t_ols = ah$t_ols, ah_b = ah$b_ah, ah_t = ah$t_ah,
             ah_rho = ah$rho, ah_rho_c = ah$rho_c, ah_corr_uv = ah$corr_uv, ah_bias = ah$stambaugh_bias)
}
res$pred <- rbindlist(lapply(list(rep = cfg$rep, post = cfg$post, pub = cfg$pub, full = cfg$full), function(p)
  rbind(pred(p, "LH"), pred(p, "BN"))), idcol = "sample")

# focused strategy comparison: does conditioning on valuation (and momentum) in real time help?
# Rule taken from GFKSW's published result (alphas concentrated in the two highest B/P spread
# states): hold the long-short portfolio when the real-time state is 4 or 5, otherwise T-bills.
# Momentum variant: additionally require a positive trailing 12-month return of the strategy itself
# (time-series momentum of the strategy, known at formation; not stock momentum and not the MOM
# factor). Signals use only information available at formation, but the rule was written down by
# us after the data existed, so no period is an untouched out-of-sample test.
logmsg("timing comparison")
tim <- S[in_period(mid, cfg$full)][, .(mid, LH, BN, LH_net, BN_net, rf, mkt, smb, hml, mom)]
tim[rt, state_rt := i.state, on = "mid"]
setorder(tim, mid)
gross_notional <- c(LH = 2, BN = 1 + cfg$bn_short)
timing_rows <- function(p) {
  out <- rbindlist(lapply(c("LH", "BN"), function(cl) {
    x <- copy(tim)
    cs <- c(0, cumsum(log1p(x[[cl]])))
    n <- nrow(x); tr <- rep(NA_real_, n)
    if (n > 12) tr[13:n] <- exp(cs[13:n] - cs[1:(n - 12)]) - 1   # chained return over months t-12 .. t-1
    x[, trail := tr]
    x[, ":="(inv_always = 1, inv_val = as.numeric(state_rt %in% cfg$timing_states),
             inv_vm = as.numeric(state_rt %in% cfg$timing_states & trail > 0))]
    x <- x[in_period(mid, p) & !is.na(state_rt) & is.finite(trail)]
    rbindlist(lapply(c("always", "val", "vm"), function(k) {
      inv <- x[[paste0("inv_", k)]]
      g <- inv * x[[cl]]
      nt <- inv * x[[paste0(cl, "_net")]] - cfg$base_cost * gross_notional[[cl]] * abs(inv - shift(inv, fill = inv[1]))
      ps <- perf_stats(g + x$rf, x$rf, x); pn <- perf_stats(nt + x$rf, x$rf, x)
      data.table(strategy = cl, rule = k, invested = mean(inv), switches = sum(abs(diff(inv))),
                 mean_ex = ps$mean_ex, t_ex = ps$t_ex, sharpe = ps$sharpe, a_capm = ps$alpha_capm, t_capm = ps$t_capm,
                 a_4f = ps$alpha_4f, t_4f = ps$t_4f, sharpe_net = pn$sharpe, months = nrow(x), g = list(g), nt = list(nt))
    }))
  }))
  out[, c("sr_diff", "sr_se", "p_sr") := {
    base <- nt[[which(rule == "always")]]
    rbindlist(lapply(seq_len(.N), function(j) if (rule[j] == "always") data.table(NA_real_, NA_real_, NA_real_) else
      as.data.table(as.list(unname(sharpe_diff_test(nt[[j]], base))))))
  }, by = strategy]
  out[!is.na(p_sr), p_holm := p.adjust(p_sr, "holm")]
  out[, c("g", "nt") := NULL]
  out
}
res$timing <- rbindlist(lapply(c("full", "post", "pub"), function(k) timing_rows(cfg[[k]])[, period := k]))

# ---- 5. extension A: beta measurement (frequency x window) --------------------------------
logmsg("extension A: beta specifications")
ext_spec <- function(spec) {
  x <- Q[[spec]][in_period(mid, cfg$full)]
  st <- make_states(Q[[spec]], cfg$full)
  fw <- forward_returns(Q[[spec]], st$states, cfg$full, cols = "LH", hs = 12L)
  dv <- fw[, mean(value[state == 5]) - mean(value[state == 1])]
  rb <- function(v) coef(lm(v - x$rf ~ x$mkt))[2]
  pB <- perf_stats(x$BNd + x$rf, x$rf, x)
  pN <- perf_stats(x$BNd_net + x$rf, x$rf, x)
  data.table(spec = spec, n_q = mean(x$n1), b1_ante = mean(x$beta1), b5_ante = mean(x$beta5),
             b1_real = rb(x$ret1), b5_real = rb(x$ret5),
             q1_sharpe = sqrt(12) * mean(x$ret1 - x$rf) / sd(x$ret1 - x$rf),
             lh = 12 * mean(x$LH), t_lh = nw_t(x$LH), bnd = 12 * mean(x$BNd), t_bnd = nw_t(x$BNd),
             bnd_sharpe = pB$sharpe, bnd_beta = pB$beta, bnd_a4 = pB$alpha_4f, bnd_t4 = pB$t_4f,
             to1 = 12 * mean(x$turnover1), to_lh = 12 * mean(x$turnover1 + x$turnover5),
             bnd_sharpe_net = pN$sharpe, val_effect = dv)
}
res$extA <- rbindlist(lapply(cfg$specs, ext_spec))
res$sr_bnd_vs_m60 <- rbindlist(lapply(setdiff(cfg$specs, "beta_M60"), function(spec) {
  a <- Q[[spec]][in_period(mid, cfg$full), .(mid, BNd)]; b <- S[in_period(mid, cfg$full), .(mid, BNd)]
  m <- merge(a, b, by = "mid"); z <- sharpe_diff_test(m$BNd.x, m$BNd.y)
  data.table(spec = spec, diff = z[["diff"]], se = z[["se"]], p = z[["p"]])
}))[, p_holm := p.adjust(p, "holm")]

# beta forecast quality. Target: realised beta of a quintile portfolio over the next 12 months
# (monthly returns). (i) own-sort RMSE (the forecast an investor in that sort relies on);
# (ii) Diebold-Mariano tests on COMMON targets: for a pair of estimators, the test assets are the
# quintiles of both sorts, so each estimator is evaluated on exactly the same portfolios.
logmsg("beta forecast accuracy")
RB <- rbindlist(lapply(cfg$specs, function(spec) {
  x <- Q[[spec]][in_period(mid, cfg$full)]; setkey(x, mid)
  rbindlist(lapply(1:5, function(k) {
    y <- x[[paste0("ret", k)]] - x$rf
    rr <- frollapply(seq_len(nrow(x)), 12, function(i) cov(y[i], x$mkt[i]) / var(x$mkt[i]), align = "left")
    data.table(sort = spec, q = k, fmid = x$fmid, real = rr)
  }))
}))
FB <- rbindlist(lapply(cfg$specs, function(s) copy(BX[[s]])[, sort := s]))
E <- merge(RB[is.finite(real)], FB, by.x = c("sort", "q", "fmid"), by.y = c("sort", "q", "mid"))
res$beta_fc <- rbindlist(lapply(cfg$specs, function(s) E[sort == s, .(spec = s, ante = mean(get(s), na.rm = TRUE), real = mean(real),
                                                                      rmse = sqrt(mean((real - get(s))^2, na.rm = TRUE))), keyby = q]))
# families: (a) each estimator vs GFKSW's M60 (Holm over six); (b) frequency at a common 36-month
# window (D36 vs M36); (c) window within the daily frequency (D3, D12 vs D36); (d) Dimson correction
pairs <- c(lapply(setdiff(cfg$specs, "beta_M60"), function(s) c(s, "beta_M60")),
           list(c("beta_D36", "beta_M36"), c("beta_D3", "beta_D36"), c("beta_D12", "beta_D36"), c("beta_DD12", "beta_D12")))
res$beta_dm <- rbindlist(lapply(pairs, function(pr) {
  z <- E[sort %in% pr & is.finite(get(pr[1])) & is.finite(get(pr[2]))]
  z[, ":="(e1 = (real - get(pr[1]))^2, e2 = (real - get(pr[2]))^2)]
  d <- z[, .(d = mean(e1 - e2)), by = fmid]
  tt <- dm_test(d$d, 12L)
  data.table(est1 = pr[1], est2 = pr[2], rmse1 = sqrt(mean(z$e1)), rmse2 = sqrt(mean(z$e2)),
             dm_t = tt[["t"]], dm_p = tt[["p"]], months = tt[["n"]])
}))
res$beta_dm[est2 == "beta_M60", p_holm := p.adjust(dm_p, "holm")]
res$beta_dm[est2 == "beta_D36" & est1 %in% c("beta_D3", "beta_D12"), p_holm := p.adjust(dm_p, "holm")]
# exploratory family: all ten forecast comparisons (with M60, frequency, window and Dimson), Holm-adjusted together
res$beta_dm[, p_holm_all := p.adjust(dm_p, "holm")]

# ---- 6. extension B: long-only implementation vs well-diversified portfolios --------------
logmsg("extension B: long-only portfolios")
base <- HOLD
lo <- list(); WTS <- list()
q1 <- vw(base[q == 1L]); WTS$LowBeta <- q1
lo$LowBeta <- port_from_weights(q1)
qd <- vw(form_holdings(cm, "beta_D12", base_col = "beta_M60")[q == 1L]); WTS[["LowBeta (daily 12m)"]] <- qd
lo[["LowBeta (daily 12m)"]] <- port_from_weights(qd)
bw <- buffer_weights(cm, "beta_M60", cfg$buffer[1], cfg$buffer[2]); WTS[["LowBeta buffer"]] <- bw
lo[["LowBeta buffer"]] <- port_from_weights(bw)
# rets_m: msf_v2 returns of all security-months (delisting months included), loaded in section 1
WTS[["LowBeta quarterly"]] <- q1[mid %% 3L == 2L]
lo[["LowBeta quarterly"]] <- port_periodic(q1[mid %% 3L == 2L, .(mid, permno, w)], rets_m, 3L)
# UCITS-compatible version: value weights capped at 5% per stock (excess redistributed pro rata), so no
# position exceeds 5% and the 5/10/40 rule cannot bind
q1c <- copy(q1)[, w := cap_weights(w, cfg$cap), by = mid]
stopifnot(q1c[, max(w)] <= cfg$cap + 1e-9, q1c[, abs(sum(w) - 1), by = mid][, max(V1)] < 1e-9)
WTS[["LowBeta capped"]] <- q1c; lo[["LowBeta capped"]] <- port_from_weights(q1c)
logmsg("minimum variance (", CORES, "cores )")
mv_months <- sort(unique(base$mid)); mv_months <- mv_months[mv_months >= cfg$full[1] - 1L]
MV <- minvar_weights(cm, mv_months, cfg$minvar_window, cfg$minvar_wmax, CORES)
saveRDS(MV, "cache/minvar_weights.rds")
WTS[["Min-variance"]] <- MV
lo[["Min-variance"]] <- port_from_weights(MV[, .(mid, permno, w, ret_next)])
# the same portfolio with the covariance matrix estimated from the daily returns of the last 12 months by
# analytical nonlinear shrinkage (Ledoit and Wolf 2020)
logmsg("minimum variance, nonlinear shrinkage on daily returns")
DU <- readRDS("cache/daily_univ.rds"); setkey(DU, mid)
MVN <- minvar_weights_nl(cm, DU, mv_months, cfg$minvar_wmax, CORES)
rm(DU); invisible(gc())
saveRDS(MVN, "cache/minvar_nl_weights.rds")
res$minvar_nl_missing <- setdiff(mv_months, MVN$mid)              # formation months without a solution (expected: none)
logmsg("NL minimum variance: months without a solution:", length(res$minvar_nl_missing))
WTS[["Min-variance NL"]] <- MVN
lo[["Min-variance NL"]] <- port_from_weights(MVN[, .(mid, permno, w, ret_next)])
res$minvar_nl_info <- MVN[, .(n = .N, days = days[1], wmax = max(w), eff_n = 1 / sum(w^2)), by = mid]
WTS[["Large-cap VW"]] <- vw(base); lo[["Large-cap VW"]] <- port_from_weights(vw(base))
WTS[["Large-cap EW"]] <- ew(base); lo[["Large-cap EW"]] <- port_from_weights(ew(base))
LO <- rbindlist(lo, idcol = "strategy", fill = TRUE)
LO[is.na(turnover), turnover := 0]
LO <- fac[LO, on = "mid"]
mkt_row <- fac[, .(strategy = "Market (FF)", mid, ret = mkt + rf, turnover = 0)]
LO <- rbind(LO, fac[mkt_row, on = "mid"], fill = TRUE)
res$LO <- LO[in_period(mid, c(cfg$long[1], cfg$full[2]))]
res$minvar_info <- MV[, .(n = .N, delta = delta[1], wmax = max(w), eff_n = 1 / sum(w^2)), by = mid]
# diversification and liquidity of the target weights at rebalancing dates, 1968-2024. Liquidity is
# descriptive only: weight-averaged dollar-volume percentile within the universe and monthly share
# turnover (volume / shares outstanding) in the formation month; 'liq_cov' is the weight share with
# volume data. No market-impact or borrowing-availability data are used or assumed.
wavg <- function(w, x) { ok <- is.finite(x); if (!any(ok)) NA_real_ else sum(w[ok] * x[ok]) / sum(w[ok]) }
res$div <- rbindlist(lapply(names(WTS), function(nm) {
  W <- WTS[[nm]][in_period(mid + 1L, cfg$full), .(mid, permno, w)]
  W[cm, ":="(szp = i.szp, dvp = i.dvp, sturn = i.sturn), on = .(permno, mid)]
  lq <- W[, .(size_pct = wavg(w, szp), dv_pct = wavg(w, dvp), sturn = wavg(w, sturn),
              liq_cov = sum(w[is.finite(dvp)]) / sum(w)), by = mid]
  cbind(strategy = nm, div_stats(W), lq[, lapply(.SD, mean), .SDcols = -"mid"])
}))

candidates <- c("LowBeta", "LowBeta buffer", "LowBeta quarterly", "LowBeta capped", "LowBeta (daily 12m)", "Min-variance",
                "Min-variance NL", "Large-cap EW")
lo_table <- function(p, cost = 0) {
  x <- LO[in_period(mid, p)]
  bm <- x[strategy == "Large-cap VW", .(mid, bench = ret - cost * turnover)]
  x <- bm[x, on = "mid"]
  x[, net := ret - cost * turnover]
  out <- x[, {
    ps <- perf_stats(net, rf, .SD, bench = bench)
    bmex <- bench - rf
    sd_test <- if (strategy[1] == "Large-cap VW") c(NA_real_, NA_real_, NA_real_) else sharpe_diff_test(net - rf, bmex)
    cbind(ps, to = 12 * mean(turnover), sr_diff = sd_test[1], sr_se = sd_test[2], p_sr = sd_test[3])
  }, by = strategy]
  out[strategy %in% candidates, p_holm := p.adjust(p_sr, "holm")]
  out
}
res$LO_full <- lo_table(cfg$full); res$LO_full_net <- lo_table(cfg$full, cfg$base_cost)
res$LO_post <- lo_table(cfg$post); res$LO_post_net <- lo_table(cfg$post, cfg$base_cost)
res$LO_pub_net <- lo_table(cfg$pub, cfg$base_cost)
res$LO_rep  <- lo_table(cfg$rep)
# cost sensitivity and break-even costs relative to the large-cap VW benchmark
res$LO_costs <- rbindlist(lapply(cfg$costs, function(cst) lo_table(cfg$full, cst)[, cost := cst]))
# (i) cost that eliminates the gross CAPM alpha; (ii) cost that equalises the net Sharpe ratio
#     with the (also cost-bearing) large-cap VW benchmark
be_tab <- function(p) {
  x <- LO[in_period(mid, p)]
  bm <- x[strategy == "Large-cap VW", .(mid, bret = ret, bto = turnover)]
  x <- bm[x, on = "mid"]
  x[, {
    a <- unname(coef(lm(I(ret - rf) ~ mkt))[1])
    be_alpha <- if (a > 0 && mean(turnover) > 0) 1e4 * a / mean(turnover) else NA_real_
    gap <- function(cc) { e <- ret - cc * turnover - rf; eb <- bret - cc * bto - rf
                          mean(e) / sd(e) - mean(eb) / sd(eb) }
    be_sr <- if (strategy[1] %in% c("Large-cap VW", "Market (FF)") || gap(0) <= 0) NA_real_ else
      if (gap(0.05) > 0) Inf else 1e4 * uniroot(gap, c(0, 0.05))$root
    .(alpha_gross = 12 * a, to = 12 * mean(turnover), be_alpha = be_alpha, be_sr = be_sr)
  }, by = strategy]
}
res$LO_breakeven <- be_tab(cfg$full)
res$LO_breakeven_post <- be_tab(cfg$post)
# long-short strategies gross vs net
ls_tab <- function(p) {
  x <- S[in_period(mid, p)]
  f <- function(name, g, n, to) {
    a <- perf_stats(g + x$rf, x$rf, x); b <- perf_stats(n + x$rf, x$rf, x)
    data.table(strategy = name, mean_g = a$mean_ex, t_g = a$t_ex, sharpe_g = a$sharpe, beta = a$beta,
               a1_g = a$alpha_capm, t1_g = a$t_capm, a4_g = a$alpha_4f, t4_g = a$t_4f,
               mean_n = b$mean_ex, sharpe_n = b$sharpe, a1_n = b$alpha_capm, t1_n = b$t_capm,
               a4_n = b$alpha_4f, t4_n = b$t_4f, to = 12 * mean(to),
               be_capm = 1e4 * (a$alpha_capm / 12) / mean(to))
  }
  rbind(f("Low - high", x$LH, x$LH_net, x$turnover1 + x$turnover5),
        f("Beta neutral (fixed 25%)", x$BN, x$BN_net, x$turnover1 + cfg$bn_short * x$turnover5),
        f("Beta neutral (ex ante)", x$BNd, x$BNd_net, x$turnover1 + x$s_dyn * x$turnover5 + x$ds),
        f("Beta neutral (shrunk beta)", x$BNs, x$BNs_net, x$turnover1 + x$s_shr * x$turnover5 + x$ds_shr))
}
res$LS_full <- ls_tab(cfg$full); res$LS_rep <- ls_tab(cfg$rep); res$LS_post <- ls_tab(cfg$post); res$LS_pub <- ls_tab(cfg$pub)
res$turnover <- S[in_period(mid, cfg$full), .(q1 = 12 * mean(turnover1), q5 = 12 * mean(turnover5))]

# ---- 7. extension C: evaluation frequency (daily vs monthly) and rolling windows -----------
logmsg("extension C: daily portfolio returns")
# Daily paths of the same monthly rebalanced portfolios (same holdings and month-end weights; weights drift
# within the month). Baseline: CRSP's daily file dsf_v2 with its delisting-day records on their own dates;
# delisting-day records dated in the first days of the next month, which CRSP's monthly file credits to the
# held month, are booked on the held month's last trading day; in the few months in which the monthly file
# excludes an in-month delisting-day return, the monthly file is followed (see daily_panel_dsf). Variants:
# the daily file's delisting returns in those months ('daily_file'); the v3 reconstruction from crsp_daily
# with the whole monthly-daily difference booked on the last trading day ('v3'); crsp_daily without any
# delisting returns ('no_delisting').
dret_cd <- readRDS("cache/daily_held_cd.rds")
cal <- fd[, .(date, mid = year(date) * 12L + month(date) - 1L)]
ports <- list("LowBeta" = q1, "HighBeta" = vw(base[q == 5L]), "Min-variance" = MV[, .(mid, permno, w, ret_next)],
              "Large-cap VW" = vw(base), "Large-cap EW" = ew(base))
recon <- list(); dvar <- list(); moved <- 0L
for (nm in names(ports)) {
  pm <- daily_panel_dsf(ports[[nm]], cal, dsh, dnx, "monthly")
  recon[[nm]] <- pm$sm[, strategy := nm]; moved <- moved + attr(pm$sm, "moved")
  # every stock-month of the baseline compounds to its monthly return
  stopifnot(pm$g[, .(r = prod(1 + ret) - 1, rm = rm[1]), by = .(mid, permno)][, max(abs((1 + r) / (1 + rm) - 1))] < 1e-10)
  dvar[[nm]] <- list(base = daily_port_from_panel(pm$g),
                     daily_file = daily_port_from_panel(daily_panel_dsf(ports[[nm]], cal, dsh, dnx, "daily")$g),
                     v3 = daily_port_from_panel(daily_panel_cd(ports[[nm]], cal, dret_cd, TRUE)),
                     no_delisting = daily_port_from_panel(daily_panel_cd(ports[[nm]], cal, dret_cd, FALSE)))
  rm(pm); invisible(gc())
}
RC <- rbindlist(recon)[in_period(mid, cfg$full)]
res$recon <- RC[, .(stock_months = .N, weight_month = sum(w) / uniqueN(RC$mid)), keyby = .(strategy, cls)]
res$recon_cases <- RC[cls %in% c("conflict", "unexplained"), .(strategy, mid, permno, w, rm, r_in, r_ex)]
res$recon_moved <- moved
# long-short positions are also rebalanced only at month ends: within the month both legs drift.
# L-H (GFKSW note 9): capital of USD 1 in Q1, USD 1 of Q5 sold short, short proceeds in T-bills, so the position
# is worth G1 - G5 + Gf (cumulative gross returns since the last rebalancing). BN: capital in Q1, short
# 25% of capital in Q5, short proceeds in T-bills: G1 - 0.25 G5 + 0.25 Gf. The month's compounded daily
# total return is therefore rf + r1 - r5 and r1 - 0.25 r5 + 0.25 rf, as in the monthly engine.
build_D <- function(v) {
  D <- rbindlist(lapply(names(dvar), function(nm) dvar[[nm]][[v]][, strategy := nm]))
  D <- fd[D, on = "date"]
  lh <- dcast(D[strategy %in% c("LowBeta", "HighBeta")], date + mid + rf + mkt ~ strategy, value.var = "ret")
  setorder(lh, date)
  lh[, ":="(G1 = cumprod(1 + LowBeta), G5 = cumprod(1 + HighBeta), Gf = cumprod(1 + rf)), by = mid]
  lh[, ":="(V_lh = Gf + G1 - G5, V_bn = G1 - cfg$bn_short * G5 + cfg$bn_short * Gf)]
  lh[, ":="(r_lh = V_lh / shift(V_lh, fill = 1) - 1, r_bn = V_bn / shift(V_bn, fill = 1) - 1), by = mid]
  lh <- rbind(lh[, .(strategy = "Low - high", date, mid, rf, mkt, ret = r_lh)],
              lh[, .(strategy = "Beta neutral (fixed 25%)", date, mid, rf, mkt, ret = r_bn)])
  D <- rbind(D[, .(strategy, date, mid, rf, mkt, ret)], lh,
             fd[, .(strategy = "Market (FF)", date, mid = year(date) * 12L + month(date) - 1L, rf, mkt, ret = mkt + rf)])
  D[in_period(mid, cfg$full) & is.finite(ret)]
}
freq_stats <- function(r, rf, k) {
  ex <- r - rf
  data.table(mean_ex = k * mean(ex), vol = sqrt(k) * sd(r), sharpe = sqrt(k) * mean(ex) / sd(ex), mdd = max_drawdown(r),
             skew = mean((ex - mean(ex))^3) / sd(ex)^3, kurt = mean((ex - mean(ex))^4) / sd(ex)^4 - 3,
             var95 = -quantile(r, 0.05), es95 = -mean(r[r <= quantile(r, 0.05)]))
}
D <- build_D("base")
# monthly market excess return: compounded market total return minus compounded T-bill return
DM <- D[, .(ret = prod(1 + ret) - 1, rf = prod(1 + rf) - 1, mkt = prod(1 + mkt + rf) - prod(1 + rf)), by = .(strategy, mid)]
res$freq <- rbind(D[, freq_stats(ret, rf, 252), by = strategy][, freq := "daily"],
                  DM[, freq_stats(ret, rf, 12), by = strategy][, freq := "monthly"])
# market beta at both frequencies with Newey-West standard errors (rule-of-thumb lags)
beta_se <- function(r, rf, mkt) { z <- suppressWarnings(nw_reg(r - rf, cbind(mkt = mkt))); z[term == "mkt", .(beta = est, se = se)] }   # market on itself: exact fit
res$freq_beta <- merge(D[, beta_se(ret, rf, mkt), by = strategy][, .(strategy, beta_daily = beta, se_daily = se)],
                       DM[, beta_se(ret, rf, mkt), by = strategy][, .(strategy, beta_monthly = beta, se_monthly = se)], by = "strategy")
stopifnot(res$freq_beta[strategy == "Market (FF)", abs(beta_daily - 1) < 1e-10 && abs(beta_monthly - 1) < 1e-10])   # market on itself
# direct test of serial correlation within the month: variance ratio of 21-day log returns (Lo-MacKinlay robust z)
res$vr <- D[, as.list(vr_test(ret, 21L)), by = strategy]
# the same statistics for the alternative daily reconstructions
res$daily_variants <- rbindlist(lapply(names(dvar[[1]]), function(v) {
  Dv <- if (v == "base") D else build_D(v)
  Dv[, c(freq_stats(ret, rf, 252)[, .(vol, sharpe, kurt, mdd)], .(beta = cov(ret - rf, mkt) / var(mkt))), by = strategy][, variant := v]
}))
# consistency: compounded daily path vs monthly engine, every portfolio with a monthly counterpart.
# Long-only portfolios and L-H (excess of the compounded daily T-bill rate) must agree to rounding error;
# BN and the market differ only through monthly versus compounded daily T-bill and factor data.
eng <- rbind(LO[in_period(mid, cfg$full) & strategy %in% c("LowBeta", "Min-variance", "Large-cap VW", "Large-cap EW", "Market (FF)"),
                .(strategy, mid, rm = ret, rf_m = rf)],
             S[in_period(mid, cfg$full), .(strategy = "HighBeta", mid, rm = ret5, rf_m = rf)],
             S[in_period(mid, cfg$full), .(strategy = "Low - high", mid, rm = LH, rf_m = 0)],
             S[in_period(mid, cfg$full), .(strategy = "Beta neutral (fixed 25%)", mid, rm = BN, rf_m = 0)])
DMx <- copy(DM)[strategy %in% c("Low - high", "Beta neutral (fixed 25%)"), ret := ret - rf]   # excess returns
chk <- merge(DMx[, .(strategy, mid, rd = ret)], eng, by = c("strategy", "mid"))
stopifnot(chk[strategy %in% c("LowBeta", "HighBeta", "Min-variance", "Large-cap VW", "Large-cap EW", "Low - high"),
              max(abs(rd - rm))] < 1e-10)
res$daily_check <- chk[, {
  a <- abs(rd - rm); j <- which.max(a)
  .(cor = cor(rd, rm), mad_bp = 1e4 * median(a), p99_bp = 1e4 * quantile(a, 0.99, names = FALSE), max_bp = 1e4 * a[j],
    max_mid = mid[j], share_10bp = mean(a > 0.001), n = .N,
    vol_d = sqrt(12) * sd(rd), vol_m = sqrt(12) * sd(rm), sr_d = sqrt(12) * mean(rd - rf_m) / sd(rd - rf_m),
    sr_m = sqrt(12) * mean(rm - rf_m) / sd(rm - rf_m), mdd_d = max_drawdown(rd), mdd_m = max_drawdown(rm),
    mean_d = 12 * mean(rd - rf_m), mean_m = 12 * mean(rm - rf_m))
}, by = strategy]
# why the monthly columns of the frequency table can differ from the long-only table: same monthly returns,
# but T-bill rates compounded from daily data instead of the monthly T-bill series
t78 <- merge(LO[in_period(mid, cfg$full) & strategy %in% c(names(ports), "Market (FF)"), .(strategy, mid, rm = ret, rf_m = rf)],
             DM[, .(strategy, mid, rd = ret, rf_d = rf)], by = c("strategy", "mid"))
res$t7t8 <- t78[, .(sr_t7 = sqrt(12) * mean(rm - rf_m) / sd(rm - rf_m), sr_t8 = sqrt(12) * mean(rd - rf_d) / sd(rd - rf_d),
                    sr_t8_monthly_rf = sqrt(12) * mean(rd - rf_m) / sd(rd - rf_m), vol_t7 = sqrt(12) * sd(rm), vol_t8 = sqrt(12) * sd(rd),
                    max_ret_diff = max(abs(rd - rm)), max_rf_diff = max(abs(rf_d - rf_m))), by = strategy]
res$daily_extremes <- D[, .SD[order(-abs(ret))][1:3], by = strategy][, .(strategy, date, ret)]
# rolling 36- and 60-month evaluation windows (monthly), low beta vs large-cap VW benchmark
logmsg("rolling windows")
roll <- function(win) {
  x <- dcast(LO[in_period(mid, cfg$full) & strategy %in% c("LowBeta", "Large-cap VW", "Min-variance")],
             mid + rf + mkt + smb + hml + mom ~ strategy, value.var = "ret")
  setnames(x, c("Large-cap VW", "Min-variance"), c("VW", "MV"))
  bn <- S[, .(mid, BN, LH)]; x <- bn[x, on = "mid"]
  setorder(x, mid)
  f <- function(i) {
    z <- x[i]
    sr <- function(v) sqrt(12) * mean(v - z$rf) / sd(v - z$rf)
    a <- function(v) 12 * unname(coef(lm((v - z$rf) ~ z$mkt))[1])
    c(sr_lb = sr(z$LowBeta), sr_vw = sr(z$VW), sr_mv = sr(z$MV), a_lb = a(z$LowBeta), a_mv = a(z$MV),
      a_bn = a(z$BN + z$rf), sr_bn = sr(z$BN + z$rf),
      act_lb = 12 * mean(z$LowBeta - z$VW),                     # gross active return vs large-cap VW
      vol_lb = sqrt(12) * sd(z$LowBeta), vol_vw = sqrt(12) * sd(z$VW), vol_mv = sqrt(12) * sd(z$MV),
      vol_bn = sqrt(12) * sd(z$BN))
  }
  idx <- lapply(win:nrow(x), function(e) (e - win + 1):e)
  out <- as.data.table(do.call(rbind, lapply(idx, f)))
  out[, ":="(mid = x$mid[win:nrow(x)], window = win)]
  out
}
res$rolling <- rbind(roll(36L), roll(60L))

res$meta <- list(
  built = Sys.time(), db = DB, runtime_min = as.numeric(difftime(Sys.time(), t_start, units = "mins")),
  panel_months = range(cm$mid), n_large = cm[large == TRUE, .N], n_univ = cm[univ == TRUE, .N],
  n_large_rep = cm[large == TRUE & mid >= cfg$rep[1] - 1L & mid < cfg$rep[2], .N / uniqueN(mid)],
  n_univ_rep = cm[univ == TRUE & mid >= cfg$rep[1] - 1L & mid < cfg$rep[2], .N / uniqueN(mid)],
  bp_cov_rep = S[in_period(mid, cfg$rep), mean((bp_cov1 + bp_cov5) / 2)],
  bp_cov_long = S[in_period(mid, cfg$long), mean((bp_cov1 + bp_cov5) / 2)],
  universe_n = S[in_period(mid, cfg$rep), mean(n1 + n2 + n3 + n4 + n5)],
  panel_built_this_run = panel_built, cache_key = cache_key, db_size = db_info$size, db_mtime = db_info$mtime,
  r_version = R.version.string,
  packages = sapply(c("data.table", "DBI", "RSQLite", "sandwich", "quadprog"), function(p) as.character(packageVersion(p))),
  session = capture.output(sessionInfo()))
saveRDS(res, "results/res.rds")
logmsg("DONE in", round(res$meta$runtime_min, 1), "minutes")
