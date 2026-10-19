# =============================================================================
# Low-Volatility Cycles: replication and extensions -- function library (v4)
# All returns are decimals. Months are integer ids: mid = 12 * year + month - 1
# (January 2013 = 2013 * 12 + 0). In portfolio series, 'fmid' is the formation month
# (characteristics measured at its end) and 'mid' = fmid + 1 the return month.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(DBI); library(RSQLite); library(sandwich); library(quadprog)
})

logmsg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n", sep = " ")

# ---- dates -------------------------------------------------------------------
to_idate <- function(x) {
  if (inherits(x, "Date")) return(as.IDate(x))
  if (is.character(x)) return(as.IDate(x))
  as.IDate(as.integer(x))                      # integer days since 1970-01-01
}
mid_of <- function(d) { d <- as.IDate(d); year(d) * 12L + month(d) - 1L }
date_of_mid <- function(m) as.Date(sprintf("%d-%02d-01", m %/% 12L, m %% 12L + 1L))   # first day of month m
eom_of_mid <- function(m) date_of_mid(m + 1L) - 1L                                     # last day of month m
ym <- function(y, m) y * 12L + m - 1L

# ---- rolling sums with partial windows at the start of each series -------------
rsum <- function(v, k) {
  if (!length(v) || k == 1L) return(v)
  out <- frollsum(c(numeric(k - 1L), v), k, algo = "exact")
  out[-seq_len(k - 1L)]
}

# ---- data loading ----------------------------------------------------------------
load_monthly <- function(con) {
  cm <- setDT(dbGetQuery(con, paste(
    "SELECT permno, date, ret, ret_excess, prc, mktcap, exchange, gvkey",
    "FROM crsp_monthly")))
  cm[, mid := mid_of(to_idate(date))][, date := NULL]
  cm <- cm[exchange %in% c("NYSE", "AMEX", "NASDAQ")]
  setkey(cm, permno, mid)
  stopifnot(!anyDuplicated(cm, by = c("permno", "mid")))
  cm
}

# Monthly returns of ALL CRSP security-months (msf_v2, CIZ format), including the delisting month.
# crsp_monthly (Tidy Finance) keeps only months whose month-end falls in a security-info record with
# trading status 'A'; the month in which a stock delists therefore drops out, although msf_v2 reports
# its return including the delisting payment (mthretflg 'DE'). On rows present in both tables the
# returns are identical (checked in build.R). Also returns the delisting flag and volume fields.
load_msf <- function(con) {
  m <- setDT(dbGetQuery(con, paste(
    "SELECT permno, mthcaldt, mthret, mthdelflg, mthretflg, mthvol, mthprcvol, shrout FROM msf_v2")))
  m[, mid := mid_of(to_idate(mthcaldt))][, mthcaldt := NULL]
  setkey(m, permno, mid)
  stopifnot(!anyDuplicated(m, by = c("permno", "mid")))
  m
}

load_factors <- function(con) {
  ff <- setDT(dbGetQuery(con, "SELECT date, mkt_excess, smb, hml, risk_free FROM factors_ff3_monthly"))
  mo <- setDT(dbGetQuery(con, "SELECT date, mom FROM factors_momentum_monthly"))
  ff[, mid := mid_of(to_idate(date))]; mo[, mid := mid_of(to_idate(date))]
  fac <- merge(ff[, .(mid, mkt = mkt_excess, smb, hml, rf = risk_free)], mo[, .(mid, mom)],
               by = "mid", all.x = TRUE)
  setkey(fac, mid)
  fac
}

load_daily_factors <- function(con) {
  fd <- setDT(dbGetQuery(con, "SELECT date, mkt_excess, risk_free FROM factors_ff3_daily"))
  fd[, date := to_idate(date)]
  setkey(fd, date)
  fd[, ":="(mkt = mkt_excess, rf = risk_free, mkt_lag = shift(mkt_excess))]
  fd[, .(date, mkt, mkt_lag, rf)]
}

# ---- book equity (Fama-French 1992 via stored be; December FYE; two years) ----
load_book_equity <- function(con, dff_file = NULL) {
  ca <- setDT(dbGetQuery(con, "SELECT gvkey, datadate, be FROM compustat_annual"))
  ca[, datadate := to_idate(datadate)]
  setorder(ca, gvkey, datadate)
  ca[, nyears := seq_len(.N), by = gvkey]               # years of Compustat history
  ca <- ca[month(datadate) == 12L & nyears >= 2L & is.finite(be) & be > 0]
  ca[, fy := year(datadate)]
  ca <- unique(ca[, .(gvkey, fy, be)], by = c("gvkey", "fy"), fromLast = TRUE)
  dff <- NULL
  if (!is.null(dff_file) && file.exists(dff_file)) dff <- read_dff_be(dff_file)
  list(comp = ca, dff = dff)
}

# Davis-Fama-French historical book equity (Ken French data library), keyed by PERMNO.
# Record: permno, first year, last year, then book equity for 1926, 1927, ... (-99.99 = NA).
read_dff_be <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- trimws(ln[grepl("^ *[0-9]", ln)])
  parts <- strsplit(ln, "[[:space:]]+")
  out <- rbindlist(lapply(parts, function(p) {
    v <- as.numeric(p)
    if (length(v) < 4) return(NULL)
    vals <- v[-(1:3)]
    data.table(permno = as.integer(v[1]), fy = 1926L + seq_along(vals) - 1L, be = vals)
  }))
  out[be > 0 & be != -99.99]
}

# ---- monthly betas: OLS of stock excess return on market excess return --------------
# windows: months; minobs: minimum valid observations. Beta at month t uses t-w+1..t.
add_monthly_betas <- function(cm, fac, windows = c(36L, 60L, 120L), minobs = c(24L, 24L, 60L)) {
  x <- cm[, .(permno, mid, y = ret_excess)]
  g <- x[, .(mid = seq.int(min(mid), max(mid))), by = permno]
  x <- x[g, on = .(permno, mid)]
  x <- fac[, .(mid, mx = mkt)][x, on = "mid"]
  setorder(x, permno, mid)
  x[, ok := is.finite(y) & is.finite(mx)]
  x[, ":="(xx = fifelse(ok, mx, 0), yy = fifelse(ok, y, 0))]
  for (j in seq_along(windows)) {
    k <- windows[j]
    x[, ":="(n = rsum(as.numeric(ok), k), sx = rsum(xx, k), sy = rsum(yy, k),
             sxx = rsum(xx * xx, k), sxy = rsum(xx * yy, k)), by = permno]
    x[, b := fifelse(n >= minobs[j], (sxy - sx * sy / n) / (sxx - sx * sx / n), NA_real_)]
    setnames(x, "b", paste0("beta_M", k))
  }
  bcols <- paste0("beta_M", windows)
  cm[x[, c("permno", "mid", bcols), with = FALSE], (bcols) := mget(paste0("i.", bcols)),
     on = .(permno, mid)]
  cm
}

# ---- daily betas from monthly-aggregated daily moments ----------------------------------
# OLS beta and Dimson (1979) beta with one lag of the market.
daily_moments <- function(con, fd) {
  d <- setDT(dbGetQuery(con, "SELECT permno, date, ret_excess FROM crsp_daily"))
  d[, date := to_idate(date)]
  d <- fd[, .(date, mx = mkt, ml = mkt_lag)][d, on = "date"]
  d <- d[is.finite(ret_excess) & is.finite(mx) & is.finite(ml)]
  d[, mid := year(date) * 12L + month(date) - 1L]
  d[, ":="(xy = mx * ret_excess, xx = mx * mx, ll = ml * ml, xl = mx * ml, ly = ml * ret_excess)]
  agg <- d[, .(n = .N, sx = sum(mx), sy = sum(ret_excess), sl = sum(ml), sxx = sum(xx),
               sxy = sum(xy), sll = sum(ll), sxl = sum(xl), sly = sum(ly)), by = .(permno, mid)]
  rm(d); invisible(gc())
  agg
}

add_daily_betas <- function(cm, agg, windows = c(3L, 12L, 36L), minobs = c(50L, 190L, 560L),
                            dimson_window = 12L) {
  g <- agg[, .(mid = seq.int(min(mid), max(mid))), by = permno]
  a <- agg[g, on = .(permno, mid)]
  for (v in c("n", "sx", "sy", "sl", "sxx", "sxy", "sll", "sxl", "sly")) set(a, which(is.na(a[[v]])), v, 0)
  setorder(a, permno, mid)
  out <- a[, .(permno, mid)]
  for (j in seq_along(windows)) {
    k <- windows[j]
    r <- a[, lapply(.SD, rsum, k = k), by = permno, .SDcols = c("n", "sx", "sy", "sxx", "sxy")]
    b <- with(r, fifelse(n >= minobs[j], (sxy - sx * sy / n) / (sxx - sx * sx / n), NA_real_))
    set(out, j = paste0("beta_D", k), value = b)
    if (k == dimson_window) {
      r <- a[, lapply(.SD, rsum, k = k), by = permno,
             .SDcols = c("n", "sx", "sy", "sl", "sxx", "sxy", "sll", "sxl", "sly")]
      cxx <- r$sxx - r$sx^2 / r$n; cll <- r$sll - r$sl^2 / r$n; cxl <- r$sxl - r$sx * r$sl / r$n
      cxy <- r$sxy - r$sx * r$sy / r$n; cly <- r$sly - r$sl * r$sy / r$n
      det <- cxx * cll - cxl^2
      b0 <- (cll * cxy - cxl * cly) / det; b1 <- (cxx * cly - cxl * cxy) / det
      set(out, j = paste0("beta_DD", k), value = fifelse(r$n >= minobs[j], b0 + b1, NA_real_))
    }
  }
  bcols <- setdiff(names(out), c("permno", "mid"))
  cm[out, (bcols) := mget(paste0("i.", bcols)), on = .(permno, mid)]
  cm
}

# ---- book-to-price: December book equity of fiscal year y-1 from June of year y ---------
add_bp <- function(cm, be) {
  cm[, fy := (mid %/% 12L) - fifelse(mid %% 12L >= 5L, 1L, 2L)]   # month index 5 = June
  cm[be$comp, be_comp := i.be, on = .(gvkey, fy)]
  if (!is.null(be$dff)) {
    cm[be$dff, be_dff := i.be, on = .(permno, fy)]
  } else cm[, be_dff := NA_real_]
  cm[, be := fifelse(is.finite(be_comp), be_comp, be_dff)]
  cm[, bp := fifelse(is.finite(be) & mktcap > 0, be / mktcap, NA_real_)]
  # sensitivity ("prior fiscal year" read literally): December book equity of calendar year y-1 in every
  # formation month of year y. Differs from the baseline only in January-May, where it uses accounting
  # data that were not yet published (look-ahead of up to five months); used only to classify states.
  cm[, fy_alt := (mid %/% 12L) - 1L]
  bcomp <- be$comp; bdff <- be$dff                                 # 'be' is also a column of cm from here on
  cm[bcomp, be_comp2 := i.be, on = .(gvkey, fy_alt = fy)]
  if (!is.null(bdff)) cm[bdff, be_dff2 := i.be, on = .(permno, fy_alt = fy)] else cm[, be_dff2 := NA_real_]
  cm[, bp_alt := { b <- fifelse(is.finite(be_comp2), be_comp2, be_dff2); fifelse(is.finite(b) & mktcap > 0, b / mktcap, NA_real_) }]
  cm[, c("be_comp", "be_dff", "be_comp2", "be_dff2", "fy_alt") := NULL]
  cm
}

# ---- universe: NYSE/AMEX/NASDAQ, $2 <= price <= $1000, top third by NYSE breakpoints ----
# Next-month returns come from msf_v2 (all security-months, delisting months included), so a stock
# held at formation keeps its return in the month it delists. Eligibility at formation additionally
# requires that the security has not delisted by the formation date: msf_v2 flags a delisting in the
# formation month itself (mthdelflg != 'N') when the last trading day is in that month; CRSP then books
# the delisting payment in the formation month and the stock cannot be bought at the formation date.
add_universe <- function(cm, msf) {
  bp_nyse <- cm[exchange == "NYSE" & is.finite(mktcap) & mktcap > 0,
                .(size_cut = quantile(mktcap, 2 / 3, type = 7)), by = mid]
  cm[bp_nyse, size_cut := i.size_cut, on = "mid"]
  cm[msf, ":="(delflg = i.mthdelflg, dvol = i.mthprcvol, svol = i.mthvol, shr = i.shrout), on = .(permno, mid)]
  cm[, delisted_t := !is.na(delflg) & delflg != "N"]
  cm[, large_all := is.finite(mktcap) & mktcap > 0 & is.finite(prc) & abs(prc) >= 2 &
       abs(prc) <= 1000 & is.finite(size_cut) & mktcap >= size_cut]
  cm[, large := large_all & !delisted_t]
  # robustness universe: large caps with a valid December-FYE book-to-price
  cm[, univ := large & is.finite(bp)]
  nx <- msf[, .(permno, mid = mid - 1L, ret_next = mthret, retflg_next = mthretflg, delflg_next = mthdelflg)]
  cm[nx, ":="(ret_next = i.ret_next, retflg_next = i.retflg_next, delflg_next = i.delflg_next), on = .(permno, mid)]
  # source of the next-month return: a crsp_monthly row exists for t+1, or the month exists only in msf_v2
  in_cm <- cm[, .(permno, mid = mid - 1L)]
  cm[, src_next := fifelse(is.finite(ret_next), "msf_only", "none")]
  cm[in_cm, src_next := fifelse(is.finite(ret_next), "crsp_monthly", "none"), on = .(permno, mid)]
  # months without any trade (msf_v2 row without a return, flag 'NT', no delisting): the position is
  # valued at its last trade price, i.e. a zero return for that month (an assumption; counted in build.R)
  cm[src_next == "none" & retflg_next %in% "NT" & delflg_next %in% "N", ":="(ret_next = 0, src_next = "no_trade")]
  # liquidity at formation: dollar volume (sum of daily price x volume) and share turnover in month t
  cm[, ":="(dvol = fifelse(is.finite(dvol) & dvol > 0, dvol, NA_real_),
            sturn = fifelse(is.finite(svol) & svol > 0 & is.finite(shr) & shr > 0, svol / (1000 * shr), NA_real_))]
  cm[, c("delflg", "svol", "shr") := NULL]
  cm
}

# ---- quintile holdings for one beta specification ------------------------------------------
form_holdings <- function(cm, bcol, base_col = NULL, nq = 5L, ucol = "large") {
  u <- cm[get(ucol) == TRUE & is.finite(get(bcol))]
  if (!is.null(base_col)) u <- u[is.finite(get(base_col))]   # common universe across specs
  u[, beta := get(bcol)]
  setorder(u, mid, beta, permno)
  u[, q := as.integer(ceiling(nq * seq_len(.N) / .N)), by = mid]
  cols <- c("mid", "permno", "q", "beta", "mktcap", "bp", intersect(c("bp_t", "bp_pool", "bp_alt_t", "src_next"), names(u)), "ret_next")
  u[, ..cols]
}

vw <- function(h) { h <- h[, .(mid, permno, w = mktcap, ret_next)]; h[, w := w / sum(w), by = mid]; h }
ew <- function(h) { h <- h[, .(mid, permno, w = 1, ret_next)]; h[, w := w / sum(w), by = mid]; h }

# value weights at formation. After the msf_v2 fix every eligible holding has a next-month return
# (checked in build.R); the filter below only matters for the documented 'v2 treatment' comparison,
# in which recovered delisting-month returns are removed and the weights renormalised.
vw_returns <- function(h) {
  h[is.finite(ret_next), .(ret = sum(mktcap * ret_next) / sum(mktcap), n = .N,
                           beta = sum(mktcap * beta) / sum(mktcap)), by = .(mid, q)]
}

# Turnover at formation month t: sum |w_t - w~_t|, w~_t = weights of t-1 drifted to end of t.
# The sum counts purchases and sales (two-sided turnover); a one-way cost c per dollar traded is
# therefore charged as c * turnover. A stock that delisted during month t has a drifted weight but no
# target weight and is counted as a sale (conservative: cash merger proceeds need no trade).
# h: mid, permno, w (target weights at formation), ret_next.
turnover <- function(h) {
  h <- h[is.finite(ret_next) & w > 0]
  h[, w := w / sum(w), by = mid]
  h[, wd := w * (1 + ret_next) / sum(w * (1 + ret_next)), by = mid]
  A <- h[, .(mid = mid + 1L, permno, wd)]
  B <- h[, .(mid, permno, w)]
  M <- merge(A, B, by = c("mid", "permno"), all = TRUE)
  M[is.na(wd), wd := 0][is.na(w), w := 0]
  to <- M[mid %in% B$mid, .(turnover = sum(abs(w - wd))), by = mid]
  to[mid == min(mid), turnover := NA_real_]          # initial purchase not counted
  setkey(to, mid)
  to
}

# portfolio return series from target weights (monthly rebalancing)
port_from_weights <- function(h) {
  h <- h[is.finite(ret_next) & w > 0]
  r <- h[, .(ret = sum(w * ret_next) / sum(w), n = .N), by = mid]
  r[turnover(h), turnover := i.turnover, on = "mid"]
  r[, mid := mid + 1L]                                  # stamp with holding (return) month
  setkey(r, mid)
  r
}

# buy-and-hold between rebalancing dates (every "every" months), drifting weights.
# weights: data.table(mid, permno, w) at rebalancing months; rets: data.table(permno, mid, ret) from
# msf_v2, so the delisting month carries the delisting return; in later months of the holding period
# the delisted stock has no return and its proceeds are reinvested pro rata in the remaining holdings.
port_periodic <- function(weights, rets, every = 3L) {
  rb <- sort(unique(weights$mid))
  out <- vector("list", length(rb)); prev <- NULL
  for (j in seq_along(rb)) {
    t0 <- rb[j]; t1 <- if (j < length(rb)) rb[j + 1L] else t0 + every
    w <- weights[mid == t0, .(permno, w)]
    to <- NA_real_
    if (!is.null(prev)) {
      m <- merge(prev, w, by = "permno", all = TRUE)
      m[is.na(wd), wd := 0][is.na(w), w := 0]
      to <- sum(abs(m$w - m$wd))
    }
    v <- w$w; names(v) <- w$permno; rr <- list()
    for (k in seq_len(t1 - t0)) {
      rk <- rets[.(as.integer(names(v)), t0 + k), on = .(permno, mid), ret]
      keep <- is.finite(rk)
      if (!any(keep)) break
      v <- v[keep]; rk <- rk[keep]
      newv <- v * (1 + rk)
      rr[[k]] <- data.table(mid = t0 + k, ret = sum(newv) / sum(v) - 1,
                            turnover = if (k == 1L) to else 0)
      v <- newv
    }
    prev <- data.table(permno = as.integer(names(v)), wd = v / sum(v))
    out[[j]] <- rbindlist(rr)
  }
  res <- rbindlist(out); setkey(res, mid); res
}

# ---- buffer (hysteresis) low-beta portfolio: enter <= p_in, exit > p_out -------------
buffer_weights <- function(cm, bcol, p_in = 0.2, p_out = 0.3, ucol = "large") {
  u <- cm[get(ucol) == TRUE & is.finite(get(bcol)), .(mid, permno, b = get(bcol), mktcap, ret_next)]
  u[, pct := frank(b, ties.method = "first") / .N, by = mid]
  months <- sort(unique(u$mid)); held <- integer(0); res <- vector("list", length(months))
  for (j in seq_along(months)) {
    ut <- u[mid == months[j]]
    keep <- ut$permno %in% held & ut$pct <= p_out
    sel <- ut[ut$pct <= p_in | keep]
    held <- sel$permno
    res[[j]] <- sel[, .(mid, permno, w = mktcap / sum(mktcap), ret_next)]
  }
  rbindlist(res)
}

# ---- Ledoit-Wolf (2004) shrinkage towards constant correlation --------------------------
cov_shrink_cc <- function(X) {
  t <- nrow(X); n <- ncol(X)
  X <- sweep(X, 2, colMeans(X))
  S <- crossprod(X) / t
  v <- diag(S); sv <- sqrt(v)
  rbar <- (sum(S / outer(sv, sv)) - n) / (n * (n - 1))
  Fm <- rbar * outer(sv, sv); diag(Fm) <- v
  Y <- X^2
  phiMat <- crossprod(Y) / t - S^2
  phi <- sum(phiMat)
  thetaMat <- crossprod(X^3, X) / t - v * S
  diag(thetaMat) <- 0
  rho <- sum(diag(phiMat)) + rbar * sum(outer(1 / sv, sv) * thetaMat)
  gam <- sum((S - Fm)^2)
  delta <- max(0, min(1, (phi - rho) / gam / t))
  list(sigma = delta * Fm + (1 - delta) * S, delta = delta)
}

# long-only global minimum variance with an upper weight bound
gmv_long_only <- function(Sigma, wmax = 0.03) {
  n <- ncol(Sigma)
  if (n * wmax < 1) wmax <- 1.5 / n
  sc <- 1 / mean(diag(Sigma))
  sol <- solve.QP(Dmat = Sigma * sc, dvec = rep(0, n),
                  Amat = cbind(rep(1, n), diag(n), -diag(n)),
                  bvec = c(1, rep(0, n), rep(-wmax, n)), meq = 1)
  w <- pmax(sol$solution, 0); w / sum(w)
}

minvar_weights <- function(cm, months, window = 60L, wmax = 0.03, cores = 1L, bcol = "beta_M60",
                           ucol = "large") {
  R <- cm[, .(permno, mid, ret)]
  U <- cm[get(ucol) == TRUE & is.finite(get(bcol)), .(mid, permno)]
  f <- function(t) {
    ids <- U[mid == t, permno]
    x <- R[permno %in% ids & mid > t - window & mid <= t]
    W <- dcast(x, mid ~ permno, value.var = "ret")
    if (nrow(W) < window) return(NULL)
    M <- as.matrix(W[, -1]); M <- M[, colSums(is.finite(M)) == window, drop = FALSE]
    if (ncol(M) < 30) return(NULL)
    sh <- cov_shrink_cc(M)
    w <- gmv_long_only(sh$sigma, wmax)
    data.table(mid = t, permno = as.integer(colnames(M)), w = w, delta = sh$delta)
  }
  res <- if (cores > 1) parallel::mclapply(months, f, mc.cores = cores) else lapply(months, f)
  out <- rbindlist(res)
  out[cm, ret_next := i.ret_next, on = .(permno, mid)]
  out
}

# ---- analytical nonlinear shrinkage (Ledoit and Wolf 2020, Annals of Statistics 48(5)) -----------
# X: T x N matrix of returns (complete). Port of the authors' published Matlab code
# (analytical_shrinkage.m), including the case N > T - 1.
nl_shrink <- function(X) {
  n <- nrow(X); p <- ncol(X)
  X <- sweep(X, 2, colMeans(X)); n <- n - 1L
  S <- crossprod(X) / n
  e <- eigen(S, symmetric = TRUE)
  lambda <- rev(e$values); u <- e$vectors[, p:1, drop = FALSE]      # ascending
  lambda <- lambda[max(1L, p - n + 1L):p]                          # non-zero eigenvalues
  m <- length(lambda)
  L <- matrix(lambda, m, m)                                         # L[i, j] = lambda_i
  h <- n^(-1/3)
  H <- h * t(L)                                                     # H[i, j] = h * lambda_j
  x <- (L - t(L)) / H
  ftilde <- (3 / 4 / sqrt(5)) * rowMeans(pmax(1 - x^2 / 5, 0) / H)
  Hftemp <- (-3 / 10 / pi) * x + (3 / 4 / sqrt(5) / pi) * (1 - x^2 / 5) * log(abs((sqrt(5) - x) / (sqrt(5) + x)))
  edge <- abs(x) == sqrt(5); Hftemp[edge] <- (-3 / 10 / pi) * x[edge]
  Hftilde <- rowMeans(Hftemp / H)
  if (p <= n) {
    cr <- p / n
    dtilde <- lambda / ((pi * cr * lambda * ftilde)^2 + (1 - cr - pi * cr * lambda * Hftilde)^2)
  } else {
    Hftilde0 <- (1 / pi) * (3 / 10 / h^2 + 3 / 4 / sqrt(5) / h * (1 - 1 / 5 / h^2) *
                              log((1 + sqrt(5) * h) / (1 - sqrt(5) * h))) * mean(1 / lambda)
    dtilde0 <- 1 / (pi * (p - n) / n * Hftilde0)
    dtilde1 <- lambda / (pi^2 * lambda^2 * (ftilde^2 + Hftilde^2))
    dtilde <- c(rep(dtilde0, p - n), dtilde1)
  }
  u %*% (dtilde * t(u))
}

# long-only minimum variance from daily returns of the trailing 12 calendar months (formation month
# included), covariance matrix by analytical nonlinear shrinkage, same universe, cap and rebalancing as
# minvar_weights(). D: data.table(permno, date, mid, ret) keyed by mid. Stocks need a return on every
# trading day of the window.
minvar_weights_nl <- function(cm, D, months, wmax = 0.03, cores = 1L, bcol = "beta_M60", ucol = "large",
                              nmonths = 12L) {
  U <- cm[get(ucol) == TRUE & is.finite(get(bcol)), .(mid, permno)]
  f <- function(t) {
    ids <- U[mid == t, permno]
    x <- D[.(seq.int(t - nmonths + 1L, t)), on = "mid", nomatch = 0L][permno %in% ids]
    if (!nrow(x)) return(NULL)
    W <- dcast(x, date ~ permno, value.var = "ret")
    M <- as.matrix(W[, -1]); M <- M[, colSums(is.finite(M)) == nrow(M), drop = FALSE]
    if (ncol(M) < 30 || nrow(M) < 200) return(NULL)
    w <- gmv_long_only(nl_shrink(M), wmax)
    data.table(mid = t, permno = as.integer(colnames(M)), w = w, days = nrow(M))
  }
  res <- if (cores > 1) parallel::mclapply(months, f, mc.cores = cores) else lapply(months, f)
  out <- rbindlist(res)
  out[cm, ret_next := i.ret_next, on = .(permno, mid)]
  out
}

# cap individual weights at 'cap' and redistribute the excess pro rata to the uncapped positions
cap_weights <- function(w, cap) {
  if (length(w) * cap < 1) return(rep(1 / length(w), length(w)))
  w <- w / sum(w)
  repeat {
    over <- w > cap + 1e-12
    if (!any(over)) break
    excess <- sum(w[over] - cap); w[over] <- cap
    free <- w < cap - 1e-12
    w[free] <- w[free] + excess * w[free] / sum(w[free])
  }
  w
}

# predictive regression with a persistent predictor: OLS and the reduced-bias estimator of Amihud and
# Hurvich (2004, JFQA 39(4)), which corrects the Stambaugh (1999) small-sample bias. y[i] is the return
# over the period that starts at date i, x0[i] the predictor at date i and x1[i] the predictor at the end
# of that period (its innovation is contemporaneous with y[i]). Annual, non-overlapping observations.
pred_ah <- function(y, x0, x1) {
  n <- length(y); stopifnot(length(x0) == n, length(x1) == n)
  ar <- summary(lm(x1 ~ x0))$coefficients
  rho <- ar[2, 1]; se_rho <- ar[2, 2]
  rho_c <- min(rho + (1 + 3 * rho) / n + 3 * (1 + 3 * rho) / n^2, 0.9999)
  theta_c <- mean(x1) - rho_c * mean(x0)
  v_c <- x1 - theta_c - rho_c * x0
  ols <- summary(lm(y ~ x0))$coefficients
  aug <- summary(lm(y ~ x0 + v_c))$coefficients
  phi <- aug[3, 1]
  var_rho_c <- (1 + 3 / n + 9 / n^2)^2 * se_rho^2
  se_c <- sqrt(aug[2, 2]^2 + phi^2 * var_rho_c)
  u <- resid(lm(y ~ x0)); v <- resid(lm(x1 ~ x0))
  data.table(n = n, b_ols = ols[2, 1], t_ols = ols[2, 3], b_ah = aug[2, 1], t_ah = aug[2, 1] / se_c,
             rho = rho, rho_c = rho_c, corr_uv = cor(u, v),
             stambaugh_bias = -(cov(u, v) / var(v)) * (1 + 3 * rho) / n)
}

# ---- statistics ---------------------------------------------------------------------------
# Inference conventions
#  * tstat(): ordinary (iid) t-statistic of a mean. Used ONLY where GFKSW's own tables are
#    replicated (Table 1, Table 2 Panel A), so that the numbers are comparable.
#  * nw_t() / nw_reg(): Newey-West (1987) HAC t-statistics with Bartlett kernel. Default lag is
#    the Newey-West (1994) rule of thumb floor(4 (T/100)^(2/9)), i.e. 6 for T = 684 months and 4
#    for T = 144. Used for every own extension.
tstat <- function(x) { x <- x[is.finite(x)]; mean(x) / (sd(x) / sqrt(length(x))) }
nw_lag <- function(n) max(1L, as.integer(floor(4 * (n / 100)^(2 / 9))))
nw_t <- function(x, lag = NULL) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 8) return(NA_real_)
  if (is.null(lag)) lag <- nw_lag(n)
  fit <- lm(x ~ 1)
  se <- sqrt(NeweyWest(fit, lag = min(lag, n - 2L), prewhite = FALSE, adjust = TRUE)[1, 1])
  mean(x) / se
}
# regression with Newey-West t-statistics; returns coefficient, standard error and t for every term
nw_reg <- function(y, X, lag = NULL, w = NULL) {
  ok <- is.finite(y) & rowSums(!is.finite(X)) == 0
  if (!is.null(w)) ok <- ok & is.finite(w) & w > 0
  d <- data.frame(y = y[ok], X[ok, , drop = FALSE])
  if (is.null(lag)) lag <- nw_lag(sum(ok))
  fit <- if (is.null(w)) lm(y ~ ., data = d) else lm(y ~ ., data = d, weights = w[ok])
  se <- sqrt(diag(NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE)))
  data.table(term = c("alpha", colnames(X)), est = unname(coef(fit)), se = unname(se),
             t = unname(coef(fit) / se), n = sum(ok))
}

max_drawdown <- function(r) { w <- cumprod(1 + r); max(1 - w / cummax(c(1, w))[-1]) }

# Ledoit-Wolf (2008) HAC test for the difference of two Sharpe ratios (excess returns):
# delta method on (mean, uncentred second moment) of both series; long-run covariance by the
# prewhitened quadratic-spectral kernel with Andrews (1991) bandwidth (sandwich::lrvar returns
# the covariance of the sample means). Returns the annualised difference, its standard error and p.
sharpe_diff_test <- function(r1, r2, k = 12) {
  ok <- is.finite(r1) & is.finite(r2); r1 <- r1[ok]; r2 <- r2[ok]
  if (length(r1) < 24 || sd(r1) == 0 || sd(r2) == 0) return(c(diff = NA_real_, se = NA_real_, p = NA_real_))
  m1 <- mean(r1); m2 <- mean(r2); g1 <- mean(r1^2); g2 <- mean(r2^2)
  d <- m1 / sqrt(g1 - m1^2) - m2 / sqrt(g2 - m2^2)
  grad <- c(g1 / (g1 - m1^2)^1.5, -g2 / (g2 - m2^2)^1.5,
            -0.5 * m1 / (g1 - m1^2)^1.5, 0.5 * m2 / (g2 - m2^2)^1.5)
  Z <- cbind(r1, r2, r1^2, r2^2)
  V <- tryCatch(lrvar(Z, type = "Andrews", prewhite = TRUE, adjust = TRUE),     # prewhitening can fail
                error = function(e) lrvar(Z, type = "Andrews", prewhite = FALSE, adjust = TRUE))  # for sparse series
  se <- sqrt(drop(crossprod(grad, V %*% grad)))
  c(diff = sqrt(k) * d, se = sqrt(k) * se, p = 2 * pnorm(-abs(d / se)))
}

# Diebold-Mariano test on a loss differential series d (Newey-West variance)
dm_test <- function(d, lag) {
  d <- d[is.finite(d)]
  t <- nw_t(d, lag)
  c(mean = mean(d), t = t, p = 2 * pnorm(-abs(t)), n = length(d))
}

# cross-sectional trimming of a characteristic: values outside the [p, 1 - p] quantiles -> NA
trim_cs <- function(x, p = 0.005) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE, names = FALSE)
  fifelse(is.finite(x) & x >= q[1] & x <= q[2], x, NA_real_)
}

# ---- inference for GFKSW's state averages of overlapping-window statistics ---------------------
# Estimand: for state k and horizon h, the average over the formation months in state k of a
# statistic computed on the window of the next h return months (the annualised chained return or the
# window's CAPM / four-factor alpha). Main procedure: moving-block bootstrap (circular blocks of 'bl'
# months) of the monthly strategy and factor returns, jointly, with the state of every formation month
# held fixed; each window statistic and each state average is recomputed on every resampled history,
# so overlap, the clustering of states in calendar time and the estimator itself are reproduced
# exactly. The bootstrap standard error describes sampling variability under stationary returns that
# are unrelated to the state calendar. Indices are drawn once in the parent process (reproducible
# for any number of cores). sim_state_inference() below checks the size of the test by simulation.
mbb_index <- function(T, b) {
  st <- sample.int(T, ceiling(T / b), replace = TRUE)
  (unlist(lapply(st, function(s) (s - 1L + 0:(b - 1L)) %% T)) + 1L)[seq_len(T)]
}
# x: the period's return months in calendar order (strategy columns 'cols' and mkt, smb, hml, mom);
# state_by_pos[i]: state of the formation month whose window starts at return month i.
# idx: T x B matrix of bootstrap indices. Returns point estimates and bootstrap standard errors.
# vp = TRUE: volatility-preserving null model (v4 check). Every series (strategy and factors) is standardised with
# its own local volatility (local_sigma); the standardised series are resampled jointly with idx and rescaled with the
# local volatility of the original calendar position, so volatile periods stay where they occurred. With the
# identity index the original data are reproduced, so the point estimates are unchanged.
state_window_stats <- function(x, state_by_pos, hs, cols, idx = NULL, cores = 1L, vp = FALSE) {
  T <- nrow(x)
  Y <- as.matrix(x[, ..cols]); F0 <- as.matrix(x[, .(mkt, smb, hml, mom)]); F4 <- cbind(1, F0)
  if (vp) {
    stdz <- function(M) { mu <- colMeans(M); S <- apply(M, 2, local_sigma); list(mu = mu, S = S, U = sweep(M, 2, mu) / S) }
    sY <- stdz(Y); sF <- stdz(F0)
  }
  once <- function(ix) {
    out <- vector("list", length(cols) * length(hs)); r <- 0L
    for (cl in cols) {
      if (vp) {
        y <- sY$mu[[cl]] + sY$S[, cl] * sY$U[ix, cl]
        X4 <- cbind(1, sweep(sF$S * sF$U[ix, , drop = FALSE], 2, sF$mu, "+"))
      } else { y <- Y[ix, cl]; X4 <- F4[ix, , drop = FALSE] }
      X1 <- X4[, 1:2]; cs <- c(0, cumsum(log1p(y)))
      for (h in hs) {
        n <- T - h + 1L; s <- state_by_pos[seq_len(n)]
        g <- exp(cs[seq_len(n) + h] - cs[seq_len(n)]); a <- if (h == 1L) g - 1 else g^(12 / h) - 1
        ac <- a4 <- rep(NA_real_, n)
        if (h > 1L) for (i in seq_len(n)) {
          j <- i:(i + h - 1L)
          ac[i] <- .lm.fit(X1[j, , drop = FALSE], y[j])$coefficients[1]
          a4[i] <- .lm.fit(X4[j, , drop = FALSE], y[j])$coefficients[1]
        }
        m <- function(v) vapply(1:5, function(k) mean(v[s == k & !is.na(s)]), 0)
        r <- r + 1L
        out[[r]] <- data.table(strategy = cl, h = h, state = 1:5, ret = m(a), a_capm = m(ac), a_4f = m(a4))
      }
    }
    rbindlist(out)
  }
  est <- once(seq_len(T))
  if (is.null(idx)) return(est)
  bs <- parallel::mclapply(seq_len(ncol(idx)), function(b) once(idx[, b])[, b := b], mc.cores = cores)
  bs <- rbindlist(bs)
  se <- bs[, .(se_ret = sd(ret), se_capm = sd(a_capm), se_4f = sd(a_4f), B = .N), keyby = .(strategy, h, state)]
  out <- est[se, on = .(strategy, h, state)]
  # the central contrast: state 5 (low-beta stocks cheapest) minus state 1 (most expensive)
  d51 <- function(z, by) z[, .(d_ret = ret[state == 5] - ret[state == 1], d_capm = a_capm[state == 5] - a_capm[state == 1],
                               d_4f = a_4f[state == 5] - a_4f[state == 1]), by = by]
  dpt <- d51(est, c("strategy", "h"))
  dse <- d51(bs, c("strategy", "h", "b"))[, .(se_d_ret = sd(d_ret), se_d_capm = sd(d_capm), se_d_4f = sd(d_4f)), by = .(strategy, h)]
  setattr(out, "d51", dpt[dse, on = .(strategy, h)])
  out
}

# Size of three tests for a state average of annualised chained h-month returns, by simulation.
# Null: iid normal monthly returns (mean mu, sd sigma), unrelated to the states; the state calendar
# 'st' (one state per formation month, calendar order) is held fixed, so the clustering of states in
# a few episodes is reproduced. Tests: (a) GFKSW-style Newey-West with h lags on the state-filtered
# sequence; (b) calendar-index HAC (state_mean_hac, L = 2h); (c) the moving-block bootstrap used for
# the paper. Returns rejection rates of a true null at the nominal 5% level by horizon and state.
sim_state_inference <- function(st, hs = c(12L, 60L), nsim = 500L, B = 299L, bl = 24L, mu = 0.005,
                                sigma = 0.045, cores = 1L) {
  T <- length(st)
  wret <- function(y, h) { cs <- c(0, cumsum(log1p(y))); n <- length(y) - h + 1L; exp(cs[(1:n) + h] - cs[1:n])^(12 / h) - 1 }
  pop <- vapply(hs, function(h) mean(wret(rnorm(5e6, mu, sigma), h)), 0)   # population mean of the statistic
  seeds <- sample.int(.Machine$integer.max, nsim)                            # one seed per simulation
  one <- function(i) {
    set.seed(seeds[i])
    y <- rnorm(T, mu, sigma)
    BI <- replicate(B, mbb_index(T, bl))
    rbindlist(lapply(seq_along(hs), function(j) {
      h <- hs[j]; a <- wret(y, h); s <- st[seq_along(a)]; m0 <- pop[j]
      bs <- apply(BI, 2, function(ix) { ab <- wret(y[ix], h); vapply(1:5, function(k) mean(ab[s == k]), 0) })
      rbindlist(lapply(1:5, function(k) {
        z <- state_mean_hac(a - m0, s == k, 2L * h)
        data.table(h = h, state = k, rej_gfksw = abs(nw_t(a[s == k] - m0, h)) > 1.96, rej_hac = abs(z[["t"]]) > 1.96,
                   rej_boot = abs((mean(a[s == k]) - m0) / sd(bs[k, ])) > 1.96)
      }))
    }))
  }
  out <- rbindlist(parallel::mclapply(seq_len(nsim), one, mc.cores = cores))
  out[, lapply(.SD, mean), keyby = .(h, state)][, nsim := nsim]
}

# Calendar-index HAC alternative (kept for documentation; not used in the paper because it
# over-rejects in the simulation above when states cluster in a few episodes).
# Standard error for GFKSW's estimand: the average, over the formation months of one state, of a
# statistic computed on the overlapping window that starts after each formation month (a chained
# h-month return or a window alpha). With a_m the statistic for formation month m (all states, in
# calendar order, complete windows only) and d_m the state indicator, theta = sum(d a) / sum(d) and
# theta_hat - theta = sum(z) / N with z_m = d_m (a_m - theta). Var(sum z) is estimated on the full
# calendar index of formation months with a Bartlett kernel and L lags, so months of other states keep
# their calendar distance (GFKSW's Newey-West on the state-filtered sequence treats the last month of
# one episode and the first month of the next as adjacent). Windows overlap for h - 1 months, so L
# must exceed h; sim_state_inference() uses L = 2h.
state_mean_hac <- function(a, d, L) {
  ok <- is.finite(a); a <- a[ok]; d <- as.numeric(d[ok])
  N <- sum(d)
  if (N < 3) return(c(est = NA_real_, se = NA_real_, t = NA_real_, n = N))
  theta <- sum(d * a) / N
  z <- d * (a - theta); T <- length(z)
  s <- sum(z^2)
  for (l in seq_len(min(L, T - 1L))) s <- s + 2 * (1 - l / (L + 1)) * sum(z[(l + 1):T] * z[1:(T - l)])
  se <- sqrt(max(s, 0)) / N
  c(est = theta, se = se, t = theta / se, n = N)
}

# Calendar-time alternative to GFKSW's overlapping-window regressions (a different estimand).
# Every formation month m in 'starts' opens a window over return months m .. m + h - 1 (only
# complete windows inside x). Stacking all window-month observations and running one OLS is
# identical to a WLS regression on calendar months with weight c_t = number of open windows
# covering month t. The estimand is therefore a window-count-weighted average of monthly
# (abnormal) returns with factor loadings common to all windows of the state, not the average of
# window-specific alphas. Every calendar month enters once, so Newey-West standard errors apply
# (cf. Jegadeesh and Titman 1993). Months that no window covers have weight zero; they stay in the
# calendar with a zero score, so that a lag of 12 always spans 12 calendar months (wls_hac).
ct_reg <- function(x, starts, h, ycol, xcols = character(0), lag = 12L) {
  m <- x$mid
  stopifnot(all(diff(m) == 1L))                                     # complete monthly calendar
  i <- match(starts, m); ok <- !is.na(i) & i + h - 1L <= length(m)
  i <- i[ok]; ok2 <- (m[i + h - 1L] - m[i]) == (h - 1L); i <- i[ok2]
  if (length(i) < 3L) return(NULL)
  cnt <- tabulate(unlist(lapply(i, function(a) a:(a + h - 1L))), nbins = length(m))
  X <- cbind(alpha = rep(1, length(m)), if (length(xcols)) as.matrix(x[, ..xcols]))
  r <- wls_hac(x[[ycol]], X, cnt, lag)
  r[, ":="(windows = length(i), months = sum(cnt > 0))]
  r
}
# Weighted least squares with a Newey-West (Bartlett) covariance computed on the full calendar: rows
# with zero weight (or missing data) are excluded from the estimate but kept in the time ordering with a
# zero score, so autocovariances at lag l pair observations l calendar months apart. Without such rows
# the result equals sandwich::NeweyWest(lm(y ~ X - 1, weights = w), lag, prewhite = FALSE, adjust = TRUE).
wls_hac <- function(y, X, w, lag) {
  ok <- is.finite(y) & rowSums(!is.finite(X)) == 0 & is.finite(w) & w > 0
  Xo <- X[ok, , drop = FALSE]; wo <- w[ok]
  A <- crossprod(Xo * wo, Xo); b <- drop(solve(A, crossprod(Xo * wo, y[ok])))
  U <- matrix(0, length(y), ncol(X)); U[ok, ] <- (wo * drop(y[ok] - Xo %*% b)) * Xo
  M <- crossprod(U)
  for (l in seq_len(min(lag, nrow(U) - 1L))) {
    G <- crossprod(U[-seq_len(l), , drop = FALSE], U[seq_len(nrow(U) - l), , drop = FALSE])
    M <- M + (1 - l / (lag + 1)) * (G + t(G))
  }
  n <- sum(ok); Ai <- solve(A)
  se <- sqrt(diag(Ai %*% M %*% Ai) * n / (n - ncol(X)))
  data.table(term = colnames(X), est = unname(b), se = unname(se), t = unname(b / se), n = n)
}

# diversification of a weight table (mid, permno, w): averages over rebalancing dates.
# Positions are counted from a weight of 1 bp (the QP solver leaves numerically tiny weights).
div_stats <- function(W) {
  W <- W[is.finite(w) & w > 0]
  W[, w := w / sum(w), by = mid]
  W[, .(n = sum(w >= 1e-4), wmax = max(w), top10 = sum(head(sort(w, decreasing = TRUE), 10)), effn = 1 / sum(w^2)),
    by = mid][, .(n = mean(n), wmax = mean(wmax), top10 = mean(top10), effn = mean(effn))]
}

# performance summary of monthly simple returns r (rf = monthly risk-free, fac = factors)
perf_stats <- function(r, rf, fac, bench = NULL, k = 12, lag = NULL) {
  ex <- r - rf
  # (suppressWarnings: the market regressed on itself is an exact fit; harmless)
  capm <- suppressWarnings(nw_reg(ex, cbind(mkt = fac$mkt), lag))
  c4 <- suppressWarnings(nw_reg(ex, cbind(mkt = fac$mkt, smb = fac$smb, hml = fac$hml, mom = fac$mom), lag))
  out <- data.table(
    mean = k * mean(r), mean_ex = k * mean(ex), t_ex = nw_t(ex, lag), vol = sqrt(k) * sd(r),
    sharpe = sqrt(k) * mean(ex) / sd(ex),
    cagr = prod(1 + r)^(k / length(r)) - 1, mdd = max_drawdown(r),
    skew = mean((ex - mean(ex))^3) / sd(ex)^3, kurt = mean((ex - mean(ex))^4) / sd(ex)^4 - 3,
    beta = capm[term == "mkt", est], alpha_capm = k * capm[term == "alpha", est],
    t_capm = capm[term == "alpha", t], alpha_4f = k * c4[term == "alpha", est],
    t_4f = c4[term == "alpha", t], b_hml = c4[term == "hml", est], b_mom = c4[term == "mom", est],
    months = length(r))
  if (!is.null(bench)) {
    a <- r - bench
    out[, ":="(active = k * mean(a), t_active = nw_t(a, lag), te = sqrt(k) * sd(a), ir = sqrt(k) * mean(a) / sd(a))]
  }
  out
}

# =============================================================================
# v4 additions: daily returns from CRSP's daily file, dependence diagnostics and inference checks
# =============================================================================

# ---- daily return reconstruction ----------------------------------------------------------------
# dsh: dsf_v2 records dated inside held stock-months (permno, date, mid, ret, del = delisting-day record);
# dnx: delisting-day records dated in the month after a held month (permno, mid = held month, ret_nd).
# The daily path of every holding runs over the trading calendar 'cal' (date, mid) of the holding month:
#  (a) daily records are used on their own dates, including delisting-day records inside the month;
#  (b) a delisting-day record dated in the first days of the following month, which CRSP's monthly file credits
#      to the held month, is booked on the held month's last trading day with a record;
#  (c) where the monthly file excludes an in-month delisting-day return ('conflict'), "monthly" (baseline)
#      replaces that return so that the stock-month compounds to the monthly return, "daily" keeps it;
#  (d) any other difference is 'unexplained' and booked on the last trading day with a record; differences of at
#      most 'tol' (rounding of daily returns) are booked there as well, so every stock-month compounds exactly.
# Days without a record earn zero: proceeds after a delisting are held as cash until the month-end
# rebalancing, as in CRSP's monthly return. Returns the day-level panel g and the stock-month table sm.
daily_panel_dsf <- function(h, cal, dsh, dnx, conflict = "monthly", tol = 1e-4) {
  hh <- h[w > 0, .(mid = mid + 1L, permno, w, rm = ret_next)]
  g <- cal[hh, on = "mid", allow.cartesian = TRUE, nomatch = 0L]
  rec <- dsh[unique(hh[, .(permno, mid)]), on = .(permno, mid), nomatch = 0L]
  cc <- copy(cal)[, cd := date]
  rec <- cc[rec, on = .(mid, date), roll = -Inf]                 # next trading day on or after the record date
  last_day <- cal[, .(last = max(date)), by = mid]
  rec[is.na(cd), cd := last_day[.SD, on = "mid", x.last]]
  moved <- rec[cd != date, .N]
  rec <- rec[, .(ret = if (all(is.na(ret))) NA_real_ else prod(1 + ret[is.finite(ret)]) - 1, del = any(del)), by = .(permno, date = cd)]
  g[rec, ":="(ret = i.ret, del = i.del), on = .(permno, date)]
  g[, rec := is.finite(ret)][rec == FALSE, ret := 0][is.na(del), del := FALSE]
  setorder(g, mid, permno, date)
  g[, k := seq_len(.N), by = .(mid, permno)]
  sm <- g[, .(w = w[1], rm = rm[1], r_in = prod(1 + ret) - 1, r_ex = prod(1 + ret[!del]) - 1, has_d = any(del),
              k_d = if (any(del)) max(k[del]) else NA_integer_, k_last = if (any(rec)) max(k[rec]) else .N,
              nrec = sum(rec), ndays = .N), by = .(mid, permno)]
  sm[dnx, r_nd := i.ret_nd, on = .(permno, mid)]
  sm[, dev := (1 + rm) / (1 + r_in) - 1]
  sm[, cls := fifelse(is.finite(dev) & abs(dev) <= tol, "match",
              fifelse(!has_d & is.finite(r_nd) & abs((1 + r_in) * (1 + r_nd) / (1 + rm) - 1) <= tol, "next_month_delisting",
              fifelse(has_d, "conflict", "unexplained")))]
  # bookings on the last trading day with a record: rounding differences of at most 'tol' ('match'), next-month
  # delisting-day returns and unexplained differences; conflicts: the delisting-day return itself is replaced
  adj <- rbind(
    sm[cls %in% c("match", "next_month_delisting", "unexplained") & 1 + r_in > 0 & dev != 0,
       .(mid, permno, k = k_last, f = (1 + rm) / (1 + r_in), set = NA_real_)],
    if (conflict == "monthly") sm[cls == "conflict", .(mid, permno, k = k_d, f = NA_real_, set = (1 + rm) / (1 + r_ex) - 1)])
  if (nrow(adj)) {
    g[adj, ":="(f = i.f, set = i.set), on = .(mid, permno, k)]
    g[is.finite(set), ret := set][is.finite(f), ret := (1 + ret) * f - 1]
    g[, c("f", "set") := NULL]
  }
  attr(sm, "moved") <- moved
  list(g = g, sm = sm)
}
# v3 reconstruction (crsp_daily, which omits delisting-day records): the whole difference between the monthly
# return and the compounded daily returns is booked on the last trading day with a record (book = TRUE), or
# ignored (book = FALSE: delisting returns left out of the daily path).
daily_panel_cd <- function(h, cal, dret, book = TRUE) {
  g <- cal[h[w > 0, .(mid = mid + 1L, permno, w, rm = ret_next)], on = "mid", allow.cartesian = TRUE, nomatch = 0L]
  g[dret, ret := i.ret, on = .(permno, date)]
  g[, rec := is.finite(ret)][rec == FALSE, ret := 0]
  setorder(g, mid, permno, date)
  g[, ":="(rd = prod(1 + ret) - 1, last_rec = if (any(rec)) max(which(rec)) else .N, k = seq_len(.N)), by = .(mid, permno)]
  if (book) g[k == last_rec & (1 + rd) > 0, ret := (1 + ret) * (1 + rm) / (1 + rd) - 1]
  g
}
# daily value of a monthly rebalanced portfolio: weights drift with the daily returns within the month
daily_port_from_panel <- function(g) {
  g[, gr := cumprod(1 + ret), by = .(mid, permno)]
  v <- g[, .(v = sum(w * gr)), by = .(mid, date)]
  v0 <- g[k == 1L, .(v0 = sum(w)), by = mid]
  setorder(v, date); v[v0, v0 := i.v0, on = "mid"]
  v[, ret := v / shift(v, fill = v0[1]) - 1, by = mid]
  v[, .(date, mid, ret)]
}
# Lo and MacKinlay (1988) variance ratio of q-period log returns with the heteroskedasticity-robust z-statistic
vr_test <- function(r, q) {
  x <- log1p(r[is.finite(r)]); n <- length(x); mu <- mean(x); e <- x - mu
  s_a <- sum(e^2) / (n - 1)
  sq <- frollsum(x, q)[q:n] - q * mu
  vr <- (sum(sq^2) / (q * (n - q + 1) * (1 - q / n))) / s_a
  den <- sum(e^2)^2
  theta <- sum(vapply(seq_len(q - 1L), function(j) (2 * (q - j) / q)^2 * n * sum(e[(j + 1):n]^2 * e[1:(n - j)]^2) / den, 0))
  c(vr = vr, z = sqrt(n) * (vr - 1) / sqrt(theta))
}

# ---- dependence diagnostics ------------------------------------------------------------------------
# Politis and White (2004) automatic block length for the circular block bootstrap of the mean, with the
# correction of Patton, Politis and White (2009)
pw_block <- function(x) {
  n <- length(x); x <- x - mean(x)
  KN <- max(5, ceiling(sqrt(log10(n)))); mmax <- ceiling(sqrt(n)) + KN
  rho <- acf(x, lag.max = mmax, plot = FALSE)$acf[-1]
  crit <- 2 * sqrt(log10(n) / n); mhat <- mmax
  for (m in 1:(mmax - KN)) if (all(abs(rho[m:(m + KN - 1)]) < crit)) { mhat <- m; break }
  M <- min(2 * mhat, mmax)
  R <- acf(x, lag.max = M, type = "covariance", plot = FALSE)$acf[, 1, 1]
  k <- 0:M; lam <- ifelse(abs(k / M) <= 0.5, 1, 2 * (1 - abs(k / M)))
  G <- 2 * sum(lam[-1] * k[-1] * R[-1]); g0 <- R[1] + 2 * sum(lam[-1] * R[-1])
  (2 * G^2 / ((4 / 3) * g0^2))^(1 / 3) * n^(1 / 3)
}
# simple (non-robust) variance ratio of overlapping q-month sums
vr_simple <- function(v, q) { s <- frollsum(v, q)[q:length(v)]; var(s) / (q * var(v)) }

# ---- block bootstraps of state averages of chained h-month returns (returns only) -----------------
win_ret <- function(y, h) {
  cs <- c(0, cumsum(log1p(y))); n <- length(y) - h + 1L
  g <- exp(cs[seq_len(n) + h] - cs[seq_len(n)]); if (h == 1L) g - 1 else g^(12 / h) - 1
}
state_means <- function(a, s) vapply(1:5, function(k) { z <- a[s == k]; if (length(z)) mean(z) else NA_real_ }, 0)
# (i) null model, as in state_window_stats(): returns resampled in circular blocks, state calendar fixed (5 x B)
boot_null <- function(y, st, h, idx) { n <- length(y) - h + 1L; s <- st[seq_len(n)]
  apply(idx, 2, function(ix) state_means(win_ret(y[ix], h), s)) }
# (ii) observed relation: circular blocks of formation months, each with its own state and window statistic
boot_pairs <- function(a, s, idx) apply(idx, 2, function(ix) state_means(a[ix], s[ix]))
# (iii) volatility-preserving null model: standardised returns are resampled in blocks and rescaled with the
# local volatility of the original calendar position, so that volatile periods stay where they occurred.
# Local volatility: two-sided moving average of squared demeaned returns over 2k + 1 months.
local_sigma <- function(y, k = 12L) { e2 <- (y - mean(y))^2; n <- length(y)
  sqrt(vapply(seq_len(n), function(t) mean(e2[max(1L, t - k):min(n, t + k)]), 0)) }
boot_vp <- function(y, st, h, idx, sig = local_sigma(y)) { n <- length(y) - h + 1L; s <- st[seq_len(n)]; m <- mean(y); u <- (y - m) / sig
  apply(idx, 2, function(ix) state_means(win_ret(m + sig * u[ix], h), s)) }
# summary of a bootstrap matrix M (5 x B) for the state means and the state-5 minus state-1 contrast:
# standard error, t-statistic, Monte Carlo standard error of the standard error (delta method with the
# kurtosis of the draws), and the share of centred draws at least as extreme as the estimate (null models)
boot_summary <- function(th, M) {
  ok <- colSums(!is.finite(M)) == 0; M <- M[, ok, drop = FALSE]; B <- ncol(M)
  one <- function(est, v) { se <- sd(v); k <- mean((v - mean(v))^4) / se^4
    data.table(est = est, se = se, t = est / se, se_mc = se * sqrt(max(k - 1, 0) / (4 * B)),
               p_null = mean(abs(v - mean(v)) >= abs(est)), lo = quantile(v, 0.025, names = FALSE), hi = quantile(v, 0.975, names = FALSE), B = B) }
  rbind(one(th[1], M[1, ])[, stat := "state1"], one(th[5], M[5, ])[, stat := "state5"], one(th[5] - th[1], M[5, ] - M[1, ])[, stat := "d51"])
}

# ---- AR(1)-GARCH(1,1) with Student-t innovations (maximum likelihood) for the size simulations ---------
garch_fit <- function(y) {
  n <- length(y)
  filt <- function(p) { mu <- p[1]; phi <- p[2]; om <- exp(p[3]); al <- plogis(p[4]); be <- plogis(p[5]) * (1 - al)
    e <- y - mu - phi * (c(mean(y), y[-n]) - mu); h <- numeric(n); h[1] <- var(y)
    for (t in 2:n) h[t] <- om + al * e[t - 1]^2 + be * h[t - 1]
    list(e = e, h = h, mu = mu, phi = phi, om = om, al = al, be = be, nu = 2.05 + exp(p[6])) }
  nll <- function(p) { f <- filt(p); sc <- sqrt(f$h * (f$nu - 2) / f$nu); -sum(dt(f$e / sc, f$nu, log = TRUE) - log(sc)) }
  o <- optim(c(mean(y), 0, log(var(y) * 0.1), qlogis(0.1), qlogis(0.85), log(6)), nll, method = "BFGS", control = list(maxit = 2000))
  f <- filt(o$par)
  list(mu = f$mu, phi = f$phi, omega = f$om, alpha = f$al, beta = f$be, nu = f$nu, sigma = sqrt(f$h), conv = o$convergence)
}
rt_std <- function(n, nu) rt(n, nu) * sqrt((nu - 2) / nu)
garch_sim <- function(g, n, burn = 200L) {
  N <- n + burn; z <- rt_std(N, g$nu); h <- numeric(N); y <- numeric(N); e <- numeric(N)
  h[1] <- g$omega / max(1e-6, 1 - g$alpha - g$beta); y[1] <- g$mu
  for (t in 2:N) { h[t] <- g$omega + g$alpha * e[t - 1]^2 + g$beta * h[t - 1]; e[t] <- sqrt(h[t]) * z[t]
                   y[t] <- g$mu + g$phi * (y[t - 1] - g$mu) + e[t] }
  y[(burn + 1):N]
}
# Size of the tests for state averages and for the state-5 minus state-1 contrast of chained h-month returns,
# by simulation under a null in which the state calendar 'st' carries no information about expected returns.
# gen(): one simulated monthly return series of length length(st). Rejection rates (nominal 5%, two-sided) are
# computed around the Monte Carlo mean of each estimator, i.e. its expectation under the DGP, which need not be
# zero for compounded returns when volatility differs across windows. crit95: size-adjusted critical value
# (95% quantile of |t| under the DGP).
sim_tests <- function(st, gen, h = 12L, nsim = 1000L, B = 299L, cores = 1L) {
  T <- length(st); n <- T - h + 1L; s <- st[seq_len(n)]
  seeds <- sample.int(.Machine$integer.max, nsim)
  one <- function(i) {
    set.seed(seeds[i]); y <- gen()
    a <- win_ret(y, h); th <- state_means(a, s)
    I1 <- replicate(B, mbb_index(T, 24L)); I3 <- replicate(B, mbb_index(T, 60L)); I2 <- replicate(B, mbb_index(n, 24L))
    se <- function(M) { M <- M[, colSums(!is.finite(M)) == 0, drop = FALSE]; c(apply(M, 1, sd), sd(M[5, ] - M[1, ])) }
    list(th = th, gf = vapply(1:5, function(k) { z <- a[s == k]; mean(z) / nw_t(z, h) }, 0),
         null24 = se(boot_null(y, st, h, I1)), null60 = se(boot_null(y, st, h, I3)),
         vp24 = se(boot_vp(y, st, h, I1)), pairs24 = se(boot_pairs(a, s, I2)))
  }
  R <- parallel::mclapply(seq_len(nsim), one, mc.cores = cores)
  TH <- do.call(rbind, lapply(R, `[[`, "th")); D <- TH[, 5] - TH[, 1]
  dev <- cbind(abs(sweep(TH, 2, colMeans(TH))), abs(D - mean(D)))
  out <- list(data.table(test = "gfksw_nw", stat = paste0("state", 1:5),
                         rej = colMeans(dev[, 1:5] / do.call(rbind, lapply(R, `[[`, "gf")) > 1.96),
                         crit95 = apply(dev[, 1:5] / do.call(rbind, lapply(R, `[[`, "gf")), 2, quantile, 0.95)))
  for (nm in c("null24", "null60", "vp24", "pairs24")) {
    tt <- dev / do.call(rbind, lapply(R, `[[`, nm))
    out[[length(out) + 1L]] <- data.table(test = nm, stat = c(paste0("state", 1:5), "d51"), rej = colMeans(tt > 1.96),
                                          crit95 = apply(tt, 2, quantile, 0.95))
  }
  out <- rbindlist(out)
  out[, ":="(nsim = nsim, B = B, mc_se = sqrt(rej * (1 - rej) / nsim), e_d51 = mean(D), sd_d51 = sd(D))]
  out
}
