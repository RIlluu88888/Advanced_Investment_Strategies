# =============================================================================
# Revision 3: applies two corrections to results/res.rds without rerunning the bootstraps and size
# simulations of build.R (build.R and R/functions.R contain the same corrections for future fresh builds):
#  (1) Table 7 (res$freq_beta): monthly market excess return compounded as prod(1 + mkt + rf) - prod(1 + rf)
#      instead of prod(1 + mkt) - 1 (build.R, section 7).
#  (2) Appendix Table 15 (res$T34_ct) and its real-time variant (res$T2_rt_ct): Newey-West covariance of the
#      calendar-time WLS on the full monthly calendar, months without weight entering with a zero score
#      (R/functions.R: ct_reg, wls_hac). Point estimates are unchanged.
# Inputs: results/res.rds and cache/fd.rds (daily factors). Before any change, the stored values are
# reproduced from the same inputs with the previous formulas (stops otherwise).
# Usage (in lowvol_v4): Rscript update_r3.R > logs/update_r3.log 2>&1
# =============================================================================
source("R/functions.R")
in_period <- function(m, p) m >= p[1] & m <= p[2]                    # as in build.R
res <- readRDS("results/res.rds"); cfg <- res$cfg
in_md5 <- unname(tools::md5sum("results/res.rds"))
stopifnot(is.null(res$meta$update_r3))                                # apply once
logmsg("input results/res.rds", in_md5)

# ---- (1) Table 7: monthly market betas -------------------------------------------------------------
beta_se <- function(r, rf, mkt) { z <- suppressWarnings(nw_reg(r - rf, cbind(mkt = mkt))); z[term == "mkt", .(beta = est, se = se)] }
fd <- readRDS("cache/fd.rds")
fm <- fd[, .(mid = year(date) * 12L + month(date) - 1L, mkt, rf)][in_period(mid, cfg$full),
         .(rf = prod(1 + rf) - 1, mkt_old = prod(1 + mkt) - 1, mkt_new = prod(1 + mkt + rf) - prod(1 + rf),
           mret = prod(1 + mkt + rf) - 1), keyby = mid]
S <- copy(res$series); setkey(S, mid)
lo <- function(s) res$LO[strategy == s & in_period(mid, cfg$full), .(strategy, mid, ret)]
# monthly compounded daily paths (build.R verifies these equal the monthly engine to < 1e-10 for the
# long-only portfolios and L-H; L-H and BN use the compounded daily T-bill return, the market its daily factor)
DM <- rbind(lo("LowBeta"), lo("Min-variance"), lo("Large-cap VW"), lo("Large-cap EW"),
            S[in_period(mid, cfg$full), .(strategy = "HighBeta", mid, ret = ret5)],
            fm[S[in_period(mid, cfg$full), .(mid, ret1, ret5)], on = "mid"][
              , .(strategy = c("Low - high", "Beta neutral (fixed 25%)"), ret = c(rf + ret1 - ret5, ret1 - cfg$bn_short * ret5 + cfg$bn_short * rf)),
              by = mid],
            fm[, .(strategy = "Market (FF)", mid, ret = mret)], use.names = TRUE)
DM <- fm[DM, on = "mid"]
old <- DM[, beta_se(ret, rf, mkt_old), keyby = strategy][res$freq_beta, on = "strategy"]
dev_old <- old[, max(abs(c(beta - beta_monthly, se - se_monthly)))]
logmsg("Table 7: stored monthly betas reproduced with the previous formula, max abs difference", signif(dev_old, 3))
stopifnot(nrow(old) == 8L, dev_old < 1e-8)
new <- DM[, beta_se(ret, rf, mkt_new), keyby = strategy]
stopifnot(abs(new[strategy == "Market (FF)", beta] - 1) < 1e-10)       # market on itself
fb <- copy(res$freq_beta)
fb[new, ":="(beta_monthly = i.beta, se_monthly = i.se), on = "strategy"]
print(merge(res$freq_beta[, .(strategy, beta_old = beta_monthly, se_old = se_monthly)], fb[, .(strategy, beta_new = beta_monthly, se_new = se_monthly)], by = "strategy"))

# ---- (2) calendar-time regressions ------------------------------------------------------------------
# previous version of ct_reg (Newey-West on the months with positive weight only), used only to verify the inputs
ct_reg_prev <- function(x, starts, h, ycol, xcols = character(0), lag = 12L) {
  m <- x$mid
  i <- match(starts, m); ok <- !is.na(i) & i + h - 1L <= length(m)
  i <- i[ok]; ok2 <- (m[i + h - 1L] - m[i]) == (h - 1L); i <- i[ok2]
  if (length(i) < 3L) return(NULL)
  cnt <- tabulate(unlist(lapply(i, function(a) a:(a + h - 1L))), nbins = length(m))
  X <- if (length(xcols)) as.matrix(x[, ..xcols]) else matrix(nrow = length(m), ncol = 0)
  y <- x[[ycol]]
  if (length(xcols)) {
    r <- nw_reg(y, X, lag = lag, w = cnt)
  } else {
    sel <- cnt > 0 & is.finite(y)
    fit <- lm(y[sel] ~ 1, weights = cnt[sel])
    se <- sqrt(NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE)[1, 1])
    r <- data.table(term = "alpha", est = unname(coef(fit)), se = se, t = unname(coef(fit)) / se, n = sum(sel))
  }
  r[, ":="(windows = length(i), months = sum(cnt > 0))]
  r
}
# same as ct_table() in build.R, with the regression function as an argument
ct_table_f <- function(S, st, p, f, cols = c("LH", "BN"), hs = c(12L, 36L, 60L)) {
  x <- S[in_period(mid, p)]; setkey(x, mid)
  rbindlist(lapply(cols, function(cl) rbindlist(lapply(hs, function(h) rbindlist(lapply(1:5, function(s) {
    starts <- st[state == s, mid]
    a <- f(x, starts, h, cl, character(0), cfg$ct_lag)
    if (is.null(a)) return(NULL)
    b <- f(x, starts, h, cl, "mkt", cfg$ct_lag)
    d <- f(x, starts, h, cl, c("mkt", "smb", "hml", "mom"), cfg$ct_lag)
    data.table(strategy = cl, h = h, state = s, mean = a$est, se_mean = a$se, t_mean = a$t,
               a_capm = b[term == "alpha", est], t_capm = b[term == "alpha", t],
               a_4f = d[term == "alpha", est], t_4f = d[term == "alpha", t],
               months = a$months, windows = a$windows)
  }))))))
}
# real-time states exactly as in build.R
rt <- S[in_period(mid, cfg$full) & is.finite(spread), .(mid, fmid, spread)]
rt[, state := NA_integer_]
for (i in seq_len(nrow(rt))) if (i > cfg$rt_min) {
  br <- quantile(rt$spread[1:(i - 1)], seq(0.2, 0.8, 0.2), names = FALSE)
  rt$state[i] <- findInterval(rt$spread[i], br) + 1L
}
rt <- rt[!is.na(state)]
stopifnot(min(rt$mid) == res$rt_start)
st_rep <- res$states_rep
prev_rep <- ct_table_f(S, st_rep, cfg$rep, ct_reg_prev); prev_rt <- ct_table_f(S, rt, cfg$full, ct_reg_prev, hs = 12L)
stopifnot(isTRUE(all.equal(prev_rep, res$T34_ct, tolerance = 1e-10, check.attributes = FALSE)),
          isTRUE(all.equal(prev_rt, res$T2_rt_ct, tolerance = 1e-10, check.attributes = FALSE)))
logmsg("calendar-time tables: stored values reproduced with the previous HAC")
new_rep <- ct_table_f(S, st_rep, cfg$rep, ct_reg); new_rt <- ct_table_f(S, rt, cfg$full, ct_reg, hs = 12L)
est_cols <- c("strategy", "h", "state", "mean", "a_capm", "a_4f", "months", "windows")
stopifnot(isTRUE(all.equal(new_rep[, ..est_cols], prev_rep[, ..est_cols], tolerance = 1e-12)),
          isTRUE(all.equal(new_rt[, ..est_cols], prev_rt[, ..est_cols], tolerance = 1e-12)))
# without gaps the new HAC equals the previous one: check on a fully covered calendar (all weights 1)
x0 <- S[in_period(mid, cfg$rep)]; X0 <- cbind(alpha = 1, as.matrix(x0[, .(mkt, smb, hml, mom)]))
a0 <- wls_hac(x0$LH, X0, rep(1, nrow(x0)), cfg$ct_lag); b0 <- nw_reg(x0$LH, X0[, -1], lag = cfg$ct_lag, w = rep(1, nrow(x0)))
stopifnot(max(abs(a0$se - b0$se)) < 1e-12)
cmp <- prev_rep[new_rep, on = .(strategy, h, state)][, .(strategy, h, state, months,
         t_mean_old = t_mean, t_mean_new = i.t_mean, t_capm_old = t_capm, t_capm_new = i.t_capm, t_4f_old = t_4f, t_4f_new = i.t_4f)]
print(cmp, digits = 3)

# ---- store ----------------------------------------------------------------------------------------
res$freq_beta <- fb; res$T34_ct <- new_rep; res$T2_rt_ct <- new_rt
res$meta$update_r3 <- list(time = Sys.time(), input_md5 = in_md5,
                           functions_md5 = unname(tools::md5sum("R/functions.R")), build_md5 = unname(tools::md5sum("build.R")),
                           changed = c("freq_beta (monthly columns)", "T34_ct (se, t)", "T2_rt_ct (se, t)"),
                           check_table7_prev = dev_old)
saveRDS(res, "results/res.rds")
logmsg("results/res.rds updated", unname(tools::md5sum("results/res.rds")))
