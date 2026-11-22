# Bayesian Weibull regression for right-censored survival data,
# fitted with a random-walk Metropolis sampler written from scratch.
#
# Data: survival::lung (NCCTG advanced lung cancer trial, n = 227 complete cases)
# Model: Weibull accelerated failure time (AFT)
#     log(T_i) = x_i' beta + sigma * W_i,   W_i ~ standard minimum extreme value
# Priors: beta_j ~ N(0, 10^2),  log(sigma) ~ N(0, 2^2)
#
# The sampler is checked against the maximum likelihood fit from survival::survreg.
#
# Run from the repository root:  Rscript R/analysis.R

suppressPackageStartupMessages({
  library(survival)
  library(coda)
  library(ggplot2)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------
d <- survival::lung
d <- d[complete.cases(d[, c("time", "status", "age", "sex", "ph.ecog")]), ]
d$event   <- as.integer(d$status == 2)          # 1 = death, 0 = censored
d$age_c   <- (d$age - 60) / 10                  # age centred at 60, per 10 years
d$female  <- as.integer(d$sex == 2)
d$ecog    <- factor(pmin(d$ph.ecog, 2), levels = 0:2,
                    labels = c("0", "1", "2+"))  # ECOG 2 and 3 grouped
d$time_m  <- d$time / 30.44                     # months

X <- model.matrix(~ age_c + female + ecog, data = d)
y <- log(d$time_m)
delta <- d$event
p <- ncol(X)
param_names <- c(colnames(X), "log_sigma")

cat("n =", nrow(d), " deaths =", sum(delta), "\n")

# ---------------------------------------------------------------------------
# 2. Log posterior
# ---------------------------------------------------------------------------
log_lik <- function(theta) {
  beta  <- theta[1:p]
  sigma <- exp(theta[p + 1])
  z <- (y - X %*% beta) / sigma
  # event:    log f(log t) = -log(sigma) + z - exp(z)
  # censored: log S(log t) = -exp(z)
  sum(delta * (-log(sigma) + z) - exp(z))
}
log_prior <- function(theta) {
  sum(dnorm(theta[1:p], 0, 10, log = TRUE)) + dnorm(theta[p + 1], 0, 2, log = TRUE)
}
log_post <- function(theta) log_lik(theta) + log_prior(theta)

# ---------------------------------------------------------------------------
# 3. Random-walk Metropolis
#    Proposal covariance = (2.38^2 / d) * inverse Hessian at the posterior mode
# ---------------------------------------------------------------------------
mode_fit <- optim(rep(0, p + 1), function(t) -log_post(t), method = "BFGS", hessian = TRUE)
prop_cov <- (2.38^2 / (p + 1)) * solve(mode_fit$hessian)
L <- t(chol(prop_cov))

run_chain <- function(n_iter, init) {
  draws <- matrix(NA_real_, n_iter, p + 1, dimnames = list(NULL, param_names))
  cur <- init; cur_lp <- log_post(cur); acc <- 0
  for (i in seq_len(n_iter)) {
    prop <- cur + as.vector(L %*% rnorm(p + 1))
    prop_lp <- log_post(prop)
    if (log(runif(1)) < prop_lp - cur_lp) { cur <- prop; cur_lp <- prop_lp; acc <- acc + 1 }
    draws[i, ] <- cur
  }
  attr(draws, "acceptance") <- acc / n_iter
  draws
}

n_chains <- 4; n_iter <- 20000; burn <- 5000
chains <- lapply(seq_len(n_chains), function(k) {
  init <- mode_fit$par + rnorm(p + 1, 0, 0.5)   # dispersed starting values
  run_chain(n_iter, init)
})
acc_rates <- sapply(chains, attr, "acceptance")
cat("Acceptance rates:", round(acc_rates, 3), "\n")

mcmc_list <- coda::mcmc.list(lapply(chains, function(ch) coda::mcmc(ch[(burn + 1):n_iter, ])))

# ---------------------------------------------------------------------------
# 4. Convergence diagnostics
# ---------------------------------------------------------------------------
rhat <- coda::gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[, 1]
ess  <- coda::effectiveSize(mcmc_list)
diag_tab <- data.frame(parameter = param_names, Rhat = round(rhat, 3), ESS = round(ess))
write.csv(diag_tab, "output/tables/convergence.csv", row.names = FALSE)
print(diag_tab)

# Trace plots
post <- do.call(rbind, lapply(seq_along(chains), function(k)
  data.frame(iter = (burn + 1):n_iter, chain = factor(k), chains[[k]][(burn + 1):n_iter, ])))
names(post)[3:ncol(post)] <- param_names
pretty <- c("(Intercept)" = "Intercept", age_c = "Age (per 10 y)", female = "Female",
            ecog1 = "ECOG 1", "ecog2+" = "ECOG 2+", log_sigma = "log(sigma)")
trace_long <- do.call(rbind, lapply(param_names, function(v)
  data.frame(iter = post$iter, chain = post$chain, value = post[[v]], parameter = unname(pretty[v]))))
thin <- trace_long[trace_long$iter %% 10 == 0, ]
p_trace <- ggplot(thin, aes(iter, value, colour = chain)) +
  geom_line(linewidth = 0.2, alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
  labs(x = "Iteration (after burn-in)", y = NULL, colour = "Chain",
       title = "Trace plots: four chains from dispersed starting values") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("output/figures/trace_plots.png", p_trace, width = 8, height = 7, dpi = 150, bg = "white")

# ---------------------------------------------------------------------------
# 5. Posterior summaries vs maximum likelihood (survreg)
#    Report time ratios exp(beta): > 1 means longer survival
# ---------------------------------------------------------------------------
all_draws <- as.matrix(mcmc_list)
mle <- survreg(Surv(time_m, event) ~ age_c + female + ecog, data = d, dist = "weibull")
mle_est <- c(coef(mle), log(mle$scale))
mle_se  <- sqrt(diag(vcov(mle)))

summ <- data.frame(
  parameter   = pretty[param_names],
  post_mean   = colMeans(all_draws),
  post_sd     = apply(all_draws, 2, sd),
  q2.5        = apply(all_draws, 2, quantile, 0.025),
  q97.5       = apply(all_draws, 2, quantile, 0.975),
  mle         = mle_est,
  mle_se      = mle_se,
  row.names = NULL
)
summ_round <- summ; summ_round[, -1] <- round(summ_round[, -1], 3)
write.csv(summ_round, "output/tables/posterior_vs_mle.csv", row.names = FALSE)
print(summ_round)

tr_idx <- 2:p
time_ratio <- data.frame(
  covariate = pretty[param_names[tr_idx]],
  TR        = round(exp(colMeans(all_draws[, tr_idx])), 2),
  lower95   = round(exp(apply(all_draws[, tr_idx], 2, quantile, 0.025)), 2),
  upper95   = round(exp(apply(all_draws[, tr_idx], 2, quantile, 0.975)), 2),
  P_TR_gt_1 = round(colMeans(all_draws[, tr_idx] > 0), 3),
  row.names = NULL
)
write.csv(time_ratio, "output/tables/time_ratios.csv", row.names = FALSE)
print(time_ratio)

# Posterior densities with MLE overlay
dens_long <- do.call(rbind, lapply(param_names, function(v)
  data.frame(value = all_draws[, v], parameter = unname(pretty[v]))))
mle_df <- data.frame(parameter = pretty[param_names], mle = mle_est)
p_dens <- ggplot(dens_long, aes(value)) +
  geom_density(fill = "#1b6ca8", alpha = 0.35, colour = "#1b6ca8") +
  geom_vline(data = mle_df, aes(xintercept = mle), colour = "#c0392b", linetype = "dashed") +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  labs(x = NULL, y = "Posterior density",
       title = "Posterior distributions (blue) and maximum likelihood estimates (red)") +
  theme_minimal(base_size = 11)
ggsave("output/figures/posterior_vs_mle.png", p_dens, width = 9, height = 5, dpi = 150, bg = "white")

# ---------------------------------------------------------------------------
# 6. Posterior predictive survival curves vs Kaplan-Meier, by ECOG
#    (for a 60-year-old man)
# ---------------------------------------------------------------------------
tgrid <- seq(0.1, 33, length.out = 200)
ecog_levels <- c("0", "1", "2+")
curve_df <- do.call(rbind, lapply(ecog_levels, function(e) {
  x0 <- c(1, 0, 0, as.numeric(e == "1"), as.numeric(e == "2+"))
  mu <- all_draws[, 1:p] %*% x0
  sig <- exp(all_draws[, p + 1])
  S <- sapply(tgrid, function(t) exp(-exp((log(t) - mu) / sig)))
  data.frame(time = tgrid, ecog = e, S = colMeans(S),
             lo = apply(S, 2, quantile, 0.025), hi = apply(S, 2, quantile, 0.975))
}))
km_fit <- survfit(Surv(time_m, event) ~ ecog, data = d[d$female == 0, ])
km_df <- data.frame(time = km_fit$time, S = km_fit$surv,
                    ecog = rep(sub("ecog=", "", names(km_fit$strata)), km_fit$strata))
p_surv <- ggplot() +
  geom_ribbon(data = curve_df, aes(time, ymin = lo, ymax = hi, fill = ecog), alpha = 0.2) +
  geom_line(data = curve_df, aes(time, S, colour = ecog), linewidth = 0.9) +
  geom_step(data = km_df, aes(time, S, colour = ecog), linetype = "dashed", linewidth = 0.4) +
  scale_colour_manual(values = c("#16a085", "#1b6ca8", "#c0392b")) +
  scale_fill_manual(values = c("#16a085", "#1b6ca8", "#c0392b")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Months", y = "Survival probability", colour = "ECOG", fill = "ECOG",
       title = "Posterior survival curves with 95% credible bands",
       subtitle = "Solid: Bayesian Weibull model for a 60-year-old man. Dashed: Kaplan-Meier, men only.") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("output/figures/posterior_survival.png", p_surv, width = 8, height = 5.5, dpi = 150, bg = "white")

# Posterior median survival time (months) by ECOG for a 60-year-old man and woman
med_tab <- do.call(rbind, lapply(c(0, 1), function(f) do.call(rbind, lapply(ecog_levels, function(e) {
  x0 <- c(1, 0, f, as.numeric(e == "1"), as.numeric(e == "2+"))
  mu <- all_draws[, 1:p] %*% x0; sig <- exp(all_draws[, p + 1])
  med <- exp(mu + sig * log(log(2)))
  data.frame(sex = ifelse(f == 1, "Female", "Male"), ecog = e,
             median_months = round(median(med), 1),
             lower95 = round(quantile(med, 0.025), 1), upper95 = round(quantile(med, 0.975), 1),
             row.names = NULL)
}))))
write.csv(med_tab, "output/tables/median_survival.csv", row.names = FALSE)
print(med_tab)

writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
cat("\nDone. Results written to output/\n")
