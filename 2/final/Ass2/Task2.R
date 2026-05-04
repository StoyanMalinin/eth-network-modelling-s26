#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

inv_logit <- function(x) 1 / (1 + exp(-x))

read_adjacency <- function(path) {
  m <- as.matrix(read.csv(path, header = FALSE))
  storage.mode(m) <- "integer"
  diag(m) <- 0L
  m
}

read_matrix <- function(path) {
  m <- as.matrix(read.csv(path, header = FALSE))
  storage.mode(m) <- "double"
  m
}

coef_row <- function(model) {
  sm <- summary(model)
  cf <- sm$coefficients
  cn <- colnames(cf)
  if (is.null(cn)) cn <- rep("", ncol(cf))

  pick <- function(primary, fallbacks = character(0)) {
    idx <- match(primary, cn)
    if (!is.na(idx)) return(idx)
    for (fb in fallbacks) {
      idx <- match(fb, cn)
      if (!is.na(idx)) return(idx)
    }
    NA_integer_
  }

  i_est <- pick("Estimate")
  i_se <- pick("Std. Error", c("Std.Error", "Std. err."))
  i_mcmc <- pick("MCMC %", c("MCMC%"))
  i_z <- pick("z value", c("z", "z-value"))
  i_p <- pick("Pr(>|z|)", c("Pr(>|t|)", "Pr(>|W|)"))

  data.frame(
    term = rownames(cf),
    estimate = unname(cf[, i_est]),
    std_error = unname(cf[, i_se]),
    mcmc_percent = if (is.na(i_mcmc)) NA_real_ else unname(cf[, i_mcmc]),
    z = if (is.na(i_z)) NA_real_ else unname(cf[, i_z]),
    p_value = if (is.na(i_p)) NA_real_ else unname(cf[, i_p]),
    row.names = NULL
  )
}

library(coda)
suppressPackageStartupMessages(library(ergm))
suppressPackageStartupMessages(library(network))

dir.create("task2_figures", showWarnings = FALSE)
dir.create("task2_results", showWarnings = FALSE)
refit <- identical(Sys.getenv("TASK2_REFIT", "0"), "1")

adj <- read_adjacency(file.path("Glasgow", "f1.csv"))
n <- nrow(adj)

demo <- read.csv(file.path("Glasgow", "demographic.csv"))
gender <- as.integer(demo$gender)

alcohol_df <- read.csv(file.path("Glasgow", "alcohol.csv"))
alcohol <- as.numeric(alcohol_df$Wave.1)

logdistance <- read_matrix(file.path("Glasgow", "logdistance.csv"))
diag(logdistance) <- 0

net <- network(adj, directed = TRUE, matrix.type = "adjacency")
set.vertex.attribute(net, "gender", gender)
set.vertex.attribute(net, "alcohol", alcohol)
set.vertex.attribute(net, "age", demo$age)

set.seed(1)
m1 <- ergm(net ~ edges + nodematch("gender"), estimate = "MPLE")
m1_tab <- coef_row(m1)
write.csv(m1_tab, file = file.path("task2_results", "part1_model_coefficients.csv"), row.names = FALSE)

b_edges <- unname(coef(m1)[["edges"]])
nm_name <- grep("^nodematch\\.gender", names(coef(m1)), value = TRUE)
b_same <- unname(coef(m1)[[nm_name]])

p_diff <- inv_logit(b_edges)
p_same <- inv_logit(b_edges + b_same)
or_same <- exp(b_same)

write.csv(
  data.frame(
    p_tie_diff_gender = p_diff,
    p_tie_same_gender = p_same,
    odds_ratio_same_vs_diff = or_same
  ),
  file = file.path("task2_results", "part1_tie_probabilities.csv"),
  row.names = FALSE
)


alpha_gw_esp <- 0.2
nsim_gof <- 20

set.seed(2)
  formula_full <-
    net ~ edges +
      nodematch("gender") +
      mutual +
      gwesp(alpha_gw_esp, fixed = TRUE, type = "OTP") +
      nodeofactor("gender", levels = 1) +
      nodeifactor("gender", levels = 1) +
      nodeocov("alcohol") +
      nodeicov("alcohol") +
      absdiff("alcohol") +
      absdiff("age") +
      edgecov(logdistance) +
	    gwidegree(0.6, fixed = TRUE) + 
	    gwodegree(0.6, fixed = TRUE)

m2_rds <- file.path("task2_results", "full_model.rds")
if (!refit && file.exists(m2_rds)) {
  m2 <- readRDS(m2_rds)
} else {
  m2_noisy <- ergm(formula_full, 
                   control = control.ergm(
                     main.method = "Stochastic-Approximation",
                     MCMLE.maxit = 60,
                     MCMC.burnin = 1e6,
                     MCMC.interval = 1e4,
                     MCMC.samplesize = 6e3,
                     seed = 2
                   ))
  write.csv(
    coef_row(m2_noisy),
    file = file.path("task2_results", "full_model_noisy_coefficients.csv"),
    row.names = FALSE
  )

  ctrl <- control.ergm(
    init = coef(m2_noisy),
    MCMLE.maxit = 60,
    MCMC.burnin = 1e6,
    MCMC.interval = 1e4,
    MCMC.samplesize = 6e3,
    seed = 2
  )
  
  m2 <- try(ergm(formula_full, control = ctrl), silent = TRUE)
  if (inherits(m2, "try-error")) {
    writeLines(as.character(m2), con = file.path("task2_results", "full_model_mcmle_error.txt"))
    stop("try-error")
  }
  saveRDS(m2, file = m2_rds)
}

m2_tab <- coef_row(m2)
write.csv(m2_tab, file = file.path("task2_results", "full_model_coefficients.csv"), row.names = FALSE)

coef_hist <- m2$coef.hist
tail_window <- min(6L, nrow(coef_hist))
tail_hist <- coef_hist[(nrow(coef_hist) - tail_window + 1L):nrow(coef_hist), , drop = FALSE]
max_abs_delta_last <- if (nrow(tail_hist) >= 2L) {
  max(abs(diff(tail_hist)))
} else {
  NA_real_
}

write.csv(
  data.frame(
    mcmle_iterations = m2$iterations,
    mcmle_maxit = m2$control$MCMLE.maxit,
    mcmle_hit_maxit = m2$iterations >= m2$control$MCMLE.maxit,
    mcmle_failure = m2$failure,
    mcmle_termination = m2$control$MCMLE.termination,
    mcmle_confidence = m2$control$MCMLE.confidence,
    mcmle_conv_min_pval = m2$control$MCMLE.conv.min.pval,
    max_abs_score = max(abs(m2$gradient)),
    max_abs_coef_change_last_window = max_abs_delta_last,
    last_step_length = tail(m2$steplen.hist, 1)
  ),
  file = file.path("task2_results", "full_model_convergence_summary.csv"),
  row.names = FALSE
)

# MCMC sample diagnostics for the last estimation round (heuristic).
# Note: m2$sample contains centered statistics (simulated - observed), so means should be near 0.
if (requireNamespace("coda", quietly = TRUE) && inherits(m2$sample, "mcmc.list") && length(m2$sample) >= 1L) {
  s <- m2$sample[[1]]
  sm <- summary(s)$statistics
  terms <- rownames(sm)
  mean_centered <- as.numeric(sm[, "Mean"])
  ts_se <- as.numeric(sm[, "Time-series SE"])
  z_mean0 <- mean_centered / ts_se
  p_mean0 <- 2 * stats::pnorm(-abs(z_mean0))

  ess <- coda::effectiveSize(s)
  acf1 <- coda::autocorr.diag(s, lags = 1)[1, ]
  geweke_z <- coda::geweke.diag(s, frac1 = 0.1, frac2 = 0.5)$z

  diag_df <- data.frame(
    term = terms,
    mean_centered = mean_centered,
    ts_se = ts_se,
    z_mean0 = z_mean0,
    p_mean0 = p_mean0,
    ess = ess[terms],
    acf_lag1 = acf1[terms],
    geweke_z = geweke_z[terms],
    row.names = NULL
  )
  write.csv(
    diag_df,
    file = file.path("task2_results", "mcmc_sample_diagnostics.csv"),
    row.names = FALSE
  )
}

grDevices::pdf(file.path("task2_figures", "mcmc_diagnostics_full_model.pdf"), width = 11, height = 8.5)
mcmc.diagnostics(m2)
auto_data <- autocorr(m2$sample)
plot.new()
text(x = 0, y = 1,
     labels = paste(capture.output(print(auto_data[[1]])), collapse = "\n"),
     adj = c(0, 1),
     family = "mono",
     cex = 0.6)
grDevices::dev.off()

set.seed(3)
gof_rds <- file.path("task2_results", "gof_full_model.rds")
if (!refit && file.exists(gof_rds)) {
  g <- readRDS(gof_rds)
} else {
  g <- gof(m2, GOF = ~ odegree + idegree + distance + triadcensus)
  saveRDS(g, file = gof_rds)
}

grDevices::png(file.path("task2_figures", "gof_full_model.png"), width = 1400, height = 900)
plot(g)
grDevices::dev.off()

gof_outside <- function(pobs, bds, stat_type) {
  lower <- bds[1, ]
  upper <- bds[2, ]
  outside <- (pobs < lower) | (pobs > upper)
  if (!any(outside)) return(NULL)
  data.frame(
    stat_type = stat_type,
    bin = names(pobs)[outside],
    pobs = unname(pobs[outside]),
    lower = unname(lower[outside]),
    upper = unname(upper[outside]),
    direction = ifelse(pobs[outside] < lower[outside], "below", "above"),
    row.names = NULL
  )
}

gof_details <- rbind(
  gof_outside(g$pobs.odeg, g$bds.odeg, "odegree"),
  gof_outside(g$pobs.ideg, g$bds.ideg, "idegree"),
  gof_outside(g$pobs.dist, g$bds.dist, "distance"),
  gof_outside(g$pobs.triadcensus, g$bds.triadcensus, "triadcensus")
)

if (!is.null(gof_details)) {
  write.csv(
    gof_details,
    file = file.path("task2_results", "gof_outside_envelope_details.csv"),
    row.names = FALSE
  )

  gof_summary <- aggregate(
    bin ~ stat_type + direction,
    data = gof_details,
    FUN = length
  )
  names(gof_summary)[names(gof_summary) == "bin"] <- "n_bins"
  gof_summary <- gof_summary[order(gof_summary$stat_type, gof_summary$direction), ]

  write.csv(
    gof_summary,
    file = file.path("task2_results", "gof_outside_envelope_summary.csv"),
    row.names = FALSE
  )
}

cat("Part (1): p(tie | diff gender) =", signif(p_diff, 4), "\n")
cat("Part (1): p(tie | same gender) =", signif(p_same, 4), "\n")
cat("Part (1): OR(same vs diff) =", signif(or_same, 4), "\n\n")

cat("Full model coefficients written to task2_results/full_model_coefficients.csv\n")
cat("GOF figure written to task2_figures/gof_full_model.png\n")
cat("MCMC diagnostics written to task2_figures/mcmc_diagnostics_full_model.pdf\n")
