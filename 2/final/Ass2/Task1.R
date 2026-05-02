#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

delta_gwodegree <- function(k, alpha) {
  base <- -expm1(-alpha)
  base^k
}

read_adjacency <- function(path) {
  m <- as.matrix(read.csv(path, header = FALSE))
  storage.mode(m) <- "integer"
  if (nrow(m) != ncol(m)) stop("Adjacency matrix must be square.")
  diag(m) <- 0L
  m
}

observed_stats <- function(adj) {
  n <- nrow(adj)
  outdeg <- rowSums(adj != 0)
  m <- sum(adj != 0)
  max_edges <- n * (n - 1)
  list(
    n = n,
    outdeg = outdeg,
    max_outdeg = max(outdeg),
    edges = m,
    max_edges = max_edges,
    density = m / max_edges
  )
}

plot_change_statistics <- function(k_values, alphas, outfile, main) {
  cols <- grDevices::hcl.colors(length(alphas), "Dark 3")
  y <- vapply(
    alphas,
    function(a) delta_gwodegree(k_values, a),
    numeric(length(k_values))
  )

  grDevices::png(outfile, width = 900, height = 600)
  on.exit(grDevices::dev.off(), add = TRUE)
  matplot(
    x = k_values,
    y = y,
    type = "l",
    lwd = 2,
    lty = 1,
    col = cols,
    xlab = "Current out-degree k (sender)",
    ylab = expression(delta[GWD](k, alpha)),
    main = main
  )
  legend(
    "topright",
    legend = paste0("alpha = ", format(alphas, digits = 3)),
    col = cols,
    lty = 1,
    lwd = 2,
    cex = 0.9,
    bty = "n"
  )
}

simulate_formula_outdegree_counts <- function(
  formula,
  coef,
  samples = 100L,
  seed = 1L
) {
  if (!requireNamespace("ergm", quietly = TRUE)) {
    stop("Package 'ergm' is required. Please install it to run Task 1 simulations.")
  }
  if (!requireNamespace("network", quietly = TRUE)) {
    stop("Package 'network' is required. Please install it to run Task 1 simulations.")
  }

  set.seed(seed)
  sims <- ergm::simulate_formula(
    formula,
    nsim = samples,
    coef = coef,
    output = "network",
    verbose = FALSE
  )

  n <- network::network.size(sims[[1]])
  counts <- matrix(0L, nrow = samples, ncol = n) # bins 0..(n-1)
  edges <- integer(samples)
  for (s in seq_len(samples)) {
    m <- as.matrix(sims[[s]], matrix.type = "adjacency")
    diag(m) <- 0
    a <- m != 0
    outdeg <- rowSums(a)
    counts[s, ] <- tabulate(outdeg + 1L, nbins = n)
    edges[s] <- sum(a)
  }

  max_edges <- n * (n - 1)
  list(counts = counts, edges = edges, density = edges / max_edges)
}

estimate_theta_edges_match_density_ergm <- function(
  net,
  theta_gwd,
  alpha,
  target_density,
  seed = 1L,
  sims_per_eval = 30L
) {
  formula <- net ~ edges + gwodegree(alpha, fixed = TRUE)

  f <- function(theta_edges) {
    sim <- simulate_formula_outdegree_counts(
      formula = formula,
      coef = c(theta_edges, theta_gwd),
      samples = sims_per_eval,
      seed = seed
    )
    mean(sim$density) - target_density
  }

  lo <- -12
  hi <- 0
  flo <- f(lo)
  fhi <- f(hi)
  tries <- 0L
  while (sign(flo) == sign(fhi) && tries < 6L) {
    lo <- lo - 4
    hi <- hi + 4
    flo <- f(lo)
    fhi <- f(hi)
    tries <- tries + 1L
  }
  if (sign(flo) == sign(fhi)) {
    stop("Failed to bracket root for theta_edges; try increasing the search interval.")
  }

  uniroot(f, lower = lo, upper = hi, tol = 1e-3)$root
}

boxplot_compare_counts <- function(
  counts_ref,
  counts_comp,
  degrees,
  label_ref,
  label_comp,
  outfile,
  title,
  observed_counts = NULL
) {
  d_idx <- degrees + 1L
  ref_mat <- counts_ref[, d_idx, drop = FALSE]
  comp_mat <- counts_comp[, d_idx, drop = FALSE]

  data_list <- vector("list", length(degrees) * 2L)
  at <- numeric(length(degrees) * 2L)
  cols <- character(length(degrees) * 2L)
  for (i in seq_along(degrees)) {
    data_list[[2L * i - 1L]] <- ref_mat[, i]
    data_list[[2L * i]] <- comp_mat[, i]
    at[2L * i - 1L] <- i - 0.18
    at[2L * i] <- i + 0.18
    cols[2L * i - 1L] <- "grey75"
    cols[2L * i] <- "steelblue3"
  }

  grDevices::png(outfile, width = 1100, height = 600)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- par(mar = c(5, 5, 4, 2) + 0.1)
  on.exit(par(op), add = TRUE)

  boxplot(
    data_list,
    at = at,
    xaxt = "n",
    col = cols,
    border = "grey30",
    outline = FALSE,
    ylab = "Number of nodes",
    main = title
  )
  axis(1, at = seq_along(degrees), labels = degrees)
  legend(
    "topright",
    legend = c(label_ref, label_comp, "Observed"),
    fill = c("grey75", "steelblue3", NA),
    border = c("grey30", "grey30", NA),
    pch = c(NA, NA, 16),
    col = c(NA, NA, "firebrick"),
    bty = "n"
  )

  if (!is.null(observed_counts)) {
    points(seq_along(degrees), observed_counts[d_idx], pch = 16, col = "firebrick")
  }
}

main <- function() {
  dir.create("task1_figures", showWarnings = FALSE)
  dir.create("task1_results", showWarnings = FALSE)

  adj <- read_adjacency(file.path("Glasgow", "f1.csv"))
  obs <- observed_stats(adj)

  write.csv(
    data.frame(node = seq_len(obs$n), outdegree = obs$outdeg),
    file = file.path("task1_results", "observed_outdegree_wave1.csv"),
    row.names = FALSE
  )

  k_values <- 0:obs$max_outdeg

  plot_change_statistics(
    k_values,
    alphas = -1,
    outfile = file.path("task1_figures", "gwodegree_change_alpha_-1.png"),
    main = "Change statistic for gwodegree (alpha = -1)"
  )
  plot_change_statistics(
    k_values,
    alphas = -0.5,
    outfile = file.path("task1_figures", "gwodegree_change_alpha_-0.5.png"),
    main = "Change statistic for gwodegree (alpha = -0.5)"
  )

  alphas_nonneg <- c(0, 0.1, 0.2, 0.3, 0.5, log(2), 1, 2, 3)
  plot_change_statistics(
    k_values,
    alphas = alphas_nonneg,
    outfile = file.path("task1_figures", "gwodegree_change_alpha_nonneg.png"),
    main = "Change statistic for gwodegree (non-negative alpha)"
  )

  p_obs <- obs$density
  theta_edges_m1 <- stats::qlogis(p_obs)

  alpha_fixed <- 0.3
  theta_gwd_vals <- c(-2, 2)

  if (!requireNamespace("network", quietly = TRUE)) stop("Package 'network' is required.")
  net <- network::network(adj, directed = TRUE, matrix.type = "adjacency")

  theta_edges_m2 <- vapply(
    seq_along(theta_gwd_vals),
    function(idx) {
      estimate_theta_edges_match_density_ergm(
        net = net,
        theta_gwd = theta_gwd_vals[idx],
        alpha = alpha_fixed,
        target_density = p_obs,
        seed = 10L + idx,
        sims_per_eval = 25L
      )
    },
    numeric(1)
  )

  params <- data.frame(
    model = c("M1_edges_only", "M2_gwd_-2", "M2_gwd_+2"),
    theta_edges = c(theta_edges_m1, theta_edges_m2[1], theta_edges_m2[2]),
    theta_gwodegree = c(0, theta_gwd_vals[1], theta_gwd_vals[2]),
    alpha = c(NA, alpha_fixed, alpha_fixed),
    target_density = c(p_obs, p_obs, p_obs)
  )
  write.csv(params, file = file.path("task1_results", "model_parameters.csv"), row.names = FALSE)

  sims_m1 <- simulate_formula_outdegree_counts(
    formula = net ~ edges,
    coef = c(theta_edges_m1),
    samples = 100L,
    seed = 2L
  )
  sims_m2_neg <- simulate_formula_outdegree_counts(
    formula = net ~ edges + gwodegree(alpha_fixed, fixed = TRUE),
    coef = c(theta_edges_m2[1], -2),
    samples = 100L,
    seed = 3L
  )
  sims_m2_pos <- simulate_formula_outdegree_counts(
    formula = net ~ edges + gwodegree(alpha_fixed, fixed = TRUE),
    coef = c(theta_edges_m2[2], 2),
    samples = 100L,
    seed = 4L
  )

  write.csv(sims_m1$counts, file = file.path("task1_results", "outdegree_counts_M1.csv"), row.names = FALSE)
  write.csv(sims_m2_neg$counts, file = file.path("task1_results", "outdegree_counts_M2_gwd_-2.csv"), row.names = FALSE)
  write.csv(sims_m2_pos$counts, file = file.path("task1_results", "outdegree_counts_M2_gwd_+2.csv"), row.names = FALSE)

  nonzero_cols <- function(mat) which(colSums(mat) > 0)
  kmax <- max(
    max(nonzero_cols(sims_m1$counts)),
    max(nonzero_cols(sims_m2_neg$counts)),
    max(nonzero_cols(sims_m2_pos$counts))
  ) - 1L
  kmax_plot <- min(kmax, 25L)
  degrees <- 0:kmax_plot

  obs_counts <- tabulate(obs$outdeg + 1L, nbins = obs$n)

  boxplot_compare_counts(
    counts_ref = sims_m1$counts,
    counts_comp = sims_m2_neg$counts,
    degrees = degrees,
    label_ref = "M1 (edges only)",
    label_comp = expression(paste("M2 (", theta[GWD], " = -2, ", alpha, " = 0.3)")),
    outfile = file.path("task1_figures", "outdegree_box_M2_gwd_-2.png"),
    title = "Out-degree distribution: M2 (theta_GWD = -2) vs M1 (matched density)",
    observed_counts = obs_counts
  )

  boxplot_compare_counts(
    counts_ref = sims_m1$counts,
    counts_comp = sims_m2_pos$counts,
    degrees = degrees,
    label_ref = "M1 (edges only)",
    label_comp = expression(paste("M2 (", theta[GWD], " = +2, ", alpha, " = 0.3)")),
    outfile = file.path("task1_figures", "outdegree_box_M2_gwd_+2.png"),
    title = "Out-degree distribution: M2 (theta_GWD = +2) vs M1 (matched density)",
    observed_counts = obs_counts
  )

  cat("Observed density:", signif(p_obs, 4), "\n")
  cat("M1 theta_edges:", signif(theta_edges_m1, 4), "\n")
  cat("M2 theta_edges (theta_gwd=-2, alpha=0.3):", signif(theta_edges_m2[1], 4), "\n")
  cat("M2 theta_edges (theta_gwd=+2, alpha=0.3):", signif(theta_edges_m2[2], 4), "\n")
  cat("Plots written to task1_figures/\n")
  cat("Tables written to task1_results/\n")
}

if (sys.nframe() == 0L) {
  main()
}
