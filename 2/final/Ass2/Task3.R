#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

read_adjacency_bin <- function(path) {
  m <- as.matrix(read.csv(path, header = FALSE))
  storage.mode(m) <- "integer"
  diag(m) <- 0L
  m != 0L
}

read_matrix <- function(path) {
  m <- as.matrix(read.csv(path, header = FALSE))
  storage.mode(m) <- "double"
  diag(m) <- 0
  m
}

jaccard_edges <- function(a, b) {
  inter <- sum(a & b)
  uni <- sum(a | b)
  if (uni == 0) return(NA_real_)
  inter / uni
}

`%||%` <- function(x, y) if (is.null(x)) y else x

extract_estimates_table <- function(ans) {
  effects_df <- as.data.frame(ans$requestedEffects)

  tbl <- data.frame(
    dependent = effects_df$name,
    type = effects_df$type,
    effect = effects_df$effectName,
    short_name = effects_df$shortName,
    interaction1 = effects_df$interaction1,
    estimate = ans$theta,
    se = ans$se,
    tconv = ans$tconv,
    stringsAsFactors = FALSE
  )

  if (!any(tbl$type == "rate") && !is.null(ans$rate) && length(ans$rate) > 0L) {
    rate_se <- if (!is.null(ans$vrate) && length(ans$vrate) == length(ans$rate)) {
      ans$vrate
    } else {
      rep(NA_real_, length(ans$rate))
    }

    rate_tbl <- data.frame(
      dependent = unique(tbl$dependent)[1] %||% NA_character_,
      type = "rate",
      effect = paste0("Rate parameter period ", seq_along(ans$rate)),
      short_name = "Rate",
      interaction1 = NA_character_,
      estimate = as.numeric(ans$rate),
      se = as.numeric(rate_se),
      tconv = NA_real_,
      stringsAsFactors = FALSE
    )
    tbl <- rbind(rate_tbl, tbl)
  }

  tbl
}

write_estimates_csv <- function(ans, path) {
  tbl <- extract_estimates_table(ans)
  write.csv(tbl, file = path, row.names = FALSE)
  invisible(tbl)
}

extract_gof_joint_p <- function(gof_obj) {
  p <- tryCatch(gof_obj$Joint$p, error = function(e) NA_real_)
  if (is.null(p) || length(p) != 1L) return(NA_real_)
  as.numeric(p)
}

write_gof_pvalues_csv <- function(gof_list, path) {
  df <- data.frame(
    stat = names(gof_list),
    p_value = vapply(gof_list, extract_gof_joint_p, numeric(1)),
    stringsAsFactors = FALSE
  )
  write.csv(df, file = path, row.names = FALSE)
  invisible(df)
}

save_gof_pngs <- function(gof_list, prefix, overwrite = FALSE) {
  dir.create("task3_figures", showWarnings = FALSE)
  for (stat in names(gof_list)) {
    out_path <- file.path("task3_figures", paste0(prefix, "_", stat, ".png"))
    if (!overwrite && file.exists(out_path)) next
    grDevices::png(out_path, width = 1600, height = 1100, res = 200)
    print(plot(gof_list[[stat]]))
    grDevices::dev.off()
  }
  invisible(TRUE)
}

geodesic_dist_undirected <- function(adj, max_dist = 8L) {
  n <- nrow(adj)
  neigh <- lapply(seq_len(n), function(v) which(adj[v, ]))

  counts <- integer(max_dist + 2L)
  names(counts) <- c(paste0("d", seq_len(max_dist)), paste0("ge", max_dist + 1L), "inf")

  for (s in seq_len(n - 1L)) {
    dist <- rep.int(-1L, n)
    dist[s] <- 0L
    queue <- integer(n)
    head <- 1L
    tail <- 1L
    queue[tail] <- s

    while (head <= tail) {
      v <- queue[head]
      head <- head + 1L
      nd <- dist[v] + 1L
      for (u in neigh[[v]]) {
        if (dist[u] < 0L) {
          dist[u] <- nd
          tail <- tail + 1L
          queue[tail] <- u
        }
      }
    }

    for (t in (s + 1L):n) {
      d <- dist[t]
      if (d < 0L) {
        counts[["inf"]] <- counts[["inf"]] + 1L
      } else if (d <= max_dist) {
        counts[[d]] <- counts[[d]] + 1L
      } else {
        counts[[paste0("ge", max_dist + 1L)]] <- counts[[paste0("ge", max_dist + 1L)]] + 1L
      }
    }
  }

  counts
}

GeodesicDistributionUndirected <- function(
  i,
  obsData,
  sims,
  period,
  groupName,
  varName,
  max_dist = 8L
) {
  dep <- obsData[[groupName]]$depvars[[varName]]
  if (!(attr(dep, "missing") | attr(dep, "structural"))) {
    x <- RSiena:::sparseMatrixExtraction0(i, obsData, sims, period, groupName, varName)
  } else {
    x <- RSiena:::sparseMatrixExtraction(i, obsData, sims, period, groupName, varName)
  }
  a <- as.matrix(x != 0)
  diag(a) <- FALSE
  und <- a | t(a)
  geodesic_dist_undirected(und, max_dist = max_dist)
}

main <- function() {
  dir.create("task3_figures", showWarnings = FALSE)
  dir.create("task3_results", showWarnings = FALSE)
  refit <- identical(Sys.getenv("TASK3_REFIT", "0"), "1")

  a1 <- read_adjacency_bin(file.path("Glasgow", "f1.csv"))
  a2 <- read_adjacency_bin(file.path("Glasgow", "f2.csv"))
  a3 <- read_adjacency_bin(file.path("Glasgow", "f3.csv"))
  n <- nrow(a1)

  j12 <- jaccard_edges(a1, a2)
  j23 <- jaccard_edges(a2, a3)
  j13 <- jaccard_edges(a1, a3)

  jacc_df <- data.frame(
    pair = c("wave1-wave2", "wave2-wave3", "wave1-wave3"),
    jaccard = c(j12, j23, j13),
    edges_t1 = c(sum(a1), sum(a2), sum(a1)),
    edges_t2 = c(sum(a2), sum(a3), sum(a3)),
    n = n
  )
  write.csv(jacc_df, file = file.path("task3_results", "jaccard_indices.csv"), row.names = FALSE)

  cat("Jaccard indices written to task3_results/jaccard_indices.csv\n")
  cat("J(w1,w2) =", signif(j12, 4), "\n")
  cat("J(w2,w3) =", signif(j23, 4), "\n")
  cat("J(w1,w3) =", signif(j13, 4), "\n\n")

  suppressPackageStartupMessages(library(RSiena))

  demo <- read.csv(file.path("Glasgow", "demographic.csv"))
  gender <- as.integer(demo$gender)

  alcohol_df <- read.csv(file.path("Glasgow", "alcohol.csv"))
  alcohol_mat <- as.matrix(alcohol_df)
  storage.mode(alcohol_mat) <- "integer"

  logdist <- read_matrix(file.path("Glasgow", "logdistance.csv"))

  net_arr <- array(0L, dim = c(n, n, 3))
  net_arr[, , 1] <- a1 * 1L
  net_arr[, , 2] <- a2 * 1L
  net_arr[, , 3] <- a3 * 1L

  friendship <- sienaDependent(net_arr, type = "oneMode")

  gender_cov <- coCovar(gender)
  logdist_cov <- coDyadCovar(logdist)
  alcohol_cov <- varCovar(alcohol_mat)

  gof_iterations <- 200
  dat1 <- sienaDataCreate(friendship, gender_cov, logdist_cov, alcohol_cov)
  eff1 <- getEffects(dat1)

  eff1 <- includeEffects(eff1, recip, transTrip, cycle3, inPop, name = "friendship")
  eff1 <- includeEffects(eff1, sameX, egoX, altX, interaction1 = "gender_cov", name = "friendship")
  eff1 <- includeEffects(eff1, simX, interaction1 = "alcohol_cov", name = "friendship")
  eff1 <- includeEffects(eff1, X, interaction1 = "logdist_cov", name = "friendship")

  alg1 <- sienaAlgorithmCreate(
    projname = file.path("task3_results", "saom_part1_network_only"),
    n3 = 3000,
    seed = 101
  )

  ans1_path <- file.path("task3_results", "ans_part1_network_only.rds")
  if (!refit && file.exists(ans1_path)) {
    ans1 <- readRDS(ans1_path)
  } else {
    ans1 <- siena07(alg1, data = dat1, effects = eff1, batch = TRUE, returnDeps = TRUE)
    saveRDS(ans1, file = ans1_path)
  }

  sink(file.path("task3_results", "summary_part1_network_only.txt"))
  print(summary(ans1))
  sink()

  write_estimates_csv(ans1, file.path("task3_results", "estimates_part1_network_only.csv"))

  gof1_path <- file.path("task3_results", "gof_part1_network_only.rds")
  gof1_pdf <- file.path("task3_figures", "gof_part1_network_only.pdf")
  if (!refit && file.exists(gof1_path) && file.exists(gof1_pdf)) {
    gof1 <- readRDS(gof1_path)
  } else {
    gof1_out <- sienaGOF(ans1, OutdegreeDistribution, varName = "friendship", iterations = gof_iterations, verbose = FALSE)
    gof1_in <- sienaGOF(ans1, IndegreeDistribution, varName = "friendship", iterations = gof_iterations, verbose = FALSE)
    gof1_tri <- sienaGOF(ans1, TriadCensus, varName = "friendship", iterations = gof_iterations, verbose = FALSE)
    gof1_geo <- sienaGOF(ans1, GeodesicDistributionUndirected, varName = "friendship", iterations = gof_iterations, verbose = FALSE)
    gof1 <- list(outdegree = gof1_out, indegree = gof1_in, geodesic = gof1_geo, triad = gof1_tri)
    saveRDS(gof1, file = gof1_path)

    grDevices::pdf(gof1_pdf, width = 11, height = 8.5)
    plot(gof1_out)
    plot(gof1_in)
    plot(gof1_tri)
    plot(gof1_geo)
    grDevices::dev.off()
  }

  write_gof_pvalues_csv(gof1, file.path("task3_results", "gof_part1_pvalues.csv"))
  save_gof_pngs(gof1, prefix = "gof_part1_network_only", overwrite = refit)

  alcohol <- sienaDependent(alcohol_mat, type = "behavior")
  dat2 <- sienaDataCreate(friendship, alcohol, gender_cov, logdist_cov)

  eff2 <- getEffects(dat2)
  eff2 <- includeEffects(eff2, recip, transTrip, cycle3, inPop, name = "friendship")
  eff2 <- includeEffects(eff2, sameX, egoX, altX, interaction1 = "gender_cov", name = "friendship")
  eff2 <- includeEffects(eff2, simX, interaction1 = "alcohol", name = "friendship")
  eff2 <- includeEffects(eff2, X, interaction1 = "logdist_cov", name = "friendship")

  eff2 <- includeEffects(eff2, linear, quadratic, name = "alcohol")
  eff2 <- includeEffects(eff2, avAlt, interaction1 = "friendship", name = "alcohol")
  eff2 <- includeEffects(eff2, indeg, interaction1 = "friendship", name = "alcohol")

  alg2 <- sienaAlgorithmCreate(
    projname = file.path("task3_results", "saom_part2_coevolution"),
    n3 = 2500,
    seed = 202
  )

  ans2_path <- file.path("task3_results", "ans_part2_coevolution.rds")
  if (!refit && file.exists(ans2_path)) {
    ans2 <- readRDS(ans2_path)
  } else {
    ans2 <- siena07(alg2, data = dat2, effects = eff2, batch = TRUE, returnDeps = TRUE)
    saveRDS(ans2, file = ans2_path)
  }

  sink(file.path("task3_results", "summary_part2_coevolution.txt"))
  print(summary(ans2))
  sink()

  write_estimates_csv(ans2, file.path("task3_results", "estimates_part2_coevolution.csv"))

  gof2_path <- file.path("task3_results", "gof_part2_coevolution.rds")
  gof2_pdf <- file.path("task3_figures", "gof_part2_coevolution.pdf")
  if (!refit && file.exists(gof2_path) && file.exists(gof2_pdf)) {
    gof2 <- readRDS(gof2_path)
  } else {
    gof2_out <- sienaGOF(ans2, OutdegreeDistribution, varName = "friendship", iterations = gof_iterations, verbose = FALSE)
    gof2_in <- sienaGOF(ans2, IndegreeDistribution, varName = "friendship", iterations = gof_iterations, verbose = FALSE)
    gof2_tri <- sienaGOF(ans2, TriadCensus, varName = "friendship", iterations = gof_iterations, verbose = FALSE)
    gof2_geo <- sienaGOF(ans2, GeodesicDistributionUndirected, varName = "friendship", iterations = gof_iterations, verbose = FALSE)
    gof2_beh <- sienaGOF(ans2, BehaviorDistribution, varName = "alcohol", iterations = gof_iterations, verbose = FALSE)
    gof2 <- list(outdegree = gof2_out, indegree = gof2_in, geodesic = gof2_geo, triad = gof2_tri, alcohol_dist = gof2_beh)
    saveRDS(gof2, file = gof2_path)

    grDevices::pdf(gof2_pdf, width = 11, height = 8.5)
    plot(gof2_out)
    plot(gof2_in)
    plot(gof2_tri)
    plot(gof2_geo)
    plot(gof2_beh)
    grDevices::dev.off()
  }

  write_gof_pvalues_csv(gof2, file.path("task3_results", "gof_part2_pvalues.csv"))
  save_gof_pngs(gof2, prefix = "gof_part2_coevolution", overwrite = refit)

  cat("SAOM estimates saved to task3_results/*.rds\n")
  cat("GOF PDFs saved to task3_figures/\n")
}

if (sys.nframe() == 0L) {
  main()
}
