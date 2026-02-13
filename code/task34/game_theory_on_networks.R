#!/usr/bin/env Rscript
# =============================================================================
# game_theory_on_networks.R
#
# Standalone script for studying how generosity in the Ultimatum Game is
# affected by network topology and evolution-rule heterogeneity.
#
# Combines:
#   - Network generation  (Gomez-Gardenez-Moreno model)
#   - Game dynamics        (Ultimatum Game payoffs via Rcpp)
#   - Strategy update      (Co-evolution of REP / UI / MOR / SP rules via Rcpp)
#   - Simulation pipeline  (Multi-instance loop over parameter grid)
#   - Plotting utilities   (Densities, mean-p-by-k, p-vs-q, rule figures)
#   - Post-hoc analysis    (Grid summary & faceted KDE grid plot)
#
# Usage:
#   Rscript game_theory_on_networks.R          # runs the default grid
#   Rscript game_theory_on_networks.R --help   # prints brief usage info
# =============================================================================

# ---- Load required packages ------------------------------------------------
suppressPackageStartupMessages({
  library(igraph)
  library(Matrix)
  library(Rcpp)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
})

# =============================================================================
# SECTION 1 — NETWORK GENERATION
# =============================================================================

# Gomez-Gardenez-Moreno model
# ---------------------------
# Model utilized by Cardillo et al. to interpolate continuously between a
# homogeneous, Erdos-Renyi-like network (degree distribution ~ uniform,
# alpha = 1) and a heterogeneous, Barabasi-Albert-like network (degree
# distribution ~ power-law, alpha = 0).
#
# Algorithm:
#   - N nodes exist from the start (so uniform linking can reach isolated nodes)
#   - Start from a fully-connected clique of size m0
#   - For each new node i (m0+1 .. N), m links depart from i:
#       * with probability alpha  -> attach uniformly at random
#       * with probability 1-alpha -> attach via preferential attachment
#         (proportional to current degree)
gomez_gardenez_moreno <- function(N, m0 = 5, m = 2, alpha = 0.5,
                                  seed = NULL, max_tries = 1000) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(N > m0, m0 >= 2, m >= 1, alpha >= 0, alpha <= 1)

  g <- make_empty_graph(n = N, directed = FALSE)

  # Initial clique among first m0 nodes
  g <- add_edges(g, as.vector(t(combn(seq_len(m0), 2))))

  # Maintain degrees manually (much faster than calling degree() repeatedly)
  deg <- integer(N)
  deg[seq_len(m0)] <- m0 - 1L

  # Stub list for O(1) preferential-attachment sampling:
  # each node id is repeated deg[id] times
  stubs <- rep.int(seq_len(m0), times = m0 - 1L)

  for (i in (m0 + 1):N) {
    for (e in seq_len(m)) {
      tries <- 0L
      repeat {
        tries <- tries + 1L
        if (tries > max_tries) break

        if (runif(1) < alpha || length(stubs) == 0) {
          # Uniform attachment
          j <- sample.int(N, 1)
        } else {
          # Preferential attachment proportional to current degree
          j <- stubs[sample.int(length(stubs), 1)]
        }

        if (j == i) next
        if (are_adjacent(g, i, j)) next

        # Accept link
        g <- add_edges(g, c(i, j))
        deg[i] <- deg[i] + 1L
        deg[j] <- deg[j] + 1L
        stubs <- c(stubs, i, j)
        break
      }
    }
  }

  # Should already be simple; keep as guardrail
  simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
}

# =============================================================================
# SECTION 2 — ULTIMATUM GAME PAYOFFS (Rcpp)
# =============================================================================

# Efficient computation of Ultimatum Game payoffs over all edges.
#   u, v : integer vectors of length M (edge list, 1-indexed)
#   p, q : numeric vectors of length N (offer / acceptance threshold)
#   returns: payoff vector of length N
Rcpp::cppFunction(code = '
  Rcpp::NumericVector ug_payoffs_general(Rcpp::IntegerVector u,
                                        Rcpp::IntegerVector v,
                                        Rcpp::NumericVector p,
                                        Rcpp::NumericVector q) {
    int M = u.size();
    int N = p.size();
    Rcpp::NumericVector payoff(N);

    for (int e = 0; e < M; ++e) {
      int i = u[e] - 1;
      int j = v[e] - 1;

      double pi = p[i], pj = p[j];
      double qi = q[i], qj = q[j];

      // i proposes to j
      if (pi >= qj) {
        payoff[i] += (1.0 - pi);
        payoff[j] += pi;
      }
      // j proposes to i
      if (pj >= qi) {
        payoff[j] += (1.0 - pj);
        payoff[i] += pj;
      }
    }
    return payoff;
  }
')

# =============================================================================
# SECTION 3 — CO-EVOLUTIONARY STRATEGY UPDATE RULES (Rcpp)
# =============================================================================

# Four local* update rules co-evolve:
#   REP (0) — Replicator: copy neighbor j's strategy with probability
#             proportional to (Pj - Pi) / (norm_factor * max(ki, kj))
#   UI  (1) — Unconditional Imitation: copy the best-performing neighbor
#   MOR (2) — Moran-like: copy a neighbor with probability proportional
#             to that neighbor's payoff (fitness-proportional selection)
#   SP  (3) — Social Penalty: the node with the minimum payoff in the
#             whole network, and its neighbors, randomize their strategy
#             and are assigned a random rule from the allowed subset.
#             (* SP uses global information but has local effect)
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

enum Rule : int { REP=0, UI=1, MOR=2, SP=3 };

// ---- REP update for a single node i ----
inline void update_rep_one(
    int i,
    const NumericVector &p,
    const NumericVector &q,
    const NumericVector &payoff,
    const IntegerVector &ptr,
    const IntegerVector &idx,
    NumericVector &p_new,
    NumericVector &q_new,
    IntegerVector &rule_new,
    const IntegerVector &rule,
    double norm_factor
){
  const int start_i = ptr[i];
  const int end_i   = ptr[i + 1];
  const int ki      = end_i - start_i;
  if (ki <= 0) return;

  int r = (int)std::floor(R::runif(0.0, (double)ki));
  if (r >= ki) r = ki - 1;
  const int j = idx[start_i + r];

  const double Pi = payoff[i];
  const double Pj = payoff[j];

  if (Pj > Pi) {
    const int kj = ptr[j + 1] - ptr[j];
    const int kmax = (ki > kj) ? ki : kj;

    double prob = (Pj - Pi) / (norm_factor * (double)kmax);
    if (prob < 0.0) prob = 0.0;
    if (prob > 1.0) prob = 1.0;

    if (R::runif(0.0, 1.0) < prob) {
      p_new[i]    = p[j];
      q_new[i]    = q[j];
      rule_new[i] = rule[j];
    }
  }
}

// ---- UI update for a single node i ----
inline void update_ui_one(
    int i,
    const NumericVector &p,
    const NumericVector &q,
    const NumericVector &payoff,
    const IntegerVector &ptr,
    const IntegerVector &idx,
    NumericVector &p_new,
    NumericVector &q_new,
    IntegerVector &rule_new,
    const IntegerVector &rule
){
  const int start_i = ptr[i];
  const int end_i   = ptr[i + 1];
  if (end_i <= start_i) return;

  double bestP = payoff[i];
  int bestJ = -1;

  for (int t = start_i; t < end_i; ++t) {
    const int j = idx[t];
    const double Pj = payoff[j];
    if (Pj > bestP) {
      bestP = Pj;
      bestJ = j;
    }
  }

  if (bestJ >= 0) {
    p_new[i]    = p[bestJ];
    q_new[i]    = q[bestJ];
    rule_new[i] = rule[bestJ];
  }
}

// ---- MOR (Moran-like) update for a single node i ----
inline void update_mor_one(
    int i,
    const NumericVector &p,
    const NumericVector &q,
    const NumericVector &payoff,
    const IntegerVector &ptr,
    const IntegerVector &idx,
    NumericVector &p_new,
    NumericVector &q_new,
    IntegerVector &rule_new,
    const IntegerVector &rule
){
  const int start_i = ptr[i];
  const int end_i   = ptr[i + 1];
  if (end_i <= start_i) return;

  double sumW = 0.0;
  for (int t = start_i; t < end_i; ++t) {
    const int j = idx[t];
    const double w = payoff[j];
    if (w > 0.0) sumW += w;
  }
  if (sumW <= 0.0) return;

  double u = R::runif(0.0, sumW);
  for (int t = start_i; t < end_i; ++t) {
    const int j = idx[t];
    const double w = payoff[j];
    if (w <= 0.0) continue;
    u -= w;
    if (u <= 0.0) {
      p_new[i]    = p[j];
      q_new[i]    = q[j];
      rule_new[i] = rule[j];
      return;
    }
  }
}

// [[Rcpp::export]]
List coevolution_update(
    NumericVector p,
    NumericVector q,
    NumericVector payoff,
    IntegerVector rule,          // 0=REP, 1=UI, 2=MOR, 3=SP
    IntegerVector ptr,
    IntegerVector idx,
    IntegerVector rule_subset,   // allowed rules when SP triggers
    double norm_factor = 2.0
){
  RNGScope scope;

  const int N = p.size();
  if (q.size() != N || payoff.size() != N || rule.size() != N) {
    stop("p, q, payoff, rule must have the same length");
  }
  if (ptr.size() != N + 1) stop("ptr must have length N+1");
  if (rule_subset.size() <= 0) stop("rule_subset must be non-empty");

  NumericVector p_new = clone(p);
  NumericVector q_new = clone(q);
  IntegerVector rule_new = clone(rule);

  // 1) Per-node update for REP/UI/MOR; SP nodes skip this step
  for (int i = 0; i < N; ++i) {
    const int ri = rule[i];

    if (ri == REP) {
      update_rep_one(i, p, q, payoff, ptr, idx, p_new, q_new, rule_new, rule, norm_factor);
    } else if (ri == UI) {
      update_ui_one(i, p, q, payoff, ptr, idx, p_new, q_new, rule_new, rule);
    } else if (ri == MOR) {
      update_mor_one(i, p, q, payoff, ptr, idx, p_new, q_new, rule_new, rule);
    } else if (ri == SP) {
      // handled in the global SP step below
    } else {
      stop("Unknown rule code at node %d: %d", i + 1, ri);
    }
  }

  // 2) Global SP step: find node with minimum payoff
  int imin = 0;
  double pmin = payoff[0];
  for (int i = 1; i < N; ++i) {
    const double Pi = payoff[i];
    if (Pi < pmin) {
      pmin = Pi;
      imin = i;
    }
  }

  // Build S = {imin} ∪ neighbors(imin)
  std::vector<unsigned char> inS(N, (unsigned char)0);
  inS[imin] = 1;

  const int start_m = ptr[imin];
  const int end_m   = ptr[imin + 1];
  for (int t = start_m; t < end_m; ++t) {
    const int j = idx[t];
    if (j >= 0 && j < N) inS[j] = 1;
  }

  // Apply SP override: nodes in S that currently have rule==SP
  // get randomized strategy and a randomly chosen rule from rule_subset
  const int mR = rule_subset.size();
  for (int i = 0; i < N; ++i) {
    if (!inS[i]) continue;
    if (rule[i] != SP) continue;

    p_new[i] = R::runif(0.0, 1.0);
    q_new[i] = R::runif(0.0, 1.0);

    int rpick = (int)std::floor(R::runif(0.0, (double)mR));
    if (rpick >= mR) rpick = mR - 1;
    rule_new[i] = rule_subset[rpick];
  }

  return List::create(_["p"] = p_new, _["q"] = q_new, _["rule"] = rule_new,
                      _["imin"] = imin, _["min_payoff"] = pmin);
}
')

# =============================================================================
# SECTION 4 — PLOTTING UTILITIES
# =============================================================================

# ---- Overlay densities (base R) ----
# Plot kernel-density estimates of multiple numeric vectors on the same axes.
plot_overlay_densities <- function(xs,
                                  main = "Overlayed densities",
                                  xlab = "value",
                                  lwd = 2,
                                  legend_pos = "topright",
                                  ...) {
  stopifnot(is.list(xs), length(xs) >= 1)
  if (is.null(names(xs)) || any(names(xs) == "")) {
    names(xs) <- if (is.null(names(xs))) paste0("x", seq_along(xs)) else {
      nm <- names(xs); nm[nm == ""] <- paste0("x", which(nm == "")); nm
    }
  }

  dens_list <- lapply(xs, density, adjust = 2, ...)

  xlim <- range(unlist(lapply(dens_list, `[[`, "x")))
  ylim <- range(unlist(lapply(dens_list, `[[`, "y")))

  cols <- seq_along(xs)

  plot(dens_list[[1]], xlim = xlim, ylim = ylim,
       col = cols[1], lwd = lwd,
       main = main, xlab = xlab)

  if (length(xs) > 1) {
    for (i in 2:length(xs)) lines(dens_list[[i]], col = cols[i], lwd = lwd)
  }

  legend(legend_pos, legend = names(xs), col = cols, lwd = lwd, bty = "n")

  invisible(dens_list)
}

# ---- Filename / path helpers ----

# Make a filesystem-safe name from arbitrary text
sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._=-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# Build a compact tag encoding all simulation parameters (used in filenames)
params_tag <- function(params) {
  snaps <- paste(params$snapshot_gens, collapse = "-")
  rules <- paste(params$rule_subset, collapse = "-")
  tag <- paste0(
    "inst=", params$n_instances,
    "__N=", params$N,
    "__m0=", params$m0,
    "__m=", params$m,
    "__alpha=", params$alpha,
    "__gens=", params$n_generations,
    "__snaps=", snaps,
    "__rules=", rules,
    "__norm=", params$norm_factor,
    "__seed0=", params$seed0
  )
  sanitize_filename(tag)
}

# Create a timestamped output directory
timestamped_plot_dir <- function(base_dir = "plots", prefix = "run") {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  outdir <- file.path(base_dir, paste0(prefix, "_", ts))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir
}

# ---- Save helpers ----

# Save a base-R plot to PNG with slightly enlarged text
save_base_plot <- function(path, width = 7, height = 5, dpi = 200, fun) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(filename = path, width = width, height = height, units = "in", res = dpi)

  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(cex = 1.10, cex.axis = 1.10, cex.lab = 1.10, cex.main = 1.10)

  fun()
  invisible(path)
}

# Save a ggplot to PNG, stripping the title for cleaner publication figures
save_ggplot_no_title <- function(path, p, width = 7, height = 5, dpi = 200) {
  p2 <- p + theme(text = element_text(size = 13)) +
    labs(title = NULL) + theme(plot.title = element_blank())
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = path, plot = p2, width = width, height = height, dpi = dpi)
  invisible(path)
}

# ---- Individual plot functions for a simulation result object (res) ----

# Pooled p / q density snapshots across all instances at each snapshot time
plot_snapshot_densities <- function(res, outdir, save = TRUE, render = TRUE, dpi = 200) {
  files <- character(0)
  tag <- params_tag(res$params)

  if (render) {
    plot_overlay_densities(res$snapshots$p,
      main = sprintf("Pooled p-distribution (%d instances)", res$params$n_instances)
    )
    plot_overlay_densities(res$snapshots$q,
      main = sprintf("Pooled q-distribution (%d instances)", res$params$n_instances)
    )
  }

  if (save) {
    f1 <- file.path(outdir, paste0("densities_p__", tag, ".png"))
    f2 <- file.path(outdir, paste0("densities_q__", tag, ".png"))

    save_base_plot(f1, dpi = dpi, fun = function() {
      plot_overlay_densities(res$snapshots$p, main = "")
    })
    save_base_plot(f2, dpi = dpi, fun = function() {
      plot_overlay_densities(res$snapshots$q, main = "")
    })
    files <- c(files, f1, f2)
  }

  files
}

# Mean offer p as a function of node degree k (final state)
plot_final_mean_p_by_k <- function(res, outdir, save = TRUE, render = TRUE, dpi = 200) {
  files <- character(0)
  tag <- params_tag(res$params)

  df <- res$summaries$mean_p_by_k
  k_vals <- df$k
  mean_p <- df$mean_p

  if (render) {
    plot(k_vals, mean_p,
      pch = 16, xlab = "degree k", ylab = "mean p at final time",
      main = sprintf("E[p | k] over %d instances", res$params$n_instances)
    )
    grid()
  }

  if (save) {
    f <- file.path(outdir, paste0("mean_p_by_k__", tag, ".png"))
    save_base_plot(f, dpi = dpi, fun = function() {
      plot(k_vals, mean_p, pch = 16,
           xlab = "degree k", ylab = "mean p at final time", main = "")
      grid()
    })
    files <- c(files, f)
  }

  files
}

# Scatter plot of final p vs q across all nodes and instances
plot_final_p_vs_q <- function(res, outdir, save = TRUE, render = TRUE, dpi = 200) {
  files <- character(0)
  tag <- params_tag(res$params)

  p <- res$final$p
  q <- res$final$q

  if (render) {
    plot(p, q,
      pch = 16, xlab = "p (final)", ylab = "q (final)",
      main = sprintf("p vs q (%d instances)", res$params$n_instances)
    )
    grid()
  }

  if (save) {
    f <- file.path(outdir, paste0("p_vs_q__", tag, ".png"))
    save_base_plot(f, dpi = dpi, fun = function() {
      plot(p, q, pch = 16,
           xlab = "p (final)", ylab = "q (final)", main = "")
      grid()
    })
    files <- c(files, f)
  }

  files
}

# Multi-rule figures: <p>_rule(t), rule abundances, neighbor-rule heatmap
# Only produced when multiple update rules co-evolve.
plot_multi_rule_figures <- function(res, outdir, save = TRUE, render = TRUE, dpi = 200) {
  files <- character(0)
  if (is.null(res$summaries$df_p_rule_t) ||
      is.null(res$summaries$df_counts) ||
      is.null(res$summaries$joint_df)) {
    return(files)
  }

  theme_set(theme_minimal(base_size = 13))
  tag <- params_tag(res$params)
  snapshot_gens <- res$params$snapshot_gens
  rule_levels <- levels(res$summaries$df_p_rule_t$rule)

  pd <- position_dodge(width = 0.12)

  # Mean offer <p> by update rule over time (± sd)
  df_p_rule_t <- res$summaries$df_p_rule_t
  g_p <- ggplot(df_p_rule_t, aes(x = t, y = mean, color = rule, group = rule)) +
    geom_line() +
    geom_point(position = pd) +
    geom_errorbar(
      aes(ymin = pmax(0, mean - sd), ymax = pmin(1, mean + sd)),
      width = 0.08, linewidth = 0.35, position = pd
    ) +
    scale_x_log10(breaks = snapshot_gens) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = sprintf("<p> by update rule over time (±sd; subset=%s; %d instances)",
        paste(rule_levels, collapse = ","), res$params$n_instances
      ),
      x = "generation t (log scale)",
      y = "<p>_rule(t)",
      color = "rule"
    )

  if (render) print(g_p)
  if (save) {
    f <- file.path(outdir, paste0("p_by_rule_over_time__", tag, ".png"))
    save_ggplot_no_title(f, g_p, dpi = dpi)
    files <- c(files, f)
  }

  # Rule abundances (fraction of nodes) over time
  df_counts <- res$summaries$df_counts
  g_ab <- ggplot(df_counts, aes(x = t, y = frac, color = rule, group = rule)) +
    geom_line() +
    geom_point() +
    scale_x_log10(breaks = snapshot_gens) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = sprintf("Rule abundances over time (subset=%s; %d instances)",
        paste(rule_levels, collapse = ","), res$params$n_instances
      ),
      x = "generation t (log scale)",
      y = "fraction of nodes",
      color = "rule"
    )

  if (render) print(g_ab)
  if (save) {
    f <- file.path(outdir, paste0("rule_abundances_over_time__", tag, ".png"))
    save_ggplot_no_title(f, g_ab, dpi = dpi)
    files <- c(files, f)
  }

  # Joint neighbor-rule distribution heatmap
  joint_df <- res$summaries$joint_df
  grid_df <- expand.grid(
    rule_i = rule_levels,
    rule_j = rule_levels,
    stringsAsFactors = FALSE
  )
  joint_df <- merge(grid_df, joint_df, by = c("rule_i", "rule_j"), all.x = TRUE)
  joint_df$prob[is.na(joint_df$prob)] <- 0
  joint_df$rule_i <- factor(joint_df$rule_i, levels = rule_levels)
  joint_df$rule_j <- factor(joint_df$rule_j, levels = rule_levels)

  p_heat <- ggplot(joint_df, aes(x = rule_i, y = rule_j, fill = prob)) +
    geom_tile() +
    scale_fill_gradient(name = "P(edge endpoint pair)") +
    labs(
      title = sprintf("Neighbor rule joint distribution (subset=%s; %d instances)",
        paste(rule_levels, collapse = ","), res$params$n_instances
      ),
      x = "rule of node i",
      y = "rule of neighbor j"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    coord_fixed() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  if (render) print(p_heat)
  if (save) {
    f <- file.path(outdir, paste0("neighbor_rule_joint__", tag, ".png"))
    save_ggplot_no_title(f, p_heat, dpi = dpi)
    files <- c(files, f)
  }

  files
}

# Master plotting function: runs all individual plotters and saves res.rds
plot_all <- function(res, base_dir = "plots", prefix = "run",
                     save = TRUE, render = TRUE, dpi = 200) {
  # Slightly larger text in base plots (rendered)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(cex = 1.10, cex.axis = 1.10, cex.lab = 1.10, cex.main = 1.10)

  outdir <- timestamped_plot_dir(base_dir = base_dir, prefix = prefix)
  files <- character(0)

  files <- c(files, plot_snapshot_densities(res, outdir, save = save, render = render, dpi = dpi))
  files <- c(files, plot_final_mean_p_by_k(res, outdir, save = save, render = render, dpi = dpi))
  files <- c(files, plot_final_p_vs_q(res, outdir, save = save, render = render, dpi = dpi))
  files <- c(files, plot_multi_rule_figures(res, outdir, save = save, render = render, dpi = dpi))

  saveRDS(res, file = file.path(outdir, "res.rds"))

  invisible(list(outdir = outdir, files = files))
}

# =============================================================================
# SECTION 5 — FULL SINGLE-CONFIGURATION SIMULATION PIPELINE
# =============================================================================

# Run the co-evolutionary Ultimatum Game pipeline for a single parameter
# configuration, looping over n_instances independent graph realizations.
# Returns a structured list with:
#   $params        — input parameters
#   $elapsed_secs  — wall-clock runtime
#   $snapshots     — pooled p, q (and rule) distributions at snapshot times
#   $final         — pooled final-state vectors (p, q, k, rule, edge rules)
#   $summaries     — derived statistics (mean_p_by_k, df_p_rule_t, etc.)
run_coevo_pipeline <- function(
  n_instances   = 100,
  N             = 1000,
  m0            = 6,
  m             = 2,
  alpha         = 1,
  n_generations = 20000,
  snapshot_gens = as.integer(c(1, 2, 5, 10, 100, 1000, 10000, 20000)),
  rule_subset   = c(3L),
  seed0         = 1000,
  norm_factor   = 2.0,
  verbose       = TRUE
) {
  # Rule mapping: integer id <-> human-readable name
  RULES <- data.frame(
    id   = c(0L, 1L, 2L, 3L),
    name = c("REP", "UI", "MOR", "SP"),
    stringsAsFactors = FALSE
  )

  rule_levels_id   <- sort(unique(as.integer(rule_subset)))
  rule_levels_name <- RULES$name[match(rule_levels_id, RULES$id)]

  rule_id_to_name <- function(x) {
    nm <- RULES$name[match(as.integer(x), RULES$id)]
    ifelse(is.na(nm), paste0("rule=", x), nm)
  }

  snapshot_gens <- as.integer(snapshot_gens)
  snap_names    <- paste0("t=", snapshot_gens)

  # Storage: pooled samples per snapshot time
  snapshot_p_all <- setNames(vector("list", length(snapshot_gens)), snap_names)
  snapshot_q_all <- setNames(vector("list", length(snapshot_gens)), snap_names)
  for (nm in snap_names) {
    snapshot_p_all[[nm]] <- numeric(0)
    snapshot_q_all[[nm]] <- numeric(0)
  }

  # Snapshot storage for rules (only when multiple rules co-evolve)
  multi_rules <- (length(rule_subset) > 1L)
  snapshot_rule_all <- setNames(vector("list", length(snapshot_gens)), snap_names)
  for (nm in snap_names) snapshot_rule_all[[nm]] <- integer(0)

  # Final pooled vectors
  p_final_all     <- numeric(0)
  q_final_all     <- numeric(0)
  k_final_all     <- integer(0)
  rule_final_all  <- integer(0)
  rule_edge_i_all <- integer(0)
  rule_edge_j_all <- integer(0)

  start_time <- Sys.time()

  for (rep in seq_len(n_instances)) {
    set.seed(seed0 + rep)

    # Generate network
    g <- gomez_gardenez_moreno(N = N, m0 = m0, m = m, alpha = alpha)

    # Edge list for UG payoff computation
    el <- igraph::as_edgelist(g, names = FALSE)
    u <- el[, 1]; v <- el[, 2]
    rm(el); gc()

    # CSR adjacency for coevolution_update
    A <- as_adjacency_matrix(g, sparse = TRUE)
    A <- as(A, "RsparseMatrix")
    ptr <- A@p
    idx <- A@j

    # Initialize strategies uniformly in [0, 1]
    p <- runif(N)
    q <- runif(N)

    # Initialize rules
    if (length(rule_subset) == 1L) {
      rule <- rep.int(as.integer(rule_subset[1L]), N)
    } else {
      rule <- sample(as.integer(rule_subset), N, replace = TRUE)
    }

    # Evolutionary dynamics
    for (n_gen in seq_len(n_generations)) {
      # Snapshot at specified generations
      if (n_gen %in% snapshot_gens) {
        nm <- paste0("t=", n_gen)
        snapshot_p_all[[nm]] <- c(snapshot_p_all[[nm]], p)
        snapshot_q_all[[nm]] <- c(snapshot_q_all[[nm]], q)
        if (multi_rules) snapshot_rule_all[[nm]] <- c(snapshot_rule_all[[nm]], rule)
      }

      # Compute payoffs and update strategies
      payoff <- ug_payoffs_general(u, v, p, q)
      res <- coevolution_update(p, q, payoff, rule, ptr, idx,
                                rule_subset, norm_factor = norm_factor)
      p <- res$p; q <- res$q; rule <- res$rule
    }

    # Collect final state
    k <- igraph::degree(g)
    p_final_all <- c(p_final_all, p)
    q_final_all <- c(q_final_all, q)
    k_final_all <- c(k_final_all, k)

    if (multi_rules) {
      rule_final_all  <- c(rule_final_all, rule)
      rule_edge_i_all <- c(rule_edge_i_all, rule[u])
      rule_edge_j_all <- c(rule_edge_j_all, rule[v])
    }

    rm(g, A, ptr, idx, u, v, p, q, rule, payoff, res)
    gc()
  }

  end_time <- Sys.time()
  elapsed_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
  if (verbose) {
    cat("loop execution time with N:", N,
        ", n_generations:", n_generations,
        ", n_instances:", n_instances,
        ", secs:", elapsed_secs, "\n")
  }

  # ---- Derived summaries ----
  summaries <- list()

  # Mean p by degree
  df_pk <- data.frame(k = k_final_all, p = p_final_all)
  mean_p_by_k <- tapply(df_pk$p, df_pk$k, mean)
  k_vals <- as.integer(names(mean_p_by_k))
  summaries$mean_p_by_k <- data.frame(k = k_vals, mean_p = as.numeric(mean_p_by_k))

  if (multi_rules) {
    # Long data.frame over snapshots (only rules in subset)
    df_snap <- do.call(rbind, lapply(seq_along(snapshot_gens), function(i) {
      tgen <- snapshot_gens[i]
      nm <- paste0("t=", tgen)
      p_vec <- snapshot_p_all[[nm]]
      r_vec <- snapshot_rule_all[[nm]]

      keep <- r_vec %in% rule_levels_id
      p_vec <- p_vec[keep]
      r_vec <- r_vec[keep]

      data.frame(
        t = rep.int(tgen, length(p_vec)),
        rule_id = r_vec,
        p = p_vec,
        stringsAsFactors = FALSE
      )
    }))
    df_snap$rule <- factor(rule_id_to_name(df_snap$rule_id), levels = rule_levels_name)
    summaries$df_snap <- df_snap

    # <p>_rule(t) with sd and se
    df_mean <- aggregate(p ~ t + rule, data = df_snap, FUN = mean)
    df_sd   <- aggregate(p ~ t + rule, data = df_snap, FUN = sd)
    df_n    <- aggregate(p ~ t + rule, data = df_snap, FUN = length)
    names(df_mean)[3] <- "mean"
    names(df_sd)[3]   <- "sd"
    names(df_n)[3]    <- "n"

    df_p_rule_t <- merge(merge(df_mean, df_sd, by = c("t", "rule"), all = TRUE),
                         df_n, by = c("t", "rule"), all = TRUE)
    df_p_rule_t$sd[is.na(df_p_rule_t$sd)] <- 0
    df_p_rule_t$se <- df_p_rule_t$sd / sqrt(df_p_rule_t$n)
    summaries$df_p_rule_t <- df_p_rule_t

    # Rule abundances (fractions) over time
    df_counts_obs <- aggregate(p ~ t + rule, data = df_snap, FUN = length)
    names(df_counts_obs)[3] <- "count"

    df_grid <- expand.grid(
      t = snapshot_gens,
      rule = factor(rule_levels_name, levels = rule_levels_name)
    )
    df_counts <- merge(df_grid, df_counts_obs, by = c("t", "rule"), all.x = TRUE)
    df_counts$count[is.na(df_counts$count)] <- 0L
    df_counts$total <- ave(df_counts$count, df_counts$t, FUN = sum)
    df_counts$frac <- ifelse(df_counts$total > 0,
                             df_counts$count / df_counts$total, 0)
    summaries$df_counts <- df_counts

    # Joint neighbor-rule distribution
    r1 <- c(rule_edge_i_all, rule_edge_j_all)
    r2 <- c(rule_edge_j_all, rule_edge_i_all)

    joint_mat <- table(
      factor(r1, levels = rule_levels_id),
      factor(r2, levels = rule_levels_id)
    )

    joint_df <- as.data.frame(joint_mat, stringsAsFactors = FALSE)
    colnames(joint_df) <- c("id_i", "id_j", "count")
    joint_df$prob <- joint_df$count / sum(joint_df$count)
    joint_df$rule_i <- factor(rule_id_to_name(joint_df$id_i), levels = rule_levels_name)
    joint_df$rule_j <- factor(rule_id_to_name(joint_df$id_j), levels = rule_levels_name)
    summaries$joint_df <- joint_df
  }

  invisible(list(
    params = list(
      n_instances = n_instances, N = N, m0 = m0, m = m, alpha = alpha,
      n_generations = n_generations, snapshot_gens = snapshot_gens,
      rule_subset = rule_subset, seed0 = seed0, norm_factor = norm_factor
    ),
    elapsed_secs = elapsed_secs,
    snapshots = list(
      p = snapshot_p_all,
      q = snapshot_q_all,
      rule = if (multi_rules) snapshot_rule_all else NULL
    ),
    final = list(
      p = p_final_all,
      q = q_final_all,
      k = k_final_all,
      rule      = if (multi_rules) rule_final_all  else NULL,
      rule_edge_i = if (multi_rules) rule_edge_i_all else NULL,
      rule_edge_j = if (multi_rules) rule_edge_j_all else NULL
    ),
    summaries = summaries
  ))
}

# =============================================================================
# SECTION 6 — POST-HOC ANALYSIS: SUMMARIZE GRID FOLDERS
# =============================================================================

# Scan a base directory of simulation output folders (named like
# N_<N>_alpha_<alpha>_rules_<rules>) and write a human-readable summary
# of all parameter combinations found.
summarize_plots_grid <- function(base = "plots_grid") {
  if (!dir.exists(base)) stop(paste("Folder not found:", base))

  pat <- "^N_(\\d+)_alpha_(\\d+(?:p\\d+|\\.\\d+)?)_rules_(\\d+(?:-\\d+)*)$"
  alpha_to_num_local <- function(a) as.numeric(gsub("p", ".", a))

  folders <- list.dirs(base, full.names = FALSE, recursive = FALSE)
  folders <- folders[folders != ""]
  folders <- sort(folders)

  parsed <- lapply(folders, function(nm) {
    m <- regexec(pat, nm, perl = TRUE)
    g <- regmatches(nm, m)[[1]]
    if (length(g) == 0) return(NULL)
    list(
      folder    = nm,
      N         = as.integer(g[2]),
      alpha_raw = g[3],
      alpha     = alpha_to_num_local(g[3]),
      rules_raw = g[4]
    )
  })

  parsed <- Filter(Negate(is.null), parsed)
  skipped <- setdiff(folders, vapply(parsed, `[[`, "", "folder"))

  if (length(parsed) == 0) stop("No matching folders found.")

  df <- data.frame(
    folder    = vapply(parsed, `[[`, "", "folder"),
    N         = vapply(parsed, `[[`, integer(1), "N"),
    alpha_raw = vapply(parsed, `[[`, "", "alpha_raw"),
    alpha     = vapply(parsed, `[[`, numeric(1), "alpha"),
    rules_raw = vapply(parsed, `[[`, "", "rules_raw"),
    stringsAsFactors = FALSE
  )

  # Unique combinations
  uniq <- unique(df[, c("N", "alpha", "rules_raw")])
  uniq <- uniq[order(uniq$N, uniq$alpha, uniq$rules_raw), ]

  out_path <- file.path(base, "combinations_summary.txt")
  con <- file(out_path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)

  writeLines(paste0("Base directory: ", normalizePath(base)), con)
  writeLines(paste0("Matched folders: ", nrow(df)), con)
  writeLines(paste0("Unique combinations (N, alpha, rules): ", nrow(uniq)), con)
  writeLines("", con)

  writeLines("=== Unique combinations ===", con)
  apply(uniq, 1, function(r) {
    writeLines(sprintf("- N=%s\talpha=%s\trules=%s",
      r[["N"]], format(as.numeric(r[["alpha"]]), trim = TRUE), r[["rules_raw"]]), con)
  })

  writeLines("\n=== Grouped by N ===", con)
  for (Nval in sort(unique(df$N))) {
    sub <- df[df$N == Nval, ]
    alphas <- sub$alpha_raw[order(alpha_to_num_local(sub$alpha_raw))]
    alphas <- unique(alphas)
    rulesets <- sort(unique(sub$rules_raw))
    writeLines(sprintf("\nN=%d", Nval), con)
    writeLines(paste0("  alphas: ", paste(alphas, collapse = ", ")), con)
    writeLines(paste0("  rules:  ", paste(rulesets, collapse = ", ")), con)
  }

  writeLines("\n=== Grouped by alpha ===", con)
  alphas_all <- unique(df$alpha_raw)
  alphas_all <- alphas_all[order(alpha_to_num_local(alphas_all))]
  for (aval in alphas_all) {
    sub <- df[df$alpha_raw == aval, ]
    Ns <- sort(unique(sub$N))
    rulesets <- sort(unique(sub$rules_raw))
    writeLines(sprintf("\nalpha=%s (raw=%s)",
      format(alpha_to_num_local(aval), trim = TRUE), aval), con)
    writeLines(paste0("  N:      ", paste(Ns, collapse = ", ")), con)
    writeLines(paste0("  rules:  ", paste(rulesets, collapse = ", ")), con)
  }

  writeLines("\n=== Grouped by rules ===", con)
  for (r in sort(unique(df$rules_raw))) {
    sub <- df[df$rules_raw == r, ]
    Ns <- sort(unique(sub$N))
    alphas <- unique(sub$alpha_raw)
    alphas <- alphas[order(alpha_to_num_local(alphas))]
    writeLines(sprintf("\nrules=%s", r), con)
    writeLines(paste0("  N:      ", paste(Ns, collapse = ", ")), con)
    writeLines(paste0("  alphas: ", paste(alphas, collapse = ", ")), con)
  }

  if (length(skipped) > 0) {
    writeLines("\n=== Skipped (did not match pattern) ===", con)
    writeLines(paste0("- ", skipped), con)
  }

  cat("Wrote:", out_path, "\n")
  invisible(out_path)
}

# =============================================================================
# SECTION 7 — POST-HOC ANALYSIS: FACETED KDE GRID PLOT
# =============================================================================

# Build a faceted grid of KDE densities of final offer p across
# (alpha x ruleset) combinations for a given N, coloring backgrounds
# by mean(p). Reads saved res.rds files produced by the simulation grid.
plot_grid_final_p <- function(
  base_dir       = "plots_grid",
  N_target       = 1000L,
  alpha_raw_vals = c("0", "0p25", "0p5", "0p75", "1"),
  rules_raw_vals = c("0", "0-1-2", "0-1-2-3", "0-2", "0-3", "1-2", "3"),
  dens_n         = 401,     # number of x points for KDE curve

  dens_adjust    = 1.5,     # >1 = smoother; try 1.2–2.5
  eps_sd         = 1e-9     # protects against sd = 0 in weights
) {
  alpha_to_num_local <- function(a) as.numeric(gsub("p", ".", a))

  make_base_folder <- function(N, alpha_raw, rules_raw) {
    paste0("N_", N, "_alpha_", alpha_raw, "_rules_", rules_raw)
  }

  # Inside each base folder, find the most recently modified timestamp subfolder
  find_latest_run_dir <- function(base_folder_path, base_folder_name) {
    if (!dir.exists(base_folder_path)) return(NA_character_)
    subs <- list.dirs(base_folder_path, full.names = TRUE, recursive = FALSE)
    if (length(subs) == 0) return(NA_character_)

    keep <- subs[basename(subs) %>% str_starts(paste0(base_folder_name, "_"))]
    if (length(keep) == 0) return(NA_character_)

    info <- file.info(keep)
    keep[which.max(info$mtime)]
  }

  load_final_p <- function(res_rds_path) {
    res <- readRDS(res_rds_path)
    p <- res$final$p
    if (is.null(p) || !is.numeric(p))
      stop("res$final$p missing or not numeric in: ", res_rds_path)
    p
  }

  make_xgrid <- function(dens_n) seq(0, 1, length.out = dens_n)

  kde_density_df <- function(p, xgrid, adjust = 1.5) {
    d <- density(p, from = 0, to = 1, n = length(xgrid), adjust = adjust, cut = 0)
    tibble(x = d$x, density = d$y)
  }

  # Output paths
  out_plot <- file.path(base_dir, sprintf("grid_final_p_N_%d.png", N_target))
  out_csv  <- file.path(base_dir, sprintf("grid_final_p_N_%d_stats.csv", N_target))

  # All combinations to plot
  combos <- expand.grid(
    alpha_raw = alpha_raw_vals,
    rules_raw = rules_raw_vals,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      N = N_target,
      alpha = alpha_to_num_local(alpha_raw),
      base_folder = pmap_chr(list(N, alpha_raw, rules_raw), make_base_folder),
      base_path = file.path(base_dir, base_folder)
    )

  # Resolve run directory and res.rds path
  combos <- combos %>%
    rowwise() %>%
    mutate(
      run_dir  = find_latest_run_dir(base_path, base_folder),
      res_path = if (!is.na(run_dir)) file.path(run_dir, "res.rds") else NA_character_,
      exists   = !is.na(res_path) && file.exists(res_path)
    ) %>%
    ungroup()

  missing <- combos %>% filter(!exists)
  if (nrow(missing) > 0) {
    message("WARNING: Missing res.rds for ", nrow(missing),
            " combinations. They will be skipped.")
    message(paste0("  - ", missing$base_folder[1:min(10, nrow(missing))],
                   collapse = "\n"))
  }

  combos_ok <- combos %>% filter(exists)
  if (nrow(combos_ok) == 0) stop("No res.rds files found for the requested subset.")

  # Load data and compute per-cell statistics
  xgrid <- make_xgrid(dens_n)

  cell_list <- lapply(seq_len(nrow(combos_ok)), function(i) {
    row <- combos_ok[i, ]
    p <- load_final_p(row$res_path)
    p <- p[p >= 0 & p <= 1]  # clamp against numeric noise

    stats <- tibble(
      N = row$N, alpha_raw = row$alpha_raw, alpha = row$alpha,
      rules_raw = row$rules_raw, mean_p = mean(p), sd_p = sd(p)
    )

    dens <- kde_density_df(p, xgrid, adjust = dens_adjust) %>%
      mutate(
        N = row$N, alpha_raw = row$alpha_raw, alpha = row$alpha,
        rules_raw = row$rules_raw, mean_p = stats$mean_p, sd_p = stats$sd_p
      )

    list(stats = stats, dens = dens)
  })

  stats_df <- bind_rows(lapply(cell_list, `[[`, "stats"))
  dens_df  <- bind_rows(lapply(cell_list, `[[`, "dens"))

  write_csv(stats_df, out_csv)

  # Sort rules (rows) by weighted mean of p across alphas (weight = 1/var)
  rule_order <- stats_df %>%
    mutate(var_p = pmax(sd_p^2, eps_sd)) %>%
    mutate(w = 1 / var_p) %>%
    group_by(rules_raw) %>%
    summarise(wmean_p = sum(w * mean_p) / sum(w), .groups = "drop") %>%
    arrange(wmean_p) %>%
    pull(rules_raw)

  # Factor levels for plotting
  alpha_labels <- rev(sprintf("%.2f", alpha_to_num_local(alpha_raw_vals)))
  dens_df <- dens_df %>%
    mutate(
      alpha_lab = factor(sprintf("%.2f", alpha), levels = alpha_labels),
      rules_raw = factor(rules_raw, levels = rule_order)
    )

  panel_bg <- stats_df %>%
    mutate(
      alpha_lab = factor(sprintf("%.2f", alpha), levels = alpha_labels),
      rules_raw = factor(rules_raw, levels = rule_order)
    )

  ymax <- max(dens_df$density, na.rm = TRUE)

  # Build faceted plot
  p_grid <- ggplot() +
    geom_rect(
      data = panel_bg,
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = mean_p),
      inherit.aes = FALSE, alpha = 0.35
    ) +
    geom_line(
      data = dens_df,
      aes(x = x, y = density,
          group = interaction(alpha_raw, rules_raw)),
      linewidth = 0.4
    ) +
    facet_grid(
      rows = vars(rules_raw),
      cols = vars(alpha_lab),
      switch = "y",
      labeller = labeller(alpha_lab = function(x) paste0("\u03B1 = ", x))
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, ymax), expand = FALSE) +
    scale_fill_viridis_c(name = "mean(p)", option = "viridis") +
    labs(title = sprintf("Final offer distribution p (N=%d)", N_target)) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 8),
      legend.position = "right"
    )

  ggsave(out_plot, p_grid, width = 12, height = 8, dpi = 200)
  message("Wrote plot: ", out_plot)
  message("Wrote stats: ", out_csv)

  invisible(list(plot = p_grid, stats = stats_df, out_plot = out_plot, out_csv = out_csv))
}

# =============================================================================
# SECTION 8 — GRID SIMULATION (configurable, runs when script is executed)
# =============================================================================

# This section sweeps over combinations of N, alpha, and rule subsets.
# Edit the grid values below to control what is run.
run_simulation_grid <- function(
  N_values      = c(100L),
  alpha_values  = c(0, 0.25, 0.5, 0.75, 1),
  rule_subsets  = list(
    c(0L),
    c(3L),
    c(0L, 2L),
    c(0L, 3L),
    c(1L, 2L),
    c(0L, 1L, 2L),
    c(0L, 1L, 2L, 3L)
  ),
  n_instances   = 100,
  n_generations = 20000,
  snapshot_gens = as.integer(c(1, 2, 5, 10, 100, 1000, 10000, 20000)),
  base_dir      = "plots_grid",
  render_plots  = FALSE
) {
  dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

  rule_label <- function(rs) paste0("rules_", paste(rs, collapse = "-"))

  run_log <- list()
  k <- 0

  for (N in N_values) {
    for (alpha in alpha_values) {
      for (rs in rule_subsets) {
        k <- k + 1
        label <- paste0(
          "N_", N,
          "_alpha_", gsub("\\.", "p", as.character(alpha)),
          "_", rule_label(rs)
        )

        message("\n[", k, "] Running: ", label)

        # Run simulation
        res <- run_coevo_pipeline(
          n_instances   = n_instances,
          N             = N,
          alpha         = alpha,
          n_generations = n_generations,
          snapshot_gens = snapshot_gens,
          rule_subset   = rs
        )

        # Save plots into an isolated folder per configuration
        outdir <- file.path(base_dir, label)
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        plot_out <- plot_all(res, base_dir = outdir, prefix = label,
                             render = render_plots)
        cat("Saved plots to:", plot_out$outdir, "\n")

        run_log[[label]] <- list(
          N = N, alpha = alpha,
          rule_subset = rs, plot_dir = plot_out$outdir
        )
      }
    }
  }

  # Build log data.frame for quick viewing
  log_df <- do.call(rbind, lapply(names(run_log), function(nm) {
    x <- run_log[[nm]]
    data.frame(
      run = nm,
      N = x$N,
      alpha = x$alpha,
      rule_subset = paste(x$rule_subset, collapse = ","),
      plot_dir = x$plot_dir,
      stringsAsFactors = FALSE
    )
  }))

  cat("\n=== Simulation grid complete ===\n")
  print(log_df)

  invisible(log_df)
}

# =============================================================================
# SECTION 9 — MAIN: run when invoked as a script
# =============================================================================

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if ("--help" %in% args || "-h" %in% args) {
    cat(
      "Usage: Rscript game_theory_on_networks.R [OPTIONS]\n",
      "\n",
      "Options:\n",
      "  --help, -h       Show this help message and exit\n",
      "  --grid           Run the full simulation grid (default)\n",
      "  --summary        Summarize existing plots_grid/ folder structure\n",
      "  --kde-grid       Build faceted KDE grid plot from saved results\n",
      "\n",
      "By default (no arguments), the full simulation grid is run.\n",
      "Edit N_values, alpha_values, rule_subsets in\n",
      "run_simulation_grid() to customize the parameter sweep.\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  if ("--summary" %in% args) {
    summarize_plots_grid("plots_grid")
    return(invisible(NULL))
  }

  if ("--kde-grid" %in% args) {
    plot_grid_final_p()
    return(invisible(NULL))
  }

  # Default: run the full grid simulation
  run_simulation_grid()
}

# Only run main() when this file is executed as a script (not sourced)
if (!interactive() && identical(sys.frame(1)$ofile, NULL) ||
    (!interactive() && !is.null(sys.frame(1)$ofile))) {
  main()
}
