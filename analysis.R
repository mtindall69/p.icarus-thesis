## Analysis of female Polyommatus icarus wing blueness:
## Genetic variation vs thermal plasticity (GxE interaction).
##
## R port of analysis.py — produces equivalent console output and plots.


suppressPackageStartupMessages({
  library(lme4)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra)
})

setwd("C:/Users/maddi/Documents/LU CLASS OF 2026/thesis")
rm(list = ls())

PLOT_DIR <- file.path("plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

PAL_TEMP <- c("Cold (18°C)" = "#3B7DD8", "Warm (26°C)" = "#E8712A")
PAL_REGION  <- c("Öland" = "#2CA02C", "Skåne" = "#9467BD")
TEMP_LABELS <- c("Cold (18°C)", "Warm (26°C)")
PAL_SEX <- c("Female" = "#D81B60", "Male" = "#1C05B3")


# ── Helpers ───────────────────────────────────────────────────────────────────

section_header <- function(title) {
  rule <- strrep("=", 70)
  cat(rule, "\n")
  cat(title, "\n")
  cat(rule, "\n")
}

## Print a Region x Temp crosstab from ind, counting distinct values of id_col
print_crosstab <- function(ind, id_col, label) {
  ct <- ind %>%
    group_by(region_label, temp_label) %>%
    summarise(n = n_distinct(.data[[id_col]]), .groups = "drop") %>%
    pivot_wider(names_from = temp_label, values_from = n, values_fill = 0) %>%
    as.data.frame()
  cat(label, "\n")
  rownames(ct) <- ct$region_label; ct$region_label <- NULL
  ct <- ct[c("Öland", "Skåne"), c("Cold (18°C)", "Warm (26°C)"), drop = FALSE]
  print(ct)
  cat("\n")
}

## Summarise a posterior sample: print mean, 95% CI, and directional probability
print_posterior <- function(samples, name) {
  m  <- mean(samples)
  ci <- quantile(samples, c(0.025, 0.975))
  p_dir <- min(mean(samples > 0), mean(samples < 0)) * 2
  cat(sprintf("    %-12s: mean=%+.4f  95%% CI=[%+.4f, %+.4f]  P(direction)=%.4f\n",
              name, m, ci[1], ci[2], 1 - p_dir / 2))
}

## Format mean [lo, hi] for a posterior sample (used in Tobit output)
fmt_ci <- function(samples) {
  m  <- mean(samples)
  ci <- quantile(samples, c(0.025, 0.975))
  list(m = m, lo = ci[1], hi = ci[2])
}

## Draw one multivariate normal sample via Cholesky decomposition
mvrnorm_chol <- function(mu, sigma) {
  L <- chol(sigma)
  z <- rnorm(length(mu))
  as.numeric(mu) + as.numeric(t(L) %*% z)
}


# ── 1. Load & preprocess ────────────────────────────────────────────────────

load_data <- function() {
  df <- read.csv("blueness.csv", stringsAsFactors = FALSE)
  df$temp_label <- ifelse(df$temp == 18, "Cold (18°C)", "Warm (26°C)")
  df$region_label  <- ifelse(df$region == "O", "Öland", "Skåne")
  df$sex_label <- ifelse(df$sex == "F", "Female", "Male") 
  df$motherID   <- as.factor(df$motherID)
  df$region = as.factor(df$region)
  df$temp = as.factor(df$temp)
  df$sex = as.factor(df$sex)
  df
}

# aggregate_to_individual <- function(df) {
#   group_cols <- c("offspringID", "motherID", "region", "region_label",
#                   "temp", "temp_label", "sex", "sex_label", "motherscore", 
#                   "daughterscore", "pupation_weight_g", "adult_weight_g", 
#                   "start_day", "pupa_day", "adult_day", "pupation_length")
#   measure_cols <- c("prop_blue", "total_mm", "blue_mm", "confidence_seg")
# 
#   agg <- df %>%
#     group_by(across(all_of(group_cols))) %>%
#     summarise(across(all_of(measure_cols), mean), .groups = "drop") %>%
#     as.data.frame()
# 
#   grand_mean <- mean(agg$proportion_blue)
#   grand_sd   <- sd(agg$proportion_blue)
#   agg$z_blue <- (agg$proportion_blue - grand_mean) / grand_sd
# 
#   cat(sprintf("Z-transformation: grand mean = %.4f, grand SD = %.4f\n",
#               grand_mean, grand_sd))
#   cat(sprintf("  z_blue range: [%.2f, %.2f]\n\n",
#               min(agg$z_blue), max(agg$z_blue)))
#   agg
# }

aggregate_to_individual <- function(df) {
  group_cols <- c("offspringID", "motherID", "region", "region_label",
                  "temp", "temp_label", "sex", "sex_label", "motherscore", 
                  "daughterscore", "pupation_weight_g", "adult_weight_g", 
                  "start_day", "pupa_day", "adult_day", "pupation_length")
  measure_cols <- c("total_mm", "blue_mm", "prop_blue")
  
  # Overall averages (all 4 wings)
  agg_all <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(across(all_of(measure_cols), mean), .groups = "drop")
  
  # Forewing averages (FL and FR only)
  agg_fw <- df %>%
    filter(wing %in% c("FL", "FR")) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      avg_fw_total_mm  = mean(total_mm,  na.rm = TRUE),
      avg_fw_blue_mm   = mean(blue_mm,   na.rm = TRUE),
      avg_fw_prop_blue = mean(prop_blue, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Hindwing averages (HL and HR only)
  agg_hw <- df %>%
    filter(wing %in% c("HL", "HR")) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      avg_hw_total_mm  = mean(total_mm,  na.rm = TRUE),
      avg_hw_blue_mm   = mean(blue_mm,   na.rm = TRUE),
      avg_hw_prop_blue = mean(prop_blue, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Join all three together
  agg <- agg_all %>%
    left_join(agg_fw, by = group_cols) %>%
    left_join(agg_hw, by = group_cols) %>%
    rename(
      avg_total_mm  = total_mm,
      avg_blue_mm   = blue_mm,
      avg_prop_blue = prop_blue
    ) %>%
    as.data.frame()
  
  agg
}

# ── 2. Sample summary ───────────────────────────────────────────────────────

print_sample_summary <- function(ind) {
  section_header("SAMPLE SUMMARY")

  cat(sprintf("Total mothers: %d\n", length(unique(ind$motherID))))
  cat(sprintf("Total offspring (individuals): %d\n\n", length(unique(ind$offspringID))))

  print_crosstab(ind, "offspringID", "Offspring per Region × Temperature:")
  print_crosstab(ind, "motherID", "Mothers represented per Region × Temperature:")

  cat("Proportional blueness summary by group:\n")
  summ <- ind %>%
    group_by(region_label, temp_label) %>%
    summarise(N = n(), Mean = mean(avg_prop_blue), SD = sd(avg_prop_blue),
              Median = median(avg_prop_blue), .groups = "drop") %>%
    arrange(region_label, temp_label) %>%
    as.data.frame()
  for (i in seq_len(nrow(summ))) {
    cat(sprintf("%-6s %-14s  %3d  %7.4f  %6.4f  %7.4f\n",
                summ$region_label[i], summ$temp_label[i],
                summ$N[i], summ$Mean[i], summ$SD[i], summ$Median[i]))
  }
  cat("\n")
}


# ── 3. Mixed models ─────────────────────────────────────────────────────────

run_mixed_models <- function(ind) {
  section_header("LINEAR MIXED MODEL: avg_blue_mm ~ avg_total_mm + temp * region + (1|motherID)")

  ind$temp_f <- relevel(factor(ind$temp), ref = "26")
  ind$region_f  <- relevel(factor(ind$region), ref = "S")

  m1 <- lmer(avg_blue_mm ~ avg_total_mm + temp_f * region_f + (1 | motherID), data = ind, REML = TRUE)
  cat("\nRandom-intercept model summary:\n")
  print(summary(m1))
  cat("\n")

  vc <- as.data.frame(VarCorr(m1))
  var_mother <- vc$vcov[vc$grp == "motherID"]
  var_resid  <- vc$vcov[vc$grp == "Residual"]
  icc <- var_mother / (var_mother + var_resid)
  cat(sprintf("Among-mother variance:  %.4f\n", var_mother))
  cat(sprintf("Residual variance:      %.4f\n", var_resid))
  cat(sprintf("ICC (mother):           %.4f  (%.1f%% of variance due to mother)\n\n",
              icc, icc * 100))

  section_header("GxE TEST: random slope for temperature")

  tryCatch({
    m2 <- lmer(avg_blue_mm ~ avg_total_mm + temp_f * region_f + (1 + temp_f | motherID),
               data = ind, REML = TRUE)
    cat("\nRandom-slope model summary:\n")
    print(summary(m2))
    cat("\n")

    m1_ml <- lmer(avg_blue_mm ~ avg_total_mm + temp_f * region_f + (1 | motherID),
                  data = ind, REML = FALSE)
    m2_ml <- lmer(avg_blue_mm ~ avg_total_mm + temp_f * region_f + (1 + temp_f | motherID),
                  data = ind, REML = FALSE)

    lr_stat <- as.numeric(-2 * (logLik(m1_ml) - logLik(m2_ml)))
    df_diff <- 2
    p_val <- pchisq(lr_stat, df_diff, lower.tail = FALSE)
    cat(sprintf("Likelihood ratio test: χ²=%.4f, df=%d, p=%.4g\n",
                lr_stat, df_diff, p_val))
    if (p_val < 0.05) {
      cat("→ Significant GxE: families differ in their response to temperature\n")
    } else {
      cat("→ No significant GxE detected at α=0.05\n")
    }
    cat("\n")
  }, error = function(e) {
    cat(sprintf("Random slope model failed: %s\n\n", e$message))
  })

  invisible(m1)
}


# ── 4. Bayesian mother-daughter regression (Gibbs sampler) ───────────────────

## Unified Gibbs sampler for the mother-daughter random-intercept model.
## When censor_thresh is non-NULL, applies Tobit data augmentation for
## left-censored observations (proportion_blue <= censor_thresh).
gibbs_mother_daughter <- function(y, motherscore, mother_idx,
                                  censor_thresh = NULL,
                                  n_iter = 15000, burnin = 5000,
                                  thin = 2, seed = 42) {
  set.seed(seed)

  n <- length(y)
  unique_mothers <- sort(unique(mother_idx))
  J <- length(unique_mothers)
  mom_map <- setNames(seq_along(unique_mothers), unique_mothers)
  mom_int <- mom_map[as.character(mother_idx)]

  X <- cbind(1, motherscore)
  p <- ncol(X)

  # Censoring setup (Tobit only)
  is_censored <- if (!is.null(censor_thresh)) y <= censor_thresh else rep(FALSE, n)
  n_censored  <- sum(is_censored)

  # Prior hyperparams
  beta_prior_mean <- rep(0, p)
  beta_prior_prec <- diag(p) / 100.0
  a_e <- 0.001; b_e <- 0.001
  a_u <- 0.001; b_u <- 0.001

  # Initialize
  beta   <- rep(0, p)
  u      <- rep(0, J)
  sig2_e <- var(y) * 0.5
  sig2_u <- 0.1
  y_star <- y

  # Storage
  n_keep <- (n_iter - burnin) %/% thin
  beta_samples   <- matrix(0, n_keep, p)
  sig2_e_samples <- numeric(n_keep)
  sig2_u_samples <- numeric(n_keep)
  y_star_cens_sum <- numeric(max(n_censored, 0))

  idx <- 0
  for (it in seq_len(n_iter)) {
    # 0) DATA AUGMENTATION (Tobit only) -- impute y* for censored obs
    if (n_censored > 0) {
      mu <- as.numeric(X %*% beta + u[mom_int])
      sd_e <- sqrt(sig2_e)
      mu_cens <- mu[is_censored]
      Phi_b <- pmin(pmax(pnorm(censor_thresh, mean = mu_cens, sd = sd_e), 1e-12), 1.0)
      u_vals <- pmin(pmax(runif(n_censored, 0, Phi_b), 1e-15), 1.0 - 1e-15)
      y_star[is_censored] <- qnorm(u_vals, mean = mu_cens, sd = sd_e)
    }

    # 1) Sample beta | rest
    y_adj <- y_star - u[mom_int]
    prec_post <- crossprod(X) / sig2_e + beta_prior_prec
    cov_post  <- solve(prec_post)
    mean_post <- cov_post %*% (crossprod(X, y_adj) / sig2_e +
                                beta_prior_prec %*% beta_prior_mean)
    beta <- as.numeric(mvrnorm_chol(mean_post, cov_post))

    # 2) Sample u_j | rest
    resid <- y_star - X %*% beta
    for (j in seq_len(J)) {
      mask <- mom_int == j
      n_j  <- sum(mask)
      if (n_j == 0) {
        u[j] <- rnorm(1, 0, sqrt(sig2_u))
        next
      }
      prec_j <- n_j / sig2_e + 1.0 / sig2_u
      var_j  <- 1.0 / prec_j
      mean_j <- var_j * (sum(resid[mask]) / sig2_e)
      u[j]   <- rnorm(1, mean_j, sqrt(var_j))
    }

    # 3) Sample sig2_e | rest
    full_resid <- y_star - X %*% beta - u[mom_int]
    shape_e <- a_e + n / 2.0
    rate_e  <- b_e + sum(full_resid^2) / 2.0
    sig2_e  <- 1.0 / rgamma(1, shape_e, rate = rate_e)

    # 4) Sample sig2_u | rest
    shape_u <- a_u + J / 2.0
    rate_u  <- b_u + sum(u^2) / 2.0
    sig2_u  <- 1.0 / rgamma(1, shape_u, rate = rate_u)

    # Store
    if (it > burnin && (it - burnin) %% thin == 0) {
      idx <- idx + 1
      beta_samples[idx, ] <- beta
      sig2_e_samples[idx] <- sig2_e
      sig2_u_samples[idx] <- sig2_u
      if (n_censored > 0) {
        y_star_cens_sum <- y_star_cens_sum + y_star[is_censored]
      }
    }
  }

  y_star_cens_mean <- if (n_censored > 0) y_star_cens_sum / n_keep else numeric(0)

  list(beta = beta_samples, sig2_e = sig2_e_samples,
       sig2_u = sig2_u_samples, unique_mothers = unique_mothers,
       is_censored = is_censored, n_censored = n_censored,
       y_star_cens_mean = y_star_cens_mean)
}

run_bayesian_mother_daughter <- function(ind) {
  section_header("BAYESIAN MOTHER-DAUGHTER REGRESSION (individual daughters, Gibbs)")
  cat("Model: avg_blue_mm_ij = α + β·motherscore_j + u_j + ε_ij\n")
  cat("       u_j ~ N(0, σ²_u),  ε_ij ~ N(0, σ²_e)\n")
  cat("       Weakly informative priors, 15000 iterations, 5000 burn-in\n\n")

  results <- list()
  for (temp_label in TEMP_LABELS) {
    sub <- ind[ind$temp_label == temp_label, ]
    res <- gibbs_mother_daughter(sub$avg_prop_blue, as.numeric(sub$motherscore), sub$motherID)

    alpha_post  <- res$beta[, 1]
    beta_post   <- res$beta[, 2]
    sig2_u_post <- res$sig2_u
    sig2_e_post <- res$sig2_e

    cat(sprintf("  %s  (N=%d daughters, %d families)\n",
                temp_label, nrow(sub), length(unique(sub$motherID))))
    print_posterior(alpha_post,  "Intercept")
    print_posterior(beta_post,   "Slope (β)")
    print_posterior(sig2_u_post, "σ²_mother")
    print_posterior(sig2_e_post, "σ²_residual")

    icc_post <- sig2_u_post / (sig2_u_post + sig2_e_post)
    ci_icc <- quantile(icc_post, c(0.025, 0.975))
    cat(sprintf("    %-12s: mean=%.4f  95%% CI=[%.4f, %.4f]\n\n",
                "ICC", mean(icc_post), ci_icc[1], ci_icc[2]))

    results[[temp_label]] <- list(
      beta_samples   = beta_post,
      alpha_samples  = alpha_post,
      sig2_u_samples = sig2_u_post,
      sig2_e_samples = sig2_e_post,
      data = sub
    )
  }

  b_cold <- results[["Cold (18°C)"]]$beta_samples
  b_warm <- results[["Warm (26°C)"]]$beta_samples
  min_len <- min(length(b_cold), length(b_warm))
  diff <- b_cold[1:min_len] - b_warm[1:min_len]
  ci_diff <- quantile(diff, c(0.025, 0.975))
  p_diff  <- min(mean(diff > 0), mean(diff < 0)) * 2
  cat(sprintf("  Slope difference (Cold − Warm): mean=%+.4f, 95%% CI=[%+.4f, %+.4f], P(overlap 0)=%.4f\n",
              mean(diff), ci_diff[1], ci_diff[2], p_diff))
  if (p_diff < 0.05) {
    cat("  → Mother-daughter relationship strength differs by temperature\n")
  } else {
    cat("  → No clear evidence that heritability differs between temperatures\n")
  }
  cat("\n")

  results
}


# ── 4b. Tobit (censored) Bayesian mother-daughter regression ─────────────────

run_tobit_mother_daughter <- function(ind) {
  section_header("TOBIT (CENSORED) MOTHER-DAUGHTER REGRESSION — floor-effect correction")
  cat("Observations with proportion_blue ≤ 0.005 are left-censored.\n")
  cat("Latent values below the floor are imputed via data augmentation.\n\n")

  for (temp_label in TEMP_LABELS) {
    sub <- ind[ind$temp_label == temp_label, ]
    n_cens <- sum(sub$avg_prop_blue <= 0.005)
    cat(sprintf("  %s: %d/%d censored at floor (%.1f%%)\n",
                temp_label, n_cens, nrow(sub), n_cens / nrow(sub) * 100))
  }
  cat("\n")

  results <- list()
  for (temp_label in TEMP_LABELS) {
    sub <- ind[ind$temp_label == temp_label, ]
    res <- gibbs_mother_daughter(sub$avg_prop_blue, as.numeric(sub$motherscore),
                                 sub$motherID, censor_thresh = 0.005)

    alpha_post  <- res$beta[, 1]
    beta_post   <- res$beta[, 2]
    sig2_u_post <- res$sig2_u
    sig2_e_post <- res$sig2_e
    icc_post    <- sig2_u_post / (sig2_u_post + sig2_e_post)

    a <- fmt_ci(alpha_post)
    b <- fmt_ci(beta_post)
    su <- fmt_ci(sig2_u_post)
    se <- fmt_ci(sig2_e_post)
    ic <- fmt_ci(icc_post)

    cat(sprintf("  %s  (N=%d, %d censored, %d families)\n",
                temp_label, nrow(sub), res$n_censored, length(unique(sub$motherID))))
    cat(sprintf("    Intercept (α): %+.4f  [%+.4f, %+.4f]\n", a$m, a$lo, a$hi))
    cat(sprintf("    Slope (β):     %+.4f  [%+.4f, %+.4f]\n", b$m, b$lo, b$hi))
    cat(sprintf("    σ²_mother:     %.4f  [%.4f, %.4f]\n", su$m, su$lo, su$hi))
    cat(sprintf("    σ²_residual:   %.4f  [%.4f, %.4f]\n", se$m, se$lo, se$hi))
    cat(sprintf("    ICC:           %.4f  [%.4f, %.4f]\n\n", ic$m, ic$lo, ic$hi))

    results[[temp_label]] <- list(
      beta_samples   = beta_post,
      alpha_samples  = alpha_post,
      sig2_u_samples = sig2_u_post,
      sig2_e_samples = sig2_e_post,
      icc_samples    = icc_post,
      is_censored    = res$is_censored,
      n_censored     = res$n_censored,
      y_star_cens_mean = res$y_star_cens_mean,
      data = sub
    )
  }

  results
}


# ── 5. Variance decomposition by temperature ────────────────────────────────

variance_decomposition_by_temp <- function(ind) {
  section_header("VARIANCE DECOMPOSITION BY TEMPERATURE")

  results <- list()
  for (temp_label in TEMP_LABELS) {
    subset_df <- ind[ind$temp_label == temp_label, ]
    subset_df$region_f <- relevel(factor(subset_df$region), ref = "S")
    m <- lmer(avg_blue_mm ~ avg_total_mm + region_f + (1 | motherID), data = subset_df, REML = TRUE)

    vc <- as.data.frame(VarCorr(m))
    var_mother <- vc$vcov[vc$grp == "motherID"]
    var_resid  <- vc$vcov[vc$grp == "Residual"]
    total <- var_mother + var_resid

    results[[temp_label]] <- list(
      among_mother = var_mother,
      within_mother = var_resid,
      total = total,
      pct_mother = var_mother / total * 100
    )

    cat(sprintf("\n  %s:\n", temp_label))
    cat(sprintf("    Among-mother variance:  %.4f (%.1f%%)\n",
                var_mother, var_mother / total * 100))
    cat(sprintf("    Residual variance:      %.4f (%.1f%%)\n",
                var_resid, var_resid / total * 100))
  }
  cat("\n")
  results
}


# ── 6. Figures ───────────────────────────────────────────────────────────────

fig1_temp_region_violin <- function(ind) {
  ind$region_label  <- factor(ind$region_label, levels = c("Öland", "Skåne"))
  ind$temp_label <- factor(ind$temp_label, levels = c("Cold (18°C)", "Warm (26°C)"))

  # Compute means and SEs for diamonds
  summ <- ind %>%
    group_by(region_label, temp_label) %>%
    summarise(m = mean(avg_prop_blue), se = sd(avg_prop_blue) / sqrt(n()), .groups = "drop")

  p <- ggplot(ind, aes(x = region_label, y = avg_prop_blue, fill = temp_label)) +
    geom_violin(alpha = 0.3, position = position_dodge(width = 0.8),
                trim = TRUE, scale = "width", colour = NA) +
    geom_jitter(aes(colour = temp_label),
                position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.12),
                alpha = 0.35, size = 0.8, show.legend = FALSE) +
    geom_pointrange(data = summ,
                    aes(x = region_label, y = m,
                        ymin = m - se, ymax = m + se, group = temp_label),
                    position = position_dodge(width = 0.8),
                    shape = 18, size = 0.6, colour = "black", show.legend = FALSE) +
    geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
    scale_fill_manual(values = PAL_TEMP) +
    scale_colour_manual(values = PAL_TEMP) +
    labs(x = "Region", y = "Proportion blue",
         title = "Wing blueness by region and temperature",
         fill = "Temperature") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.85, 0.9),
          legend.background = element_rect(colour = "grey80"))

  ggsave(file.path(PLOT_DIR, "fig1_temp_region_violin_r.png"), p,
         width = 7, height = 5, dpi = 200)
  cat("Saved: plots/fig1_temp_region_violin_r.png\n")
}


fig2_reaction_norms <- function(ind) {
  fam_means <- ind %>%
    group_by(motherID, region_label, temp_label) %>%
    summarise(mean_blue = mean(avg_prop_blue), n = n_distinct(offspringID), .groups = "drop")

  fam_both <- fam_means %>%
    group_by(motherID) %>%
    filter(n_distinct(temp_label) == 2) %>%
    ungroup()

  grand <- fam_both %>%
    group_by(region_label, temp_label) %>%
    summarise(mean_blue = mean(mean_blue), .groups = "drop")

  p <- ggplot() +
    geom_line(data = fam_both,
              aes(x = temp_label, y = mean_blue, group = motherID, colour = region_label),
              alpha = 0.45, linewidth = 0.5) +
    geom_point(data = fam_both,
               aes(x = temp_label, y = mean_blue, group = motherID, colour = region_label),
               alpha = 0.45, size = 1.5) +
    geom_line(data = grand,
              aes(x = temp_label, y = mean_blue, group = region_label, colour = region_label),
              linewidth = 1.8) +
    geom_point(data = grand,
               aes(x = temp_label, y = mean_blue, colour = region_label),
               size = 4, shape = 15, stroke = 1) +
    geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
    scale_colour_manual(values = PAL_REGION) +
    scale_x_discrete(expand = expansion(mult = 0.1)) +
    labs(x = "Temperature treatment",
         y = "Wing blueness (proportion, family mean)",
         title = "Reaction norms: family-level response to temperature",
         colour = "region") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.85, 0.85),
          legend.background = element_rect(colour = "grey80"))

  ggsave(file.path(PLOT_DIR, "fig2_reaction_norms_r.png"), p,
         width = 6, height = 5, dpi = 200)
  cat("Saved: plots/fig2_reaction_norms_r.png\n")

  n_fam <- length(unique(fam_both$motherID))
  cat(sprintf("  Families with offspring in both temperatures: %d\n", n_fam))
}


fig3_bayesian_mother_daughter <- function(bayes_results) {
  plots <- list()

  for (col in seq_along(TEMP_LABELS)) {
    temp_label <- TEMP_LABELS[col]
    res <- bayes_results[[temp_label]]
    sub <- res$data
    alpha_post <- res$alpha_samples
    beta_post  <- res$beta_samples

    # Posterior regression line
    x_grid <- seq(1.5, 5.5, length.out = 200)
    n_draw <- min(500, length(alpha_post))
    set.seed(123)
    draw_idx <- sample(length(alpha_post), n_draw, replace = FALSE)
    y_draws <- outer(alpha_post[draw_idx], x_grid, function(a, x) a) +
               outer(beta_post[draw_idx], x_grid, function(b, x) b * x)
    y_mean <- colMeans(y_draws)
    y_lo <- apply(y_draws, 2, quantile, 0.025)
    y_hi <- apply(y_draws, 2, quantile, 0.975)

    line_df <- data.frame(x = x_grid, y = y_mean, lo = y_lo, hi = y_hi)

    beta_mean <- mean(beta_post)
    beta_ci   <- quantile(beta_post, c(0.025, 0.975))
    ann_label <- sprintf("β = %.3f\n95%% CI [%.3f, %.3f]",
                         beta_mean, beta_ci[1], beta_ci[2])

    p_scatter <- ggplot(sub, aes(x = motherscore, y = avg_prop_blue)) +
      geom_point(aes(colour = region_label), alpha = 0.4, size = 1.2) +
      geom_ribbon(data = line_df, aes(x = x, ymin = lo, ymax = hi),
                  inherit.aes = FALSE, fill = "gray", alpha = 0.25) +
      geom_line(data = line_df, aes(x = x, y = y),
                inherit.aes = FALSE, colour = "black", linewidth = 1) +
      geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
      annotate("label", x = 1.7, y = max(sub$avg_prop_blue) * 0.95,
               label = ann_label, hjust = 0, vjust = 1, size = 3,
               fill = "white", label.size = 0.3) +
      scale_colour_manual(values = PAL_REGION) +
      scale_x_continuous(breaks = 2:5) +
      labs(x = "Mother blueness score",
           title = temp_label,
           colour = "region") +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, colour = PAL_TEMP[temp_label]))

    if (col == 1) {
      p_scatter <- p_scatter + labs(y = "Daughter wing blueness (%)") +
        theme(legend.position = "inside",
              legend.position.inside = c(0.8, 0.15),
              legend.background = element_rect(colour = "grey80"),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 9))
    } else {
      p_scatter <- p_scatter + labs(y = "") +
        theme(legend.position = "none")
    }

    # Density plot
    dens_df <- data.frame(beta = beta_post)
    p_density <- ggplot(dens_df, aes(x = beta)) +
      geom_histogram(aes(y = after_stat(density)), bins = 50,
                     fill = PAL_TEMP[temp_label], alpha = 0.6,
                     colour = "white", linewidth = 0.2) +
      geom_density(colour = "black", linewidth = 0.8) +
      geom_vline(xintercept = 0, colour = "gray", linetype = "dashed", linewidth = 0.7) +
      geom_vline(xintercept = beta_mean, colour = "red", linewidth = 0.8, alpha = 0.8) +
      labs(x = "Slope β (mother → daughter)", title = "Posterior of β") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))

    if (col == 1) {
      p_density <- p_density + labs(y = "Posterior density")
    } else {
      p_density <- p_density + labs(y = "")
    }

    plots[[paste0("scatter_", col)]] <- p_scatter
    plots[[paste0("density_", col)]] <- p_density
  }

  combined <- arrangeGrob(
    plots$scatter_1, plots$scatter_2,
    plots$density_1, plots$density_2,
    ncol = 2, nrow = 2,
    heights = c(2.5, 1),
    top = textGrob("Bayesian mother-daughter regression (all individual daughters)",
                   gp = gpar(fontface = "bold", fontsize = 14))
  )

  combined <- gTree(children = gList(
    rectGrob(gp = gpar(fill = "white", col = NA)),
    combined
  ))
  ggsave(file.path(PLOT_DIR, "fig3_mother_daughter_bayesian_r.png"), combined,
         width = 11, height = 9, dpi = 200)
  cat("Saved: plots/fig3_mother_daughter_bayesian_r.png\n")
}


fig4_developmental_traits <- function(ind) {
  ind$region_label  <- factor(ind$region_label, levels = c("Öland", "Skåne"))
  ind$temp_label <- factor(ind$temp_label, levels = c("Cold (18°C)", "Warm (26°C)"))

  p1 <- ggplot(ind, aes(x = region_label, y = adult_weight_g, fill = temp_label)) +
    geom_boxplot(outlier.size = 1, linewidth = 0.5) +
    scale_fill_manual(values = PAL_TEMP, guide = "none") +
    labs(x = "Region", y = "Adult weight (g)",
         title = "A) Adult body size") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  p2 <- ggplot(ind, aes(x = region_label, y = pupation_length, fill = temp_label)) +
    geom_boxplot(outlier.size = 1, linewidth = 0.5) +
    scale_fill_manual(values = PAL_TEMP) +
    labs(x = "Region", y = "Pupation duration (days)",
         title = "B) Development time", fill = "Temperature") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.85, 0.85),
          legend.background = element_rect(colour = "grey80"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9))

  p3 <- ggplot(ind, aes(x = avg_total_mm, y = avg_blue_mm, colour = temp_label)) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, alpha = 0.8) +
    geom_hline(yintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.5) +
    scale_colour_manual(values = PAL_TEMP) +
    labs(x = "Total wing area (mm²)", y = "Wing blue area (mm²)",
         title = "C) Wing size vs. blueness", colour = "Temperature") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.8, 0.85),
          legend.background = element_rect(colour = "grey80"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9))

  combined <- grid.arrange(p1, p2, p3, ncol = 3)
  ggsave(file.path(PLOT_DIR, "fig4_developmental_traits_r.png"), combined,
         width = 14, height = 4.5, dpi = 200)
  cat("Saved: plots/fig4_developmental_traits_r.png\n")
}


fig5_variance_decomposition <- function(var_results) {
  temps <- sort(names(var_results))
  among  <- sapply(temps, function(t) var_results[[t]]$among_mother)
  within <- sapply(temps, function(t) var_results[[t]]$within_mother)
  pct    <- sapply(temps, function(t) var_results[[t]]$pct_mother)

  bar_df <- data.frame(
    temp = rep(temps, 2),
    component = rep(c("Among-mother (genetic + maternal)",
                      "Within-mother (residual)"), each = length(temps)),
    value = c(among, within),
    stringsAsFactors = FALSE
  )
  bar_df$component <- factor(bar_df$component,
                             levels = c("Within-mother (residual)",
                                        "Among-mother (genetic + maternal)"))

  label_df <- data.frame(
    temp = rep(temps, 2),
    y = c(among / 2, among + within / 2),
    label = c(sprintf("%.1f%%", pct), sprintf("%.1f%%", 100 - pct)),
    stringsAsFactors = FALSE
  )

  p <- ggplot(bar_df, aes(x = temp, y = value, fill = component)) +
    geom_col(width = 0.5, alpha = 0.85) +
    geom_text(data = label_df, aes(x = temp, y = y, label = label),
              inherit.aes = FALSE,
              colour = "white", fontface = "bold", size = 4) +
    scale_fill_manual(values = c("Within-mother (residual)" = "#AAAAAA",
                                 "Among-mother (genetic + maternal)" = "#2CA02C")) +
    labs(x = "", y = "Variance in proportion blueness",
         title = "Variance decomposition by temperature",
         fill = "") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.65, 0.85),
          legend.background = element_rect(colour = "grey80"),
          legend.title = element_blank(),
          legend.text = element_text(size = 9))

  ggsave(file.path(PLOT_DIR, "fig5_variance_decomposition_r.png"), p,
         width = 5, height = 5, dpi = 200)
  cat("Saved: plots/fig5_variance_decomposition_r.png\n")
}


fig6_floor_effect_tobit <- function(ind, naive_results, tobit_results) {
  censor_thresh <- 0.005

  # Panel A: Observed distribution
  p_hist <- ggplot(ind, aes(x = avg_prop_blue, fill = temp_label)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, alpha = 0.5, position = "identity",
                   colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = censor_thresh, colour = "red",
               linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values = PAL_TEMP) +
    labs(x = "Proportion blue (observed)", y = "Density",
         title = "A) Floor effect: pile-up at zero", fill = "") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.7, 0.85),
          legend.background = element_rect(colour = "grey80"),
          legend.text = element_text(size = 8))

  # Panel B: Latent distribution from Tobit
  latent_list <- lapply(TEMP_LABELS, function(tl) {
    tres <- tobit_results[[tl]]
    data.frame(
      value = c(tres$data$avg_prop_blue[!tres$is_censored], tres$y_star_cens_mean),
      temp_label = tl
    )
  })
  latent_df <- do.call(rbind, latent_list)

  p_latent <- ggplot(latent_df, aes(x = value, fill = temp_label)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, alpha = 0.5, position = "identity",
                   colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "gray", linetype = "dotted", linewidth = 0.7) +
    geom_vline(xintercept = censor_thresh, colour = "red",
               linetype = "dashed", linewidth = 1, alpha = 0.5) +
    scale_fill_manual(values = PAL_TEMP) +
    labs(x = "Latent blueness (Tobit-imputed)", y = "Density",
         title = "B) Latent scale: censored obs spread below zero", fill = "") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
          legend.position = "inside",
          legend.position.inside = c(0.95, 0.95),
          legend.justification = c(1, 1),
          legend.background = element_rect(colour = "grey80"),
          legend.text = element_text(size = 8))

  # Panel C: Forest plot -- ICC comparison
  forest_data <- list()
  y_pos <- 0
  for (tl in TEMP_LABELS) {
    nr <- naive_results[[tl]]
    naive_icc <- nr$sig2_u_samples / (nr$sig2_u_samples + nr$sig2_e_samples)
    tobit_icc <- tobit_results[[tl]]$icc_samples

    for (info in list(
      list(label = paste0(tl, " — Naive"), samples = naive_icc,
           colour = "#999999", shape = 16),
      list(label = paste0(tl, " — Tobit (corrected)"), samples = tobit_icc,
           colour = PAL_TEMP[tl], shape = 18)
    )) {
      s <- fmt_ci(info$samples)
      forest_data[[length(forest_data) + 1]] <- data.frame(
        label = info$label, m = s$m, lo = s$lo, hi = s$hi,
        colour = info$colour, shape = info$shape, y = y_pos,
        stringsAsFactors = FALSE
      )
      y_pos <- y_pos + 1
    }
    y_pos <- y_pos + 0.7
  }
  fdf <- do.call(rbind, forest_data)
  fdf$label <- factor(fdf$label, levels = rev(fdf$label))

  p_forest <- ggplot(fdf, aes(x = m, y = label)) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3,
                   colour = fdf$colour, linewidth = 1) +
    geom_point(size = 4, colour = fdf$colour, shape = fdf$shape) +
    xlim(0, 0.85) +
    labs(x = "ICC — proportion of variance due to mother\n(with 95% credible interval)",
         y = "",
         title = "C) Floor-effect correction on ICC") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

  combined <- grid.arrange(p_hist, p_latent, p_forest, ncol = 3,
                           widths = c(1, 1, 1.2))
  ggsave(file.path(PLOT_DIR, "fig6_floor_effect_tobit_r.png"), combined,
         width = 14, height = 5, dpi = 200)
  cat("Saved: plots/fig6_floor_effect_tobit_r.png\n")
}


# ── Main ─────────────────────────────────────────────────────────────────────

main <- function() {
  df  <- load_data()
  ind <- aggregate_to_individual(df)
  #subset ind for sexes
  ind <- subset(ind, sex=="F")

  print_sample_summary(ind)
  run_mixed_models(ind)
  naive_results <- run_bayesian_mother_daughter(ind)
  tobit_results <- run_tobit_mother_daughter(ind)
  var_results   <- variance_decomposition_by_temp(ind)

  section_header("GENERATING FIGURES")
  fig1_temp_region_violin(ind)
  fig2_reaction_norms(ind)
  fig3_bayesian_mother_daughter(naive_results)
  fig4_developmental_traits(ind)
  fig5_variance_decomposition(var_results)
  fig6_floor_effect_tobit(ind, naive_results, tobit_results)
  cat("\nDone. All figures saved to plots/\n")
}

main()
