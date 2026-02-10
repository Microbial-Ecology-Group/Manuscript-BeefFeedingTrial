# ============================================================
# Negative Binomial Loop — Microbiome Version (glmmTMB)
# ============================================================

library(phyloseq)
library(glmmTMB)
library(dplyr)

# Subset to fecal samples only (if needed)
traj.ps <- subset_samples(data_micro.ps, sample_type != "Meat Rinsate")
traj.ps <- prune_taxa(taxa_sums(traj.ps) > 0, traj.ps)

# Metadata
metadata <- as(sample_data(traj.ps), "data.frame")
metadata$SampleID <- rownames(metadata)

# Aggregate to taxonomic level (pick one)
# If you want ASVs: remove this line
group_traj.ps <- tax_glom(traj.ps, taxrank = "Genus")
group_traj.ps <- prune_taxa(taxa_sums(group_traj.ps) > 0, group_traj.ps)

# Extract count data and merge with metadata
counts <- as.data.frame(t(otu_table(group_traj.ps)))
counts$SampleID <- rownames(counts)
full_df <- merge(counts, metadata, by = "SampleID")

# Add modeling variables (same as HFT)
full_df <- full_df[full_df$pool_timepoint != "G", ]
full_df$col_time <- factor(full_df$pool_timepoint, levels = c("A", "B", "C", "D", "E", "F"))
full_df$treatment <- factor(full_df$trt, levels = c("No Claim", "RWA"))
full_df$treatment_order <- factor(full_df$treatment_order, levels = c("AB", "BA"))

# Quick diagnostic summary
nonzero_features <- colSums(counts > 0)
cat("\n### Nonzero count summary across microbiome features ###\n")
print(summary(nonzero_features))
cat("Number of features with >0 counts:", sum(nonzero_features > 0), "\n\n")

# Initialize storage
results_list <- list()
features <- colnames(counts)[1:(ncol(counts)-1)]  # exclude SampleID

# ============================================================
# Main loop (glmmTMB)
# ============================================================
for (feature in features) {
  cat("Trying:", feature, "\n")
  
  # Subset this feature + predictors
  df_tmp <- full_df[, c(feature, "col_time", "treatment", "treatment_order", "participant_id")]
  colnames(df_tmp)[1] <- "y"
  df_tmp$y <- as.numeric(df_tmp$y)
  
  # Diagnostics
  cat("  Nonzero count check:", sum(df_tmp$y > 0, na.rm = TRUE), "nonzero samples\n")
  
  # Skip unusable taxa
  if (any(is.na(df_tmp$y)) || all(df_tmp$y == 0)) {
    cat("❌ Skipped:", feature, "(all zero or NA)\n")
    results_list[[feature]] <- data.frame(
      Feature = feature, AIC_full = NA, AIC_reduced = NA,
      BIC_full = NA, BIC_reduced = NA, DF = NA,
      LRT_Stat = NA, LRT_pvalue = NA
    )
    next
  }
  
  # Fit models
  suppressWarnings({
    full_model <- try(
      glmmTMB(
        y ~ col_time * treatment * treatment_order + (1 | participant_id),
        data = df_tmp,
        family = nbinom2()
      ),
      silent = TRUE
    )
    
    reduced_model <- try(
      glmmTMB(
        y ~ col_time * treatment + (1 | participant_id),
        data = df_tmp,
        family = nbinom2()
      ),
      silent = TRUE
    )
  })
  
  # Handle failures
  if (inherits(full_model, "try-error") || inherits(reduced_model, "try-error")) {
    cat("❌ Model failed:", feature, "\n")
    results_list[[feature]] <- data.frame(
      Feature = feature, AIC_full = NA, AIC_reduced = NA,
      BIC_full = NA, BIC_reduced = NA, DF = NA,
      LRT_Stat = NA, LRT_pvalue = NA
    )
    next
  }
  
  # Extract AIC/BIC
  aic_full <- tryCatch(AIC(full_model), error = function(e) NA)
  aic_reduced <- tryCatch(AIC(reduced_model), error = function(e) NA)
  bic_full <- tryCatch(BIC(full_model), error = function(e) NA)
  bic_reduced <- tryCatch(BIC(reduced_model), error = function(e) NA)
  
  # Likelihood ratio test
  lrt <- suppressWarnings(try(anova(reduced_model, full_model), silent = TRUE))
  
  if (inherits(lrt, "try-error") || is.null(lrt) || nrow(lrt) < 2) {
    p_val <- df_val <- lrt_stat <- NA
  } else {
    p_val <- lrt$`Pr(>Chisq)`[2]
    df_val <- lrt$Df[2]
    lrt_stat <- lrt$Chisq[2]
  }
  
  # Store results
  results_list[[feature]] <- data.frame(
    Feature      = feature,
    AIC_full     = aic_full,
    AIC_reduced  = aic_reduced,
    BIC_full     = bic_full,
    BIC_reduced  = bic_reduced,
    DF           = df_val,
    LRT_Stat     = lrt_stat,
    LRT_pvalue   = p_val
  )
  cat("✔ SUCCESS (recorded):", feature, "\n")
}

# Combine results
results_df <- bind_rows(results_list)

# ============================================================
# Clean models + BH/BY
# ============================================================

results_df$BH <- p.adjust(results_df$LRT_pvalue, method = "BH")
results_df$BY <- p.adjust(results_df$LRT_pvalue, method = "BY")

clean_models <- subset(
  results_df,
  !is.na(AIC_full) &
    !is.na(AIC_reduced) &
    is.finite(AIC_full) &
    is.finite(AIC_reduced) &
    abs(AIC_full) < 1e5 &
    abs(AIC_reduced) < 1e5
)

clean_models$delta_AIC <- clean_models$AIC_full - clean_models$AIC_reduced

# Save clean results
write.csv(clean_models, "Microbiome_NB_Clean_Results.csv", row.names = FALSE)
cat("\n✅ Clean model file created:", nrow(clean_models), "valid models.\n")

# ============================================================
# MICROBIOME NB MODEL SUMMARY — ONE BLOCK
# ============================================================

# Total features tested
total_features <- nrow(results_df)

# Successful fits (AIC_full not NA)
successful <- sum(!is.na(results_df$AIC_full))

# Clean valid models
clean_n <- nrow(clean_models)

# ΔAIC categories
clean_models$abs_delta <- abs(clean_models$delta_AIC)

n_negligible <- sum(clean_models$abs_delta < 2, na.rm = TRUE)
n_small      <- sum(clean_models$abs_delta >= 2 & clean_models$abs_delta < 10, na.rm = TRUE)
n_moderate   <- sum(clean_models$abs_delta >= 10 & clean_models$abs_delta < 100, na.rm = TRUE)
n_extreme    <- sum(clean_models$abs_delta >= 100, na.rm = TRUE)

# Print summary
cat("\n### MICROBIOME NEGATIVE BINOMIAL SUMMARY ###\n")
cat("Total features tested:", total_features, "\n")
cat("Models successfully fit:", successful, "\n")
cat("Clean models with valid AIC:", clean_n, "\n\n")

cat("ΔAIC < 2 (negligible):", n_negligible, "\n")
cat("ΔAIC 2–10 (small):", n_small, "\n")
cat("ΔAIC 10–100 (moderate):", n_moderate, "\n")
cat("ΔAIC > 100 (extreme):", n_extreme, "\n")

