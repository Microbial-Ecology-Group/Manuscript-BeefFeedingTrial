#Negative Binomial using GLMMTMB (AMR)

#MASS does not take into account repeated measures

# ============================================================
# Negative Binomial Loop — Final Complete Diagnostic Version
# ============================================================

library(phyloseq)
library(glmmTMB)
library(dplyr)

# Subset to fecal samples only
traj.ps <- subset_samples(AMR_data.ps, Sample_Type != "Beef Sample")
traj.ps <- prune_taxa(taxa_sums(traj.ps) > 0, traj.ps)

# Metadata
metadata <- as(sample_data(traj.ps), "data.frame")
metadata$SampleID <- rownames(metadata)

# Aggregate to gene group level
group_traj.ps <- tax_glom(traj.ps, taxrank = "group")
group_traj.ps <- prune_taxa(taxa_sums(group_traj.ps) > 0, group_traj.ps)

# Extract count data and merge
counts <- as.data.frame(t(otu_table(group_traj.ps)))
counts$SampleID <- rownames(counts)
full_df <- merge(counts, metadata, by = "SampleID")

# Add modeling variables
full_df <- full_df[full_df$pool_timepoint != "G", ]
full_df$col_time <- factor(full_df$pool_timepoint, levels = c("A", "B", "C", "D", "E", "F"))
full_df$treatment <- factor(full_df$trt_id, levels = c("No Claim", "RWA"))
full_df$treatment_order <- factor(full_df$treatment_order, levels = c("AB", "BA"))

# Quick diagnostic summary before loop
nonzero_features <- colSums(counts > 0)
cat("\n### Nonzero count summary across all samples ###\n")
print(summary(nonzero_features))
cat("Number of features with >0 counts:", sum(nonzero_features > 0), "\n\n")

# Initialize storage
results_list <- list()
features <- colnames(counts)[1:535]

# ============================================================
# Main loop (glmmTMB version)
# ============================================================
for (feature in features) {
  cat("Trying:", feature, "\n")
  
  # Build temporary dataframe for this feature
  df_tmp <- full_df[, c(feature, "col_time", "treatment", "treatment_order", "Participant_ID")]
  colnames(df_tmp)[1] <- "y"
  df_tmp$y <- as.numeric(df_tmp$y)
  
  # Diagnostic printout for counts
  cat("  Nonzero count check:", sum(df_tmp$y > 0, na.rm = TRUE), "nonzero samples\n")
  
  # Skip if all zero or NA
  if (any(is.na(df_tmp$y)) || all(df_tmp$y == 0)) {
    cat("❌ Skipped:", feature, "(all zero or NA)\n")
    results_list[[feature]] <- data.frame(
      Feature = feature, AIC_full = NA, AIC_reduced = NA,
      BIC_full = NA, BIC_reduced = NA, DF = NA,
      LRT_Stat = NA, LRT_pvalue = NA
    )
    next
  }
  
  # Fit models (glmmTMB)
  suppressWarnings({
    full_model <- try(
      glmmTMB(
        y ~ col_time * treatment * treatment_order + (1 | Participant_ID),
        data = df_tmp,
        family = nbinom2()
      ),
      silent = TRUE
    )
    
    reduced_model <- try(
      glmmTMB(
        y ~ col_time * treatment + (1 | Participant_ID),
        data = df_tmp,
        family = nbinom2()
      ),
      silent = TRUE
    )
  })
  
  # Handle failed fits
  if (inherits(full_model, "try-error") || inherits(reduced_model, "try-error")) {
    cat("❌ Model failed:", feature, "\n")
    results_list[[feature]] <- data.frame(
      Feature = feature, AIC_full = NA, AIC_reduced = NA,
      BIC_full = NA, BIC_reduced = NA, DF = NA,
      LRT_Stat = NA, LRT_pvalue = NA
    )
    next
  }
  
  # AIC/BIC extraction (glmmTMB-compatible)
  aic_full <- tryCatch(AIC(full_model), error = function(e) NA)
  aic_reduced <- tryCatch(AIC(reduced_model), error = function(e) NA)
  bic_full <- tryCatch(BIC(full_model), error = function(e) NA)
  bic_reduced <- tryCatch(BIC(reduced_model), error = function(e) NA)
  
  # Likelihood-ratio test (glmmTMB)
  lrt <- suppressWarnings(try(anova(reduced_model, full_model), silent = TRUE))
  
  if (inherits(lrt, "try-error") || is.null(lrt) || nrow(lrt) < 2) {
    p_val <- df_val <- lrt_stat <- NA
  } else {
    # glmmTMB stores results in row 2
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





#### Results ####

#All of this below is a mess. You really just need the summary counts and a clean csv

#Lets see what we got. 
# Combine results
results_df <- bind_rows(results_list)

# Summary counts
total_features <- nrow(results_df)
successful <- sum(!is.na(results_df$AIC_full))
identical_aic <- sum(round(results_df$AIC_full, 4) == round(results_df$AIC_reduced, 4), na.rm = TRUE)
mean_delta_aic <- mean(abs(results_df$AIC_full - results_df$AIC_reduced), na.rm = TRUE)
median_delta_aic <- median(abs(results_df$AIC_full - results_df$AIC_reduced), na.rm = TRUE)

cat("\n### NEGATIVE BINOMIAL MODEL SUMMARY ###\n")
cat("Total gene groups tested:", total_features, "\n")
cat("Models successfully fit:", successful, "\n")
cat("Models with identical AIC (full vs reduced):", identical_aic, "\n")
cat("Mean ΔAIC:", mean_delta_aic, "\n")
cat("Median ΔAIC:", median_delta_aic, "\n")

# Optional: inspect distribution
hist(results_df$AIC_full - results_df$AIC_reduced,
     main = "ΔAIC Distribution (Full - Reduced)",
     xlab = "ΔAIC", breaks = 40)

# ------------------------------------------------------------
# Identify non-identical AIC models
# ------------------------------------------------------------

# Add ΔAIC column for convenience
results_df <- results_df %>%
  mutate(delta_AIC = AIC_full - AIC_reduced)

# Subset only models where ΔAIC != 0 (within rounding tolerance)
diff_models <- results_df %>%
  filter(!is.na(AIC_full), !is.na(AIC_reduced)) %>%
  filter(round(delta_AIC, 4) != 0)

# How many?
nrow(diff_models)

# View the gene groups and their AICs
diff_models_sorted <- diff_models %>%
  arrange(desc(abs(delta_AIC))) %>%
  dplyr::select(Feature, AIC_reduced, AIC_full, delta_AIC, BIC_reduced, BIC_full, LRT_pvalue)

# Preview top results
head(diff_models_sorted, 37) #37 is all

# Optional: save to file for deeper inspection
write.csv(diff_models_sorted, "NB_diff_AIC_models25Nov2025.csv", row.names = FALSE)

# Sort and classify ΔAIC differences
diff_summary <- diff_models %>%
  mutate(
    abs_delta = abs(delta_AIC),
    category = case_when(
      abs_delta < 2 ~ "Negligible (<2)",
      abs_delta >= 2 & abs_delta < 10 ~ "Small (2–10)",
      abs_delta >= 10 & abs_delta < 100 ~ "Moderate (10–100)",
      abs_delta >= 100 ~ "Extreme (>100)",
      TRUE ~ "Other"
    ),
    better_model = case_when(
      AIC_reduced < AIC_full ~ "Reduced",
      AIC_reduced > AIC_full ~ "Full",
      TRUE ~ "Equal"
    )
  )

# Summary counts
cat("\n### ΔAIC CATEGORY SUMMARY ###\n")
print(table(diff_summary$category))
cat("\n### WHICH MODEL FIT BETTER ###\n")
print(table(diff_summary$better_model))

# Optional: view details for just the non-negligible ones
non_negligible <- diff_summary %>% filter(abs_delta >= 2) %>%
  arrange(desc(abs_delta)) %>%
  dplyr::select(Feature, AIC_reduced, AIC_full, delta_AIC, better_model, category, LRT_pvalue)

head(non_negligible, 37)

# ============================================================
# Combine, adjust, and save
# ============================================================
model_comparison_df <- do.call(rbind, results_list)
rownames(model_comparison_df) <- NULL
model_comparison_df$BH <- p.adjust(model_comparison_df$LRT_pvalue, method = "BH")
model_comparison_df$BY <- p.adjust(model_comparison_df$LRT_pvalue, method = "BY")

write.csv(model_comparison_df, "NB_Model_Comparison_Results_11182025.csv", row.names = FALSE)
cat("\n✅ Finished —", nrow(model_comparison_df), "models recorded.\n")

#Add cleaned models- these dont have failed runs. 
clean_models <- subset(model_comparison_df, !is.na(AIC_full) & abs(AIC_full) < 1e5)

clean_models$delta_AIC <- clean_models$AIC_full - clean_models$AIC_reduced

write.csv(clean_models,
          "NB_Model_Comparison_Clean_111825.csv",
          row.names = FALSE)


# ============================================================
# Create clean model table with BH and BY adjusted p-values
# ============================================================

# Add BH and BY to the full model output
model_comparison_df$BH <- p.adjust(model_comparison_df$LRT_pvalue, method = "BH")
model_comparison_df$BY <- p.adjust(model_comparison_df$LRT_pvalue, method = "BY")

# Filter out failed or nonsensical models
clean_models <- subset(
  model_comparison_df,
  !is.na(AIC_full) & !is.na(AIC_reduced) &
    abs(AIC_full) < 1e5 & abs(AIC_reduced) < 1e5
)

# Add ΔAIC
clean_models$delta_AIC <- clean_models$AIC_full - clean_models$AIC_reduced

# Save clean model table
write.csv(
  clean_models,
  "NB_Model_25NOV25_Clean_WithBHBY.csv",
  row.names = FALSE
)

cat("\n✅ Clean model file created with BH and BY corrections.\n")
cat("   Total clean models:", nrow(clean_models), "\n")



# Count rows in clean_models object in memory
nrow(clean_models)

# Count rows in CSV after re-reading it
clean_csv <- read.csv("NB_Model_Comparison_Clean_WithBHBY.csv")
nrow(clean_csv)

# Which feature is missing?
setdiff(clean_models$Feature, clean_csv$Feature)