# =============================================================================
# UNIFIED CONVERGENCE ANALYSIS: Microbiome & Resistome
# Comparing distance to beef reference over dietary intervention
# =============================================================================

library(vegan)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)

# =============================================================================
# FUNCTION DEFINITIONS
# =============================================================================

#' Calculate mean distance to reference samples (direct method)
#' Works with any pairwise dissimilarity metric
mean_dist_to_reference <- function(dist_mat, sample_ids, reference_ids) {
  distances <- sapply(sample_ids, function(s) {
    if(s %in% reference_ids) {
      other_ref <- setdiff(reference_ids, s)
      if(length(other_ref) > 0) mean(dist_mat[s, other_ref]) else NA
    } else {
      mean(dist_mat[s, reference_ids])
    }
  })
  names(distances) <- sample_ids
  return(distances)
}


#' Calculate distance to group centroid using PCoA embedding
#' Returns distances and PCoA coordinates for plotting
centroid_dist_pcoa <- function(dist_obj, metadata, group_var, reference_group, n_axes = NULL) {
  
  # Perform PCoA
  pcoa <- cmdscale(dist_obj, k = nrow(as.matrix(dist_obj)) - 1, eig = TRUE)
  
  # Determine axes to use
  if(is.null(n_axes)) {
    n_axes <- sum(pcoa$eig > 0)
  }
  n_axes <- min(n_axes, ncol(pcoa$points))
  
  coords <- pcoa$points[, 1:n_axes]
  
  # Variance explained
  var_explained <- pcoa$eig / sum(abs(pcoa$eig)) * 100
  total_var <- sum(var_explained[1:n_axes])
  
  # Calculate centroid - use sample_id if available, otherwise rownames
  if("sample_id" %in% colnames(metadata)) {
    reference_samples <- metadata$sample_id[metadata[[group_var]] == reference_group]
  } else {
    reference_samples <- rownames(metadata)[metadata[[group_var]] == reference_group]
  }
  centroid <- colMeans(coords[reference_samples, , drop = FALSE])
  
  # Euclidean distance to centroid
  distances <- apply(coords, 1, function(x) sqrt(sum((x - centroid)^2)))
  
  return(list(
    distances = distances,
    coords = coords,
    centroid = centroid,
    var_explained = var_explained,
    total_var = total_var,
    n_axes = n_axes
  ))
}


#' Calculate multiple distance metrics
calculate_distances <- function(ps, include_unifrac = TRUE) {
  
  otu <- as(otu_table(ps), "matrix")
  if(taxa_are_rows(ps)) otu <- t(otu)
  
  distances <- list(
    bray = vegdist(otu, method = "bray"),
    jaccard = vegdist(otu, method = "jaccard", binary = TRUE)
  )
  
  if(include_unifrac) {
    if(!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
      cat("Calculating UniFrac distances...\n")
      distances$unifrac <- UniFrac(ps, weighted = FALSE)
      distances$wunifrac <- UniFrac(ps, weighted = TRUE)
    } else {
      warning("No phylogenetic tree found. Skipping UniFrac.")
    }
  }
  
  return(distances)
}


#' Main analysis function for convergence toward beef
analyze_beef_convergence <- function(ps, treatment_name,
                                     human_timepoints = c("Start", "Midpoint", "End"),
                                     beef_timepoint,
                                     timepoint_var = "timepoint",
                                     include_unifrac = TRUE,
                                     data_type = "Microbiome") {
  
  cat("\n=============================================================================\n")
  cat("ANALYZING:", data_type, "-", treatment_name, "\n")
  cat("Beef reference:", beef_timepoint, "\n")
  cat("=============================================================================\n")
  
  # Get metadata and preserve rownames
  meta <- as(sample_data(ps), "data.frame")
  meta$sample_id <- rownames(meta)
  
  # Identify samples
  beef_samples <- meta$sample_id[meta[[timepoint_var]] == beef_timepoint]
  human_samples <- meta$sample_id[meta[[timepoint_var]] %in% human_timepoints]
  all_samples <- c(human_samples, beef_samples)
  
  cat("\nSamples: Beef =", length(beef_samples), ", Human =", length(human_samples), "\n")
  
  # Calculate distances
  dist_list <- calculate_distances(ps, include_unifrac = include_unifrac)
  
  # ---------------------------------------------------------------------------
  # Calculate distances to beef (both methods)
  # ---------------------------------------------------------------------------
  
  for(metric_name in names(dist_list)) {
    dist_mat <- as.matrix(dist_list[[metric_name]])
    
    # Method 1: Mean distance
    col_mean <- paste0(metric_name, "_mean_to_beef")
    distances_mean <- mean_dist_to_reference(dist_mat, all_samples, beef_samples)
    meta[[col_mean]] <- distances_mean[meta$sample_id]
    
    # Method 2: Centroid distance
    col_centroid <- paste0(metric_name, "_centroid_to_beef")
    pcoa_result <- centroid_dist_pcoa(dist_list[[metric_name]], meta, timepoint_var, beef_timepoint)
    meta[[col_centroid]] <- pcoa_result$distances[meta$sample_id]
  }
  
  # Store PCoA for Bray-Curtis (for ordination plot)
  pcoa_bray <- centroid_dist_pcoa(dist_list$bray, meta, timepoint_var, beef_timepoint)
  
  # ---------------------------------------------------------------------------
  # Subset to human samples (preserve rownames through tibble conversion)
  # ---------------------------------------------------------------------------
  meta_human <- meta[meta[[timepoint_var]] %in% human_timepoints, ]
  meta_human[[timepoint_var]] <- factor(meta_human[[timepoint_var]], levels = human_timepoints)
  rownames(meta_human) <- meta_human$sample_id
  
  # Diagnostic: Check that distance columns have values
  cat("\nDiagnostic - Human samples with non-NA distances:\n")
  for(metric_name in names(dist_list)) {
    col_mean <- paste0(metric_name, "_mean_to_beef")
    n_valid <- sum(!is.na(meta_human[[col_mean]]))
    cat("  ", metric_name, ":", n_valid, "/", nrow(meta_human), "\n")
  }
  
  dist_list_human <- lapply(dist_list, function(d) {
    as.dist(as.matrix(d)[human_samples, human_samples])
  })
  
  # ---------------------------------------------------------------------------
  # Statistical tests (run for BOTH methods)
  # ---------------------------------------------------------------------------
  results <- list()
  
  for(metric_name in names(dist_list)) {
    results[[metric_name]] <- list()
    
    for(method in c("mean", "centroid")) {
      col_name <- paste0(metric_name, "_", method, "_to_beef")
      
      # Kruskal-Wallis
      kw_test <- kruskal.test(
        as.formula(paste(col_name, "~", timepoint_var)), 
        data = meta_human
      )
      
      # Linear trend
      meta_human$time_numeric <- as.numeric(meta_human[[timepoint_var]])
      lm_test <- lm(as.formula(paste(col_name, "~ time_numeric")), data = meta_human)
      lm_summary <- summary(lm_test)
      
      # Store results for this method
      results[[metric_name]][[method]] <- list(
        kruskal = kw_test,
        linear_model = lm_test,
        slope = unname(coef(lm_test)[2]),
        slope_p = coef(lm_summary)[2, 4]
      )
    }
    
    # Mantel test (same for both methods - tests overall temporal structure)
    timepoint_rank <- as.numeric(meta_human[[timepoint_var]])
    timepoint_dist <- dist(timepoint_rank)
    mantel_test <- mantel(dist_list_human[[metric_name]], timepoint_dist, 
                          method = "spearman", permutations = 999)
    results[[metric_name]]$mantel <- mantel_test
  }
  
  # Calculate beef within-group variation
  beef_within <- list()
  for(metric_name in names(dist_list)) {
    col_name <- paste0(metric_name, "_mean_to_beef")
    beef_within[[metric_name]] <- mean(meta[[col_name]][meta[[timepoint_var]] == beef_timepoint], na.rm = TRUE)
  }
  
  return(list(
    treatment = treatment_name,
    data_type = data_type,
    beef_timepoint = beef_timepoint,
    timepoint_var = timepoint_var,
    meta = meta,
    meta_human = meta_human,
    dist_list = dist_list,
    beef_samples = beef_samples,
    human_samples = human_samples,
    human_timepoints = human_timepoints,
    results = results,
    beef_within = beef_within,
    pcoa_bray = pcoa_bray
  ))
}


# =============================================================================
# PLOTTING FUNCTIONS (Panel-ready)
# =============================================================================

#' Boxplot WITHOUT beef samples (no reference line)
plot_convergence_boxplot <- function(analysis_result, metric = "bray", 
                                     method = "mean", show_title = TRUE) {
  
  col_name <- paste0(metric, "_", method, "_to_beef")
  meta_human <- analysis_result$meta_human
  timepoint_var <- analysis_result$timepoint_var
  human_tps <- analysis_result$human_timepoints
  
  # Color palette for human timepoints
  n_human <- length(human_tps)
  human_colors <- c("#559e83", "#1b85b8", "#5a5255")[1:n_human]
  names(human_colors) <- human_tps
  
  # Metric labels
  metric_labels <- c(
    bray = "Bray-Curtis",
    jaccard = "Jaccard", 
    unifrac = "Unweighted UniFrac",
    wunifrac = "Weighted UniFrac"
  )
  
  # Method label for y-axis
  method_label <- ifelse(method == "mean", "Mean", "Centroid")
  
  p <- ggplot(meta_human, aes(x = .data[[timepoint_var]], y = .data[[col_name]], 
                              fill = .data[[timepoint_var]])) +
    # Boxplots
    #geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    stat_summary(
      fun.data = function(x) {
        data.frame(
          ymin = quantile(x, 0.05),
          lower = quantile(x, 0.25),
          middle = quantile(x, 0.50),
          upper = quantile(x, 0.75),
          ymax = quantile(x, 0.95)
        )
      },
      geom = "boxplot",
      alpha = 0.7
    ) +
    # Individual points
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    # Mean trajectory line
    stat_summary(fun = mean, geom = "line", aes(group = 1), 
                 linewidth = 1.2, color = "black") +
    stat_summary(fun = mean, geom = "point", size = 4, color = "black") +
    # Styling
    scale_fill_manual(values = human_colors) +
    labs(
      x = "Timepoint",
      y = paste0(method_label, " ", metric_labels[metric], " Distance to Beef")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 10, 5, 10)
    )
  
  if(show_title) {
    p <- p + labs(
      title = paste0(analysis_result$treatment, ": ", analysis_result$data_type)
    ) +
      theme(plot.title = element_text(face = "bold", size = 11))
  }
  
  return(p)
}


#' Treatment comparison line plot (for combining microbiome + resistome)
plot_treatment_comparison <- function(results_list, metric = "bray", method = "mean",
                                      y_limits = NULL, show_title = TRUE) {
  
  col_name <- paste0(metric, "_", method, "_to_beef")
  
  # Combine data from all results using base R to avoid .data pronoun issues
  meta_list <- lapply(results_list, function(r) {
    df <- r$meta_human
    data.frame(
      treatment = r$treatment,
      data_type = r$data_type,
      timepoint = df[[r$timepoint_var]],
      distance = df[[col_name]],
      stringsAsFactors = FALSE
    )
  })
  meta_combined <- do.call(rbind, meta_list)
  
  # Get unique data type for title
  data_type <- unique(meta_combined$data_type)[1]
  
  # Color palette
  treatments <- unique(meta_combined$treatment)
  treatment_palette <- c("#ae5a41", "#7b68ee")[1:length(treatments)]
  names(treatment_palette) <- treatments
  
  # Y-axis label based on method
  method_label <- ifelse(method == "mean", "Mean", "Centroid")
  metric_label <- c(bray = "Bray-Curtis", jaccard = "Jaccard", 
                    unifrac = "Unweighted UniFrac", wunifrac = "Weighted UniFrac")[metric]
  
  p <- ggplot(meta_combined, aes(x = timepoint, y = distance,
                                 color = treatment, group = treatment)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
    stat_summary(fun = mean, geom = "point", size = 4) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.12) +
    scale_color_manual(values = treatment_palette) +
    labs(
      x = "Timepoint",
      y = paste0(method_label, " ", metric_label, " Distance to Beef"),
      color = "Treatment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(5, 10, 5, 10)
    )
  
  if(!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  if(show_title) {
    p <- p + labs(
      title = paste0(data_type, ": Convergence by Treatment"),
      subtitle = "Lower values = more similar to beef"
    ) +
      theme(plot.title = element_text(face = "bold", size = 11))
  }
  
  return(p)
}


#' PCoA ordination with centroid
plot_ordination_centroid <- function(analysis_result, show_title = TRUE) {
  
  pcoa <- analysis_result$pcoa_bray
  meta <- analysis_result$meta
  timepoint_var <- analysis_result$timepoint_var
  beef_tp <- analysis_result$beef_timepoint
  human_tps <- analysis_result$human_timepoints
  
  # Prepare plot data - use sample_id column for reliable indexing
  all_tps <- c(human_tps, beef_tp)
  plot_samples <- meta$sample_id[meta[[timepoint_var]] %in% all_tps]
  
  pcoa_df <- data.frame(
    PC1 = pcoa$coords[plot_samples, 1],
    PC2 = pcoa$coords[plot_samples, 2],
    timepoint = meta[[timepoint_var]][match(plot_samples, meta$sample_id)],
    sample_type = ifelse(meta[[timepoint_var]][match(plot_samples, meta$sample_id)] == beef_tp, "Beef", "Human")
  )
  pcoa_df$timepoint <- factor(pcoa_df$timepoint, levels = all_tps)
  
  # Centroid
  centroid_df <- data.frame(
    PC1 = pcoa$centroid[1],
    PC2 = pcoa$centroid[2]
  )
  
  # Colors
  n_human <- length(human_tps)
  human_colors <- c("#559e83", "#1b85b8", "#5a5255")[1:n_human]
  names(human_colors) <- human_tps
  all_colors <- c(human_colors, setNames("#8B4513", beef_tp))
  
  # Variance explained
  var1 <- round(pcoa$var_explained[1], 1)
  var2 <- round(pcoa$var_explained[2], 1)
  
  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
    # Lines to centroid (human samples only)
    geom_segment(data = pcoa_df %>% filter(sample_type == "Human"),
                 aes(xend = centroid_df$PC1, yend = centroid_df$PC2, color = timepoint),
                 alpha = 0.3, linewidth = 0.5) +
    # Points
    geom_point(aes(color = timepoint, shape = sample_type), size = 3, alpha = 0.8) +
    # Centroid
    geom_point(data = centroid_df, aes(x = PC1, y = PC2),
               color = "#8B4513", size = 6, shape = 4, stroke = 2) +
    # Styling
    scale_color_manual(values = all_colors) +
    scale_shape_manual(values = c("Human" = 16, "Beef" = 17)) +
    labs(
      x = paste0("PC1 (", var1, "%)"),
      y = paste0("PC2 (", var2, "%)"),
      color = "Timepoint",
      shape = "Sample Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      plot.margin = margin(5, 10, 5, 10)
    )
  
  if(show_title) {
    p <- p + labs(
      title = paste0(analysis_result$treatment, ": ", analysis_result$data_type, " PCoA"),
      subtitle = "Lines show distance to beef centroid"
    ) +
      theme(plot.title = element_text(face = "bold", size = 11))
  }
  
  return(p)
}


#' Build summary table (includes both mean and centroid methods)
build_summary_table <- function(results_list) {
  
  summary_df <- data.frame(
    Data_Type = character(),
    Metric = character(),
    Method = character(),
    Treatment = character(),
    Slope = numeric(),
    Slope_p = numeric(),
    Direction = character(),
    KW_p = numeric(),
    Mantel_r = numeric(),
    Mantel_p = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(result in results_list) {
    for(metric in names(result$results)) {
      # Get Mantel test (same for both methods)
      mantel_res <- result$results[[metric]]$mantel
      
      for(method in c("mean", "centroid")) {
        res <- result$results[[metric]][[method]]
        
        direction <- ifelse(res$slope < 0 & res$slope_p < 0.05, "CONVERGENCE",
                            ifelse(res$slope > 0 & res$slope_p < 0.05, "DIVERGENCE", "STABLE"))
        
        summary_df <- rbind(summary_df, data.frame(
          Data_Type = result$data_type,
          Metric = metric,
          Method = method,
          Treatment = result$treatment,
          Slope = round(res$slope, 4),
          Slope_p = round(res$slope_p, 4),
          Direction = direction,
          KW_p = round(res$kruskal$p.value, 4),
          Mantel_r = round(mantel_res$statistic, 4),
          Mantel_p = round(mantel_res$signif, 4),
          row.names = NULL
        ))
      }
    }
  }
  
  return(summary_df)
}


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

# Define common parameters
human_timepoints <- c("Start", "Midpoint", "End")

# =============================================================================
# CSS NORMALIZATION ON FULL DATASETS FIRST
# =============================================================================

cat("\n\n########## CSS NORMALIZATION ##########\n")

# Microbiome: CSS normalize the full dataset, then subset
cat("Normalizing full microbiome dataset...\n")
data_micro.ps.css <- phyloseq_transform_css(data_micro.ps, log = FALSE)

# Resistome: Rename "Mid" to "Midpoint" for consistency BEFORE normalization
cat("Normalizing full resistome dataset...\n")
sample_data(AMR_data.ps)$col_time <- gsub("^Mid$", "Midpoint", sample_data(AMR_data.ps)$col_time)
AMR_data.ps.css <- phyloseq_transform_css(AMR_data.ps, log = FALSE)

AMR_group.ps.css <- tax_glom(AMR_data.ps.css, "group")

# =============================================================================
# MICROBIOME ANALYSIS
# =============================================================================

cat("\n\n########## MICROBIOME ANALYSIS ##########\n")

# No Claim - subset from CSS-normalized data
noclaim_micro.ps <- subset_samples(data_micro.ps.css, trt == "No Claim" & 
                                     timepoint %in% c(human_timepoints, "No_Claim"))
noclaim_micro.ps <- prune_taxa(taxa_sums(noclaim_micro.ps) > 0, noclaim_micro.ps)

results_micro_noclaim <- analyze_beef_convergence(
  ps = noclaim_micro.ps,
  treatment_name = "No Claim",
  human_timepoints = human_timepoints,
  beef_timepoint = "No_Claim",
  timepoint_var = "timepoint",
  include_unifrac = TRUE,
  data_type = "Microbiome"
)

# RWA - subset from CSS-normalized data
rwa_micro.ps <- subset_samples(data_micro.ps.css, trt == "RWA" & 
                                 timepoint %in% c(human_timepoints, "RWA"))
rwa_micro.ps <- prune_taxa(taxa_sums(rwa_micro.ps) > 0, rwa_micro.ps)

results_micro_rwa <- analyze_beef_convergence(
  ps = rwa_micro.ps,
  treatment_name = "RWA",
  human_timepoints = human_timepoints,
  beef_timepoint = "RWA",
  timepoint_var = "timepoint",
  include_unifrac = TRUE,
  data_type = "Microbiome"
)

# =============================================================================
# RESISTOME ANALYSIS
# =============================================================================

cat("\n\n########## RESISTOME ANALYSIS ##########\n")

# No Claim Beef - subset from CSS-normalized data
noclaim_amr.ps <- subset_samples(AMR_group.ps.css, clean_treatment == "No Claim" & 
                                   col_time %in% c(human_timepoints, "No Claim Beef"))
noclaim_amr.ps <- prune_taxa(taxa_sums(noclaim_amr.ps) > 0, noclaim_amr.ps)

results_amr_noclaim <- analyze_beef_convergence(
  ps = noclaim_amr.ps,
  treatment_name = "No Claim",
  human_timepoints = human_timepoints,
  beef_timepoint = "No Claim Beef",
  timepoint_var = "col_time",
  include_unifrac = FALSE,  # No tree for resistome
  data_type = "Resistome"
)

# RWA Beef - subset from CSS-normalized data
rwa_amr.ps <- subset_samples(AMR_group.ps.css, clean_treatment == "RWA" & 
                               col_time %in% c(human_timepoints, "RWA Beef"))
rwa_amr.ps <- prune_taxa(taxa_sums(rwa_amr.ps) > 0, rwa_amr.ps)

results_amr_rwa <- analyze_beef_convergence(
  ps = rwa_amr.ps,
  treatment_name = "RWA",
  human_timepoints = human_timepoints,
  beef_timepoint = "RWA Beef",
  timepoint_var = "col_time",
  include_unifrac = FALSE,
  data_type = "Resistome"
)

# =============================================================================
# SUMMARY TABLES
# =============================================================================

microbiome_summary <- build_summary_table(list(results_micro_noclaim, results_micro_rwa))
resistome_summary <- build_summary_table(list(results_amr_noclaim, results_amr_rwa))

cat("\n\n=== MICROBIOME SUMMARY ===\n")
print(microbiome_summary)

cat("\n\n=== RESISTOME SUMMARY ===\n")
print(resistome_summary)

# Combined table
full_summary <- rbind(microbiome_summary, resistome_summary)
cat("\n\n=== COMBINED SUMMARY ===\n")
print(full_summary)

# Create filtered tables for each method (for easier reading)
summary_mean <- full_summary[full_summary$Method == "mean", ]
summary_centroid <- full_summary[full_summary$Method == "centroid", ]

cat("\n\n=== MEAN DISTANCE METHOD ONLY ===\n")
print(summary_mean)

cat("\n\n=== CENTROID METHOD ONLY ===\n")
print(summary_centroid)


# =============================================================================
# PANEL FIGURES
# =============================================================================

# ---------------------------------------------------------------------------
# Figure 1: Treatment comparison - Microbiome (left) vs Resistome (right)
# WITH MATCHED Y-AXES
# ---------------------------------------------------------------------------

# Get y-axis ranges for BOTH datasets to set common limits
micro_values <- c(results_micro_noclaim$meta_human$bray_mean_to_beef,
                  results_micro_rwa$meta_human$bray_mean_to_beef)
amr_values <- c(results_amr_noclaim$meta_human$bray_mean_to_beef,
                results_amr_rwa$meta_human$bray_mean_to_beef)

# Find global range across both datasets
all_values <- c(micro_values, amr_values)
global_min <- min(all_values, na.rm = TRUE)
global_max <- max(all_values, na.rm = TRUE)
# Add some padding
y_padding <- (global_max - global_min) * 0.05
global_limits <- c(global_min - y_padding, global_max + y_padding)

p_micro_compare <- plot_treatment_comparison(
  list(results_micro_noclaim, results_micro_rwa),
  y_limits = global_limits,
  show_title = TRUE
) + 
  labs(title = "A) Microbiome",
       subtitle = "Lower values = more similar to beef") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(plot.margin = margin(10, 15, 10, 10),
        axis.title.y = element_text(size = 11))

p_amr_compare <- plot_treatment_comparison(
  list(results_amr_noclaim, results_amr_rwa),
  y_limits = global_limits,
  show_title = TRUE
) + 
  labs(title = "B) Resistome",
       subtitle = "Lower values = more similar to beef") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(plot.margin = margin(10, 10, 10, 15),
        axis.title.y = element_text(size = 11))

# Combine with shared legend
p_combined_comparison <- p_micro_compare + p_amr_compare +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_combined_comparison)
ggsave("Figures/Succession/Combined_convergence_comparison.png", 
       p_combined_comparison, width = 12, height = 5, dpi = 150)


# ---------------------------------------------------------------------------
# Figure 1b: Same comparison but using CENTROID method
# ---------------------------------------------------------------------------

p_micro_compare_centroid <- plot_treatment_comparison(
  list(results_micro_noclaim, results_micro_rwa),
  method = "centroid",
  y_limits = NULL,  # Let it auto-scale for centroid
  show_title = TRUE
) + 
  labs(title = "A) Microbiome (Centroid)",
       subtitle = "Lower values = more similar to beef") +
  theme(plot.margin = margin(10, 15, 10, 10),
        axis.title.y = element_text(size = 11))

p_amr_compare_centroid <- plot_treatment_comparison(
  list(results_amr_noclaim, results_amr_rwa),
  method = "centroid",
  y_limits = NULL,
  show_title = TRUE
) + 
  labs(title = "B) Resistome (Centroid)",
       subtitle = "Lower values = more similar to beef") +
  theme(plot.margin = margin(10, 10, 10, 15),
        axis.title.y = element_text(size = 11))

p_combined_comparison_centroid <- p_micro_compare_centroid + p_amr_compare_centroid +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_combined_comparison_centroid)
ggsave("Figures/Succession/Combined_convergence_comparison_centroid.png", 
       p_combined_comparison_centroid, width = 12, height = 5, dpi = 150)


# ---------------------------------------------------------------------------
# Figure 2: Boxplots - No beef samples, organized by data type and treatment
# ---------------------------------------------------------------------------

p_micro_nc_box <- plot_convergence_boxplot(results_micro_noclaim) + 
  labs(title = "A) Microbiome - No Claim")
p_micro_rwa_box <- plot_convergence_boxplot(results_micro_rwa) + 
  labs(title = "B) Microbiome - RWA")
p_amr_nc_box <- plot_convergence_boxplot(results_amr_noclaim) + 
  labs(title = "C) Resistome - No Claim")
p_amr_rwa_box <- plot_convergence_boxplot(results_amr_rwa) + 
  labs(title = "D) Resistome - RWA")

p_all_boxplots <- (p_micro_nc_box + p_micro_rwa_box) / (p_amr_nc_box + p_amr_rwa_box)

print(p_all_boxplots)
ggsave("Figures/Succession/Combined_convergence_boxplots.png", 
       p_all_boxplots, width = 12, height = 10, dpi = 150)


# ---------------------------------------------------------------------------
# Figure 2b: Boxplots using CENTROID method
# ---------------------------------------------------------------------------

p_micro_nc_box_cent <- plot_convergence_boxplot(results_micro_noclaim, method = "centroid") + 
  labs(title = "A) Microbiome - No Claim")
p_micro_rwa_box_cent <- plot_convergence_boxplot(results_micro_rwa, method = "centroid") + 
  labs(title = "B) Microbiome - RWA")
p_amr_nc_box_cent <- plot_convergence_boxplot(results_amr_noclaim, method = "centroid") + 
  labs(title = "C) Resistome - No Claim")
p_amr_rwa_box_cent <- plot_convergence_boxplot(results_amr_rwa, method = "centroid") + 
  labs(title = "D) Resistome - RWA")

p_all_boxplots_centroid <- (p_micro_nc_box_cent + p_micro_rwa_box_cent) / (p_amr_nc_box_cent + p_amr_rwa_box_cent)

print(p_all_boxplots_centroid)
ggsave("Figures/Succession/Combined_convergence_boxplots_centroid.png", 
       p_all_boxplots_centroid, width = 12, height = 10, dpi = 150)


# ---------------------------------------------------------------------------
# Figure 3: PCoA Ordinations with centroids
# ---------------------------------------------------------------------------

p_micro_nc_pcoa <- plot_ordination_centroid(results_micro_noclaim) + 
  labs(title = "A) Microbiome - No Claim") +
  theme(legend.position = "none")
p_micro_rwa_pcoa <- plot_ordination_centroid(results_micro_rwa) + 
  labs(title = "B) Microbiome - RWA") +
  theme(legend.position = "none")
p_amr_nc_pcoa <- plot_ordination_centroid(results_amr_noclaim) + 
  labs(title = "C) Resistome - No Claim") +
  theme(legend.position = "none")
p_amr_rwa_pcoa <- plot_ordination_centroid(results_amr_rwa) + 
  labs(title = "D) Resistome - RWA")

# Extract legend from last plot
legend_pcoa <- get_legend(p_amr_rwa_pcoa)
p_amr_rwa_pcoa <- p_amr_rwa_pcoa + theme(legend.position = "none")

p_all_pcoa <- (p_micro_nc_pcoa + p_micro_rwa_pcoa) / (p_amr_nc_pcoa + p_amr_rwa_pcoa)

# Add shared legend at bottom
p_all_pcoa_legend <- plot_grid(
  p_all_pcoa,
  legend_pcoa,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

print(p_all_pcoa_legend)
ggsave("Figures/Succession/Combined_convergence_PCoA.png", 
       p_all_pcoa_legend, width = 12, height = 11, dpi = 150)


# ---------------------------------------------------------------------------
# Figure 4: Method comparison (Mean distance vs Centroid)
# ---------------------------------------------------------------------------

plot_method_comparison <- function(result, metric = "bray") {
  
  col_mean <- paste0(metric, "_mean_to_beef")
  col_centroid <- paste0(metric, "_centroid_to_beef")
  timepoint_var <- result$timepoint_var
  
  meta_long <- result$meta_human %>%
    select(!!sym(timepoint_var), !!sym(col_mean), !!sym(col_centroid)) %>%
    pivot_longer(cols = c(!!sym(col_mean), !!sym(col_centroid)),
                 names_to = "Method",
                 values_to = "Distance") %>%
    mutate(Method = recode(Method,
                           !!col_mean := "Mean Distance",
                           !!col_centroid := "PCoA Centroid"))
  
  p <- ggplot(meta_long, aes(x = .data[[timepoint_var]], y = Distance, 
                             color = Method, group = Method)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15) +
    scale_color_manual(values = c("Mean Distance" = "#d95f02", "PCoA Centroid" = "#7570b3")) +
    labs(
      x = "Timepoint",
      y = "Distance to Beef",
      title = paste0(result$treatment, ": ", result$data_type),
      color = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "bottom"
    )
  
  return(p)
}

p_method_micro_nc <- plot_method_comparison(results_micro_noclaim) + theme(legend.position = "none")
p_method_micro_rwa <- plot_method_comparison(results_micro_rwa) + theme(legend.position = "none")
p_method_amr_nc <- plot_method_comparison(results_amr_noclaim) + theme(legend.position = "none")
p_method_amr_rwa <- plot_method_comparison(results_amr_rwa)

legend_method <- get_legend(p_method_amr_rwa)
p_method_amr_rwa <- p_method_amr_rwa + theme(legend.position = "none")

p_method_all <- (p_method_micro_nc + p_method_micro_rwa) / (p_method_amr_nc + p_method_amr_rwa)

p_method_all_legend <- plot_grid(
  p_method_all,
  legend_method,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

print(p_method_all_legend)
ggsave("Figures/Succession/Combined_method_comparison.png", 
       p_method_all_legend, width = 10, height = 9, dpi = 150)


# ---------------------------------------------------------------------------
# Figure 5: All metrics comparison (Microbiome only, has UniFrac)
# ---------------------------------------------------------------------------

plot_all_metrics_panel <- function(result) {
  
  metrics <- names(result$results)
  metric_labels <- c(bray = "Bray-Curtis", jaccard = "Jaccard", 
                     unifrac = "Unweighted UniFrac", wunifrac = "Weighted UniFrac")
  
  plots <- lapply(metrics, function(m) {
    plot_convergence_boxplot(result, metric = m, show_title = FALSE) +
      labs(title = metric_labels[m]) +
      theme(plot.title = element_text(size = 10, face = "bold"))
  })
  
  combined <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = paste0(result$treatment, ": ", result$data_type, " - All Metrics"),
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )
  
  return(combined)
}

p_micro_nc_allmetrics <- plot_all_metrics_panel(results_micro_noclaim)
print(p_micro_nc_allmetrics)
ggsave("Figures/Succession/Microbiome_NoClaimallmetrics.png", 
       p_micro_nc_allmetrics, width = 10, height = 8, dpi = 150)

p_micro_rwa_allmetrics <- plot_all_metrics_panel(results_micro_rwa)
print(p_micro_rwa_allmetrics)
ggsave("Figures/Succession/Microbiome_RWA_allmetrics.png", 
       p_micro_rwa_allmetrics, width = 10, height = 8, dpi = 150)


# =============================================================================
# SAVE RESULTS
# =============================================================================

write.csv(full_summary, "Figures/Succession/Convergence_analysis_summary_all.csv", row.names = FALSE)
write.csv(summary_mean, "Figures/Succession/Convergence_analysis_summary_mean_distance.csv", row.names = FALSE)
write.csv(summary_centroid, "Figures/Succession/Convergence_analysis_summary_centroid.csv", row.names = FALSE)

cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("Figures saved to: Figures/Succession/\n")
cat("Tables saved to: Tables/\n")
cat("  - Convergence_analysis_summary_all.csv (both methods)\n")
cat("  - Convergence_analysis_summary_mean_distance.csv\n")
cat("  - Convergence_analysis_summary_centroid.csv\n")