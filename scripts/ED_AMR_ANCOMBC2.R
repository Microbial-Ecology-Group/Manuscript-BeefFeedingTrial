# HFT ANCOMBC
set.seed(123)
library(ANCOMBC)
library(microbiome)
library(dplyr)
library(stringr)
library(cowplot)

###
##### Start ANCOMBC - Fecal samples, comparing by Combined_order #####
###

# only keep fecal samples and remove Final_washout
fecal_data.ps <- subset_samples(AMR_data.ps, Sample_Type == "Fecal")
clean_fecal_data.ps <- subset_samples(fecal_data.ps, col_time != "Washout")
clean_fecal_data.ps <- prune_taxa(taxa_sums(clean_fecal_data.ps) > 0, clean_fecal_data.ps)
any(taxa_sums(clean_fecal_data.ps)==0) # QUADRUPLE CHECKING - nope good.


sample_data(clean_fecal_data.ps)$Combined_order <- factor(sample_data(clean_fecal_data.ps)$Combined_order, levels = c("Start_No_Claim","Midpoint_No_Claim",
                                                                                                          "End_No_Claim","Start_RWA","Midpoint_RWA",
                                                                                                          "End_RWA"))
sample_data(clean_fecal_data.ps)$Participant_ID <- as.factor(sample_data(clean_fecal_data.ps)$Participant_ID)
sample_data(clean_fecal_data.ps)$Timepoint


# Trim taxa
#filtered_fecal_data_noSNP_group <- preDA(fecal_data_noSNP_group, min.samples=2, min.reads = 1) # 10% of samples (77/154)
# 203 taxa left

# Change the most specific 
colnames(phyloseq::tax_table(clean_fecal_data.ps))[4] <- "Species"


# This doesn't work by Timepoint, has to be "Start, Mid, End" variables
ANCOM_feces_output_full_interaction = ancombc2(data = clean_fecal_data.ps, assay_name = "counts", tax_level = "Species",
                              fix_formula = "col_time * trt_id * treatment_order", rand_formula = "(1 | Participant_ID)",
                              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, alpha = 0.05,n_cl= 8,
                              group = "trt_id", struc_zero = TRUE, neg_lb = TRUE,verbose = TRUE,
                              global = TRUE, pairwise = FALSE)


res_primary <- ANCOM_feces_output_full_interaction$res

# Check column names
colnames(res_primary)

# Find the three-way interaction columns
grep("col_time.*trt_id.*treatment_order", colnames(res_primary), value = TRUE)

# Extract q-values (adjusted p-values) for the interaction
# Column name will be something like: q_col_timeMid:trt_idRWA:treatment_orderOrder2
q_cols_interaction <- grep("^q_.*:.*:.*treatment_order", colnames(res_primary), value = TRUE)
print(q_cols_interaction)

# Summarize significance for each interaction term
for (col in q_cols_interaction) {
  q_vals <- res_primary[[col]]
  total <- sum(!is.na(q_vals))
  sig <- sum(q_vals < 0.05, na.rm = TRUE)
  cat(col, ":", sig, "/", total, "significant\n")
}

ANCOM_feces_output_double_interaction_trt_order = ancombc2(data = clean_fecal_data.ps, assay_name = "counts", tax_level = "Species",
                                               fix_formula = "col_time + trt_id * treatment_order", rand_formula = "(1 | Participant_ID)",
                                               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, alpha = 0.05,n_cl= 8,
                                               group = "trt_id", struc_zero = TRUE, neg_lb = TRUE,verbose = TRUE,
                                               global = TRUE, pairwise = FALSE)


# Model for pairwise comparisons by Combined_order
ANCOM_feces_output = ancombc2(data = clean_fecal_data.ps, assay_name = "counts", tax_level = "Species", 
                                  fix_formula = "Combined_order", rand_formula = "(1 | Participant_ID)",
                                  p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                                  group = "Combined_order", struc_zero = TRUE, neg_lb = TRUE,
                                  alpha = 0.05, n_cl = 10, verbose = TRUE,
                                  global = TRUE, pairwise = TRUE,
                                  dunnet = TRUE, trend = FALSE)

# Extract output for pairwise comparisons
res_pair = ANCOM_feces_output$res_pair

df_Combined_order = res_pair %>%
  dplyr::select(taxon, contains("Combined_order")) 


df_fig_Combined_order <- df_Combined_order %>%
  # Step 1: Apply the condition for the lfc_Combined_order columns
  dplyr::mutate(dplyr::across(starts_with("lfc_Combined_order"), ~ifelse(get(gsub("lfc_", "diff_", cur_column())) == 1, .x,.x)),
                # Step 2: Round the lfc_Combined_order columns and create _rounded columns
                across(starts_with("lfc_Combined_order"), ~round(.x, 3), .names = "{.col}_rounded"),
                # Step 3: Assign colors based on diff_Combined_order and passed_ss columns
                across(
                  starts_with("diff_Combined_order"),
                  ~ dplyr::case_when(
                    .x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 1 ~ "darkred",
                    .x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 0 ~ "#blue",
                    .x == 0 ~ "#333333"
                  ),
                  .names = "{.col}_color")) %>%
  # Step 4: Pivot the _rounded columns into long format
  pivot_longer(cols = contains("_rounded"), names_to = "group", values_to = "value", names_prefix = "lfc_") %>%
  # Step 5: Pivot the _color columns into long format
  pivot_longer(cols = contains("_color"), names_to = "color_group", values_to = "color", names_prefix = "diff_") %>%
  # Step 6: Filter to ensure the group and color_group match
  filter(gsub("_rounded", "", group) == gsub("_color", "", color_group)) %>%
  # Step 7: Select the relevant columns and arrange by taxon
  dplyr::select(taxon, group, value, color) %>%
  arrange(taxon)


# FOR MOLLY ###
# These are the groups you can pick from: unique(df_fig_Combined_order$group)
# add the ones for RWA
unique(df_fig_Combined_order$group)

df_fig_Combined_order_filtered <- df_fig_Combined_order %>%
  dplyr::filter(str_ends(group, "_rounded")) %>%
  dplyr::mutate(group = dplyr::case_when(
    group == "Combined_orderMidpoint_No_Claim_rounded" ~ "No Claim: Start to Midpoint",
    group == "Combined_orderEnd_No_Claim_Combined_orderMidpoint_No_Claim_rounded" ~ "No Claim: Midpoint to End",
    group == "Combined_orderEnd_No_Claim_rounded" ~ "No Claim: Start to End",
    group == "Combined_orderMidpoint_RWA_Combined_orderStart_RWA_rounded" ~ "RWA: Start to Midpoint",
    group == "Combined_orderEnd_RWA_Combined_orderMidpoint_RWA_rounded" ~ "RWA: Midpoint to End",
    group == "Combined_orderEnd_RWA_Combined_orderStart_RWA_rounded" ~ "RWA: Start to End",    
    TRUE ~ group)) %>%
  dplyr::filter(group %in% c("No Claim: Start to Midpoint", 
                             "No Claim: Midpoint to End", 
                             "No Claim: Start to End",
                             "RWA: Start to Midpoint", 
                             "RWA: Midpoint to End", 
                             "RWA: Start to End" )) %>%
  droplevels()




levels(as.factor(df_fig_Combined_order_filtered$group))

#### feces Dotplot - but with bias corrected abundances #######
log_table_feces <- ANCOM_feces_output$bias_correct_log_table

# Replace NA values in the log-transformed data with 0
log_table_feces[is.na(log_table_feces)] <- 0

# Step 1: Exponentiate the log-transformed data to get values on the original scale
pseudo_counts_feces <- exp(log_table_feces)

# Step 2: Scale the pseudo-counts (e.g., normalize each sample to a total count)
# For example, you can scale the pseudo-counts to a total library size or another baseline.

# Calculate the total sum of counts for each sample before scaling
total_counts_feces <- apply(pseudo_counts_feces, 2, sum)

# Define a library size (or use the sum of pseudo_counts as the new total count)
desired_total_count_feces <- 100000  # For example, set each sample to have a total of 10,000 counts

# Scale the pseudo-counts to the desired total count per sample
scaled_counts_feces <- t(t(pseudo_counts_feces) * (desired_total_count_feces / total_counts_feces))

# No_Claimert scaled_counts_feces to a data frame
scaled_counts_feces.df <- as.data.frame(scaled_counts_feces)

# Add a column for taxa names
scaled_counts_feces.df$taxon <- rownames(scaled_counts_feces.df)

# No_Claimert from wide to long format using pivot_longer
scaled_counts_feces.long <- scaled_counts_feces.df %>%
  pivot_longer(cols = -taxon, names_to = "Sample", values_to = "BiasAdj_counts")

# Make melted data
clean_fecal_data.melt <- psmelt(clean_fecal_data.ps)

# Perform the merge with 'taxon' from scaled_counts_feces.long and 'Species' from clean_fecal_data.melt
BiasAdj_clean_fecal_data.melt <- clean_fecal_data.melt %>%
  dplyr::left_join(scaled_counts_feces.long, by = c("Species" = "taxon", "Sample"))

# Step 1: Assign "Midpoint to End" and other labels
Species_avg_counts_feces_by_timepoint <- BiasAdj_clean_fecal_data.melt %>%
  dplyr::group_by(Species, Combined_order) %>%
  dplyr::summarise(Avg_BiasAdj_Count = mean(BiasAdj_counts, na.rm = FALSE), .groups = "drop") %>%
  
  # Assign "Midpoint to End" and other labels
  dplyr::mutate(Combined_order_label = dplyr::case_when(
    Combined_order == "Midpoint_No_Claim" ~ "No Claim: Start to Midpoint",
    Combined_order == "End_No_Claim" ~ "No Claim: Midpoint to End",   # First label for End_No_Claim
    Combined_order == "Midpoint_RWA" ~ "RWA: Start to Midpoint",
    Combined_order == "End_RWA" ~ "RWA: Midpoint to End",     # First label for End_RWA
    TRUE ~ as.character(Combined_order)  # Keep other timepoint values unchanged
  ))

# Step 2: Duplicate rows for End_No_Claim and End_RWA, then modify the Combined_order_label in the duplicated rows
duplicated_rows <- Species_avg_counts_feces_by_timepoint %>%
  dplyr::filter(Combined_order %in% c("End_No_Claim", "End_RWA")) %>%
  dplyr::mutate(Combined_order_label = dplyr::case_when(
    Combined_order == "End_No_Claim" ~ "No Claim: Start to End",
    Combined_order == "End_RWA" ~ "RWA: Start to End",
    TRUE ~ Combined_order_label
  ))

# Step 3: Combine the original and duplicated rows
Species_avg_counts_feces_by_timepoint <- Species_avg_counts_feces_by_timepoint %>%
  dplyr::bind_rows(duplicated_rows)


# Check the unique values in the final column
unique(Species_avg_counts_feces_by_timepoint$Combined_order_label)


df_fig_Combined_order_filtered_merged_data <- df_fig_Combined_order_filtered %>%
  left_join(Species_avg_counts_feces_by_timepoint, by = c("group" = "Combined_order_label", "taxon" = "Species"))

factor_order_Combined_order = c("No Claim: Start to Midpoint", 
                                "No Claim: Midpoint to End", 
                                "No Claim: Start to End",
                                "RWA: Start to Midpoint", 
                                "RWA: Midpoint to End", 
                                "RWA: Start to End")

df_fig_Combined_order_filtered_merged_data$group = factor(df_fig_Combined_order_filtered_merged_data$group, levels = factor_order_Combined_order)


plot <- ggplot(df_fig_Combined_order_filtered_merged_data, aes(x = value, y = reorder(taxon, value), color = color)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Add vertical dashed line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, alpha = 0.7)) +  # Set the size based on Avg_BiasAdj_Count
  scale_color_identity() +  # Use color directly from the dataset
  scale_size_continuous(range = c(1, 10)) +  # Adjust the range of the dot sizes
  labs(x = "Log-Fold Change", y = "Taxon", title = "logFC of AMR groups in feces over time") +
  theme_minimal() +  # Set the theme to minimal
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.7),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),  # Remove strip background
    strip.text = element_text(size = 12, colour = "black"), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    legend.position = "none",  # Hide the legend as the color is used directly
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),  # Make panel background transparent
    panel.border = element_rect(colour = "black", fill = NA, linewidth = .2),  # Add border around facets
    panel.spacing = unit(1, "lines")  # Add space between facets
  ) +
  facet_wrap(~group, ncol = 3, scales = "fixed")  # Facet by group, with fixed y-axis scales
plot

#
#### Experimental #########
## The Part below is still experimental to add the taxa labels. It works, but everything is crowded,
# you'll have to do some modifying to the figure above too.
##
#
# Or we can run ANCOMBC with a higher prevalence threshold to reduce the total number of taxa. 
# Extract taxonomic data from the phyloseq object
tax_data_feces <- as.data.frame(tax_table(clean_fecal_data.ps))

# Join the taxonomic data with your data frame
df_with_taxonomy_feces <- df_fig_Combined_order_filtered_merged_data %>%
  left_join(tax_data_feces, by = c("taxon" = "Species"),relationship = "many-to-many")  # Assuming 'Species' is your key column to join

# Reorder the taxon based on class, mechanism, and species
df_with_taxonomy_feces <- df_with_taxonomy_feces %>%
  arrange(class, mechanism, taxon) %>%  # Sort the data by class, mechanism, and Species
  dplyr::mutate(class = factor(class, levels = rev(unique(class))))  # Reorder the 'taxon' factor based on the new order


df_with_taxonomy_feces$group = factor(df_with_taxonomy_feces$group, levels = factor_order_Combined_order)


## Trying to add taxonomy
main_species_plot_feces <- ggplot(df_with_taxonomy_feces, aes(x = value, y = class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Add vertical dashed line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  scale_color_identity() +  # Use color directly from the dataset
  scale_size_continuous(range = c(1, 6)) +  # Increase the size range to create more contrast
  labs(x = "Log-Fold Change", y = "Taxon") + #  title = "FECAL - logFC of AMR groups in CEF group over time"
  theme_minimal() +  # Set the theme to minimal
  theme(
    axis.text.y =element_text(size = 10, colour = "black"),  # Hide y-axis species labels (to show in taxonomy plot)
    axis.text.x = element_text(size = 14, angle = 45, vjust = 0.5),
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),  # Remove strip background
    strip.text = element_text(size = 12, colour = "black"), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Centered title
    legend.position = "none",  # Hide the legend as the color is used directly
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),  # Make panel background transparent
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),  # Add border around facets
    panel.spacing = unit(1, "lines")  # Add space between facets
  ) +
  facet_wrap(~group, ncol = 3, scales = "fixed")

# Here's basically the figure you need:
main_species_plot_feces

# Save the final plot
ggsave("Figures/MM_DA_feces_by_Combined_order.jpg", plot = main_species_plot_feces, width = 28, height = 22, dpi = 300)




#Mollys updated fig with different colors, bigger text, etc. 
main_species_plot_feces <- ggplot(df_with_taxonomy_feces, aes(x = value, y = class)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),
             position = position_dodge(width = 0.5)) +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = NULL) +  # ⬅️ Remove Y-axis title
  theme_minimal(base_size = 22) +  # ⬅️ Increase global text size
  theme(
    # Axis text & titles
    axis.text.y   = element_text(size = 22, colour = "black"),
    axis.text.x   = element_text(size = 24, angle = 45, vjust = 0.5, colour = "black"),
    axis.title.x  = element_text(size = 26, colour = "black"),
    axis.title.y  = element_blank(),
    
    # Facet strip text
    strip.background = element_blank(),
    strip.text = element_text(size = 26, colour = "black"),
    
    # Plot title
    plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
    
    # Legend
    legend.position = "none",
    
    # Panel & layout
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 3, scales = "fixed")

main_species_plot_feces

