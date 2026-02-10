# HFT ANCOMBC MICROBIOME
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
fecal_data.ps <- subset_samples(data_micro.ps, sample_type == "Fecal")
clean_fecal_data.ps <- subset_samples(fecal_data.ps, Combined_order != "Final_Washout")
clean_fecal_data.ps <- prune_taxa(taxa_sums(clean_fecal_data.ps) > 0, clean_fecal_data.ps)
any(taxa_sums(clean_fecal_data.ps)==0) # QUADRUPLE CHECKING - nope good.

# 🔹 AGGLOMERATE TO GENUS (only change #1)
clean_fecal_data.ps <- tax_glom(clean_fecal_data.ps, taxrank = "Genus")

sample_data(clean_fecal_data.ps)$Combined_order <- factor(
  sample_data(clean_fecal_data.ps)$Combined_order,
  levels = c("Start_No_Claim","Midpoint_No_Claim","End_No_Claim",
           "Start_RWA","Midpoint_RWA","End_RWA")
)
sample_data(clean_fecal_data.ps)$participant_id <- as.factor(sample_data(clean_fecal_data.ps)$participant_id)

# 🔹 Change to Genus (only change #2)
colnames(phyloseq::tax_table(clean_fecal_data.ps))[6] <- "Genus"

# 🔹 Update tax_level to Genus (only change #3)
ANCOM_feces_output = ancombc2(data = clean_fecal_data.ps, assay_name = "counts", tax_level = "Genus", 
                              fix_formula = "Combined_order", rand_formula = "(1 | participant_id)",
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
                    .x == 1 & get(gsub("diff_", "passed_ss_", cur_column())) == 0 ~ "blue",
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
  dplyr::mutate(group = gsub(" ", "_", group)) %>%    # 🔹 NEW LINE: standardize names
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
                             "RWA: Start to End")) %>%
  droplevels()

factor_order_Combined_order = c(
  "No Claim: Start to Midpoint", 
  "No Claim: Midpoint to End", 
  "No Claim: Start to End",
  "RWA: Start to Midpoint", 
  "RWA: Midpoint to End", 
  "RWA: Start to End"
)

df_fig_Combined_order_filtered$group <- factor(
  df_fig_Combined_order_filtered$group,
  levels = factor_order_Combined_order
)



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

# Perform the merge with 'taxon' from scaled_counts_feces.long and ' Genus' from clean_fecal_data.melt
BiasAdj_clean_fecal_data.melt <- clean_fecal_data.melt %>%
  dplyr::left_join(scaled_counts_feces.long, by = c("Genus" = "taxon", "Sample"))

BiasAdj_clean_fecal_data.melt

## Quick check on why some taxa don't have BiasAdj_counts
# 118 genera
scaled_counts_feces.long %>%
  distinct(taxon) %>%
  nrow()
# 118 genera with counts >0 
scaled_counts_feces.long %>%
  filter(BiasAdj_counts > 0) %>%
  distinct(taxon) %>%
  nrow()

# 0 counts < 0 
scaled_counts_feces.long %>%
  filter(!(BiasAdj_counts > 0)) %>%
  distinct(taxon) %>%
  nrow()

# 253 total genera
BiasAdj_clean_fecal_data.melt %>%
  distinct(Genus) %>%
  nrow()

# 111 genera with counts > 0
BiasAdj_clean_fecal_data.melt %>%
  filter(BiasAdj_counts > 0) %>%
  distinct(Genus) %>%
  nrow()


# scaled_counts_feces.long %>%
#   filter(taxon =="e387d8b5cf3d6c61e65019accb884109")
# scaled_counts_feces.long %>%
#   filter(taxon =="Mitsuokella")


# Step 1: Assign "Midpoint to End" and other labels
 Genus_avg_counts_feces_by_timepoint <- BiasAdj_clean_fecal_data.melt %>%
  dplyr::group_by(Genus, Combined_order) %>%
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
duplicated_rows <-  Genus_avg_counts_feces_by_timepoint %>%
  dplyr::filter(Combined_order %in% c("End_No_Claim", "End_RWA")) %>%
  dplyr::mutate(Combined_order_label = dplyr::case_when(
    Combined_order == "End_No_Claim" ~ "No Claim: Start to End",
    Combined_order == "End_RWA" ~ "RWA: Start to End",
    TRUE ~ Combined_order_label
  ))

# Step 3: Combine the original and duplicated rows
 Genus_avg_counts_feces_by_timepoint <-  Genus_avg_counts_feces_by_timepoint %>%
  dplyr::bind_rows(duplicated_rows)


# Check the unique values in the final column
unique( Genus_avg_counts_feces_by_timepoint$Combined_order_label)


df_fig_Combined_order_filtered_merged_data <- df_fig_Combined_order_filtered %>%
  left_join( Genus_avg_counts_feces_by_timepoint, by = c("group" = "Combined_order_label", "taxon" = "Genus"))

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


#No data generated, so check numbers to confirm and write it up. 
# n_genera_start <- ntaxa(tax_glom(fecal_data.ps, "Genus"))
# n_genera_start
# 
# summary(ANCOM_feces_output)
# 
# n_total <- nrow(ANCOM_feces_output$res_global)
# n_total
# 
# sig_global  <- sum(ANCOM_feces_output$res_global$diff_abn)
# sig_pairwise <- sum(sapply(ANCOM_feces_output$res_pair, function(x)
#   if(is.data.frame(x)) sum(x$diff_abn, na.rm = TRUE) else 0))
# 
# cat("Total taxa tested:", n_total,
#     "\nSignificant (global):", sig_global,
#     "\nSignificant (pairwise):", sig_pairwise, "\n")
# 
# which_passed <- subset(ANCOM_feces_output$res_global, passed_ss == 1 & diff_abn == 1)
# head(which_passed)







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
  left_join(tax_data_feces, by = c("taxon" = "Genus"),relationship = "many-to-many")  # Assuming ' Genus' is your key column to join

# Reorder the taxon based on Phylum, Class, and  Genus
df_with_taxonomy_feces <- df_with_taxonomy_feces %>%
  arrange(Phylum, Class, taxon) %>%  # Sort the data by Phylum, Class, and  Genus
  dplyr::mutate(Class = factor(Class, levels = rev(unique(Class))))  # Reorder the 'taxon' factor based on the new order


unique(df_with_taxonomy_feces$Class)

df_with_taxonomy_feces$group = factor(df_with_taxonomy_feces$group, levels = factor_order_Combined_order)

is.na(df_with_taxonomy_feces$Phylum)

## These rows have NAs in the phylum column
# Easiest thing would be to remove them, but there's not that many unique taxa
# with errors so I'll just fix them manually below
# df_with_taxonomy_feces %>%
#   filter(is.na(Phylum))

##
### New edit to fix the taxonomy on these weird taxa. For now I fix it manually but we can
# track down the cause of this later. Likely to do with weird characters in the annotations.
##

# target_firmicutes_taxa <- c(
#   "Bacteria_Firmicutes_Clostridia_Christensenellales_Christensenellaceae_uncultured",
#   "Bacteria_Firmicutes_Clostridia_Lachnospirales_Lachnospiraceae_uncultured",
#   "Bacteria_Firmicutes_Clostridia_Oscillospirales_Oscillospiraceae_uncultured",
#   "Bacteria_Firmicutes_Clostridia_Oscillospirales_Ruminococcaceae_uncultured",
#   "Bacteria_Firmicutes_Clostridia_Peptococcales_Peptococcaceae_uncultured"
# )
# 
# df_with_taxonomy_feces_edited <- df_with_taxonomy_feces %>%
#   dplyr::mutate(
#     # work with characters so we don't mangle factor codes
#     Phylum = as.character(Phylum),
#     Class  = as.character(Class),
#     
#     Phylum = case_when(
#       taxon %in% target_firmicutes_taxa ~ "Firmicutes",
#       taxon == "Bacteria_Actinobacteriota_Coriobacteriia_Coriobacteriales_Coriobacteriales_Incertae_Sedis_uncultured" ~ "Actinobacteriota",
#       taxon == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodospirillales_uncultured_uncultured" ~ "Proteobacteria",
#       TRUE ~ Phylum
#     ),
#     
#     Class = case_when(
#       taxon %in% target_firmicutes_taxa ~ "Clostridia",
#       taxon == "Bacteria_Actinobacteriota_Coriobacteriia_Coriobacteriales_Coriobacteriales_Incertae_Sedis_uncultured" ~ "Coriobacteriia",
#       taxon == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodospirillales_uncultured_uncultured" ~ "Alphaproteobacteria",
#       TRUE ~ Class
#     )
#   )
# 




## Trying to add taxonomy
main_Genus_plot_feces <- ggplot(df_with_taxonomy_feces, aes(x = value, y = forcats::fct_rev(Phylum))) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +  # Add vertical dashed line at 0
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),  # Different shapes for each taxon
             position = position_dodge(width = 0.5)) +  # Dodge points within the same class for clarity
  scale_color_identity() +  # Use color directly from the dataset
  scale_size_continuous(range = c(1, 6)) +  # Increase the size range to create more contrast
  labs(x = "Log-Fold Change", y = "Taxon") + #  title = "FECAL - logFC of AMR groups in CEF group over time"
  theme_minimal() +  # Set the theme to minimal
  theme(
    axis.text.y =element_text(size = 10, colour = "black"),  # Hide y-axis  Genus labels (to show in taxonomy plot)
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
main_Genus_plot_feces

ggsave("figures/DA_feces_micro_Genus_by_time.jpg", plot = main_Genus_plot_feces, width = 20, height = 11, dpi = 300)


#Same fig as above, but bigger text to match AMR. 
main_Genus_plot_feces <- ggplot(df_with_taxonomy_feces, aes(x = value, y = forcats::fct_rev(Phylum))) +
  geom_vline(xintercept = 0, color = "grey70", linetype = "solid") +
  geom_point(aes(size = Avg_BiasAdj_Count, color = color),
             position = position_dodge(width = 0.5)) +
  scale_color_identity() +
  scale_size_continuous(range = c(1, 6)) +
  labs(x = "Log-Fold Change", y = NULL) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.y   = element_text(size = 22, colour = "black"),
    axis.text.x   = element_text(size = 24, angle = 45, vjust = 0.5, colour = "black"),
    axis.title.x  = element_text(size = 26, colour = "black"),
    axis.title.y  = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 26, colour = "black"),
    plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~group, ncol = 3, scales = "fixed")

main_Genus_plot_feces


# 
# # Create the taxonomy plot data
# taxonomy_plot_data_feces <- df_with_taxonomy_feces %>%
#   distinct(Phylum, Class, taxon)
# 
# 
# taxonomy_plot_data_feces <- taxonomy_plot_data_feces %>%
#   dplyr::mutate(Phylum = as.character(Phylum)) %>%   # ensure character
#   dplyr::group_by(Phylum) %>%
#   dplyr::mutate(
#     label_Phylum = dplyr::if_else(
#       row_number() == 1,
#       Phylum,
#       ""           # same type: character
#     )
#   ) %>%
#   ungroup()
# 
# # Modify the data to create new columns with the "label_" prefix
# taxonomy_plot_data_feces <- taxonomy_plot_data_feces %>%
#   dplyr::group_by(Phylum) %>%
#   dplyr::mutate(label_Phylum = ifelse(row_number() == 1, Phylum, "")) %>%  # Create 'label_Phylum' with only the first occurrence of each class
#   ungroup() 
# 
# # Create the updated taxonomy plot
# taxonomy_plot_feces <- ggplot(taxonomy_plot_data_feces) +
#   geom_text(aes(x = 1, y = Class, label = label_Phylum), hjust = 1, size = 3.5) +  # Use 'label_Phylum' for class
#   #scale_x_continuous(limits = c(0.5, 2.5),  labels = c("Class")) +
#   scale_y_discrete(limits = levels(taxonomy_plot_data_feces$Class)) +  # Ensure y-axis matches the taxon order
#   theme_void() +
#   theme(
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     plot.margin = unit(c(0, 0, 0, 0), "cm")  # Reduce margin on the right of taxonomy_plot
#   )
# 
# 
# combined_plot_feces  <- plot_grid(
#   taxonomy_plot_feces  + theme(plot.margin = unit(c(1, 0, 0, 0.5), "cm")),  # Remove margins from the taxonomy plot
#   main_Genus_plot_feces  + theme(plot.margin = unit(c(1, 2, 0, 0), "cm")),  # Remove margins from the main plot
#   ncol = 2, 
#   rel_widths = c(0.8,1.8),  # Adjust the width ratio as needed
#   align = "h", 
#   axis = "tb"
# )
# 
# combined_plot_feces
# # Save the final plot
# ggsave("figures/DA_feces_micro_Genus_by_time.jpg", plot = combined_plot_feces, width = 20, height = 11, dpi = 300)
# 
# 
