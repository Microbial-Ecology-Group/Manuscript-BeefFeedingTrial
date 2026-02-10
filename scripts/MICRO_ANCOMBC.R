# Microbiome ANCOMBC
library(ANCOMBC)
#starting object
data_micro.ps




#
##
### Feces - full model
##
#

# Subset beef, feces, and remove midpoint samples and washout
#feces_nowashout_micro.ps <- subset_samples(data_micro.ps, Combined_order != "Final_Washout" & sample_type != "Beef sample", participant_id != "135", participant_id != "115", participant_id != "117", participant_id != "118")
feces_nowashout_micro.ps <- subset_samples(data_micro.ps, Combined_order != "Final_Washout" & sample_type != "Meat Rinsate" & pool_timepoint != "G")

#feces_nowashout_micro.ps <- subset_samples(data_micro.ps, Combined_order != "Final_Washout" & sample_type != "Beef sample" & pool_timepoint != "G" &participant_id != "101" & participant_id != "113" &participant_id != "114" &participant_id != "115" & participant_id != "117" & participant_id != "118" &participant_id != "122" &participant_id != "123" &participant_id != "135")

#View(sample_data(feces_nowashout_micro.ps))

feces_nowashout_micro.ps <- prune_taxa(taxa_sums(feces_nowashout_micro.ps) > 0, feces_nowashout_micro.ps)


sample_data(feces_nowashout_micro.ps)$participant_id <- as.character(sample_data(feces_nowashout_micro.ps)$participant_id)

colnames(sample_data(feces_nowashout_micro.ps))
sample_data(feces_nowashout_micro.ps)$trt
sample_data(feces_nowashout_micro.ps)$treatment_order
sample_data(feces_nowashout_micro.ps)$timepoint
sample_data(feces_nowashout_micro.ps)$participant_id


#maybe try prev cut reduction
ANCOM_feces_output.prev0.13_ASV_2way_norand = ancombc2(data = feces_nowashout_micro.ps, assay_name = "counts", tax_level = "Genus",
                                                          fix_formula = "Fecal_Treatment * col_time",  rand_formula = "(1 | participant_id)",
                                                          p_adj_method = "holm", 
                                                          group = "Fecal_Treatment", prv_cut = 0.13, lib_cut = 0, s0_perc = 0.0, struc_zero = TRUE, neg_lb = TRUE,
                                                          alpha = 0.05, n_cl = 4, verbose = TRUE,
                                                          global = FALSE, pairwise = FALSE,
                                                          dunnet = FALSE, trend = FALSE)

write.csv(ANCOM_feces_output.prev0.13_ASV_2way_norand$res, "ANCOMBC_microbiome_feces_ASVtoGenus_prev0.13_2way_wRandom_MM.csv")

## Testing the interaction

ANCOM_micro_feces_output_full_interaction = ancombc2(data = feces_nowashout_micro.ps, assay_name = "counts", tax_level = "Genus",
                                                       fix_formula = "trt * timepoint * treatment_order",  rand_formula = "(1 | participant_id)",
                                                       p_adj_method = "holm", 
                                                       group = "trt", prv_cut = 0.13, lib_cut = 0, s0_perc = 0.0, struc_zero = TRUE, neg_lb = TRUE,
                                                       alpha = 0.05, n_cl = 4, verbose = TRUE,
                                                       global = FALSE, pairwise = FALSE,
                                                       dunnet = FALSE, trend = FALSE)

res_primary <- ANCOM_micro_feces_output_full_interaction$res

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


#
##
### Feces v Beef - Combined order 
##
#

# Subset beef, feces, and remove midpoint samples and washout
feces_v_beef_start_end_micro.ps <- subset_samples(data_micro.ps, Combined_order != "Midpoint_CONV" &  Combined_order != "Midpoint_RWA" &  Combined_order != "Final_Washout")

# Run model
ANCOM_feces_v_beef_output = ancombc2(data = feces_v_beef_start_end_micro.ps, assay_name = "counts", tax_level = "Genus",
                              fix_formula = "Combined_order",
                              p_adj_method = "holm", 
                              group = "Combined_order", prv_cut = 0.0, lib_cut = 0, s0_perc = 0.0, struc_zero = TRUE, neg_lb = TRUE,
                              alpha = 0.05, n_cl = 2, verbose = TRUE,
                              global = FALSE, pairwise = FALSE,
                              dunnet = FALSE, trend = FALSE)

# output model results
write.csv(ANCOM_feces_v_beef_output$res, "ANCOMBC_microbiome_ASVtoGenus_feces_v_beef_Combined_order_MM.csv")

