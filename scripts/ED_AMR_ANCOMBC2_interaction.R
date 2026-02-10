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
clean_fecal_data.ps <- subset_samples(fecal_data.ps, Combined_order != "Final_washout")
clean_fecal_data.ps <- subset_samples(clean_fecal_data.ps, col_time  != "Washout")
clean_fecal_data.ps <- prune_taxa(taxa_sums(clean_fecal_data.ps) > 0, clean_fecal_data.ps)
any(taxa_sums(clean_fecal_data.ps)==0) # QUADRUPLE CHECKING - nope good.


sample_data(clean_fecal_data.ps)$Combined_order <- factor(sample_data(clean_fecal_data.ps)$Combined_order, levels = c("Start_CONV","Midpoint_CONV",
                                                                                                          "End_CONV","Start_RWA","Midpoint_RWA","End_RWA"))
sample_data(clean_fecal_data.ps)$trt_id        
sample_data(clean_fecal_data.ps)$col_time 
sample_data(clean_fecal_data.ps)$Participant_ID <- as.factor(sample_data(clean_fecal_data.ps)$Participant_ID)

# One way is to aggregate to group level first
#group_clean_fecal_data.ps <- tax_glom(clean_fecal_data.ps, taxrank = "Species")
# Remove taxa with no counts
#group_clean_fecal_data.ps <- prune_taxa(taxa_sums(group_clean_fecal_data.ps) > 0, group_clean_fecal_data.ps)

# Trim taxa
#filtered_fecal_data_noSNP_group <- preDA(fecal_data_noSNP_group, min.samples=2, min.reads = 1) # 10% of samples (77/154)
# 203 taxa left
#tax_table(clean_fecal_data.ps)$group
# Change the most specific 
#colnames(phyloseq::tax_table(clean_fecal_data.ps))[4] <- "Species"


# ANCOM_feces_output = ancombc2(data = clean_fecal_data.ps, assay_name = "counts", tax_level = "Species",
#                               fix_formula = "Combined_order + participant_id + treatment_order",
#                               p_adj_method = "holm", 
#                               group = "Combined_order", struc_zero = TRUE, neg_lb = TRUE,verbose = TRUE,
#                               global = TRUE, pairwise = TRUE, n_cl = 3)


ANCOM_feces_output = ancombc2(data = clean_fecal_data.ps, assay_name = "counts", tax_level = "group", 
                                  fix_formula = "trt_id * col_time", rand_formula = "(1 | Participant_ID)",
                                  p_adj_method = "holm", prv_cut = 0.0, lib_cut = 0, s0_perc = 0.0, struc_zero = TRUE, neg_lb = TRUE,
                                  group = "trt_id", alpha = 0.05, n_cl = 2, verbose = TRUE,
                                  global = TRUE, pairwise = FALSE,
                                  dunnet = FALSE, trend = FALSE)

str(ANCOM_feces_output)
write.csv(ANCOM_feces_output$res, "ANCOMBC_AMR_feces_output_nofilter_w2wayInteraction_MM.csv")

##Now Beef vs Feces

# Subset beef, feces, and remove midpoint samples and washout
feces_v_beef_start_end_AMR.ps <- subset_samples(AMR_data.ps, Combined_order != "Midpoint_CONV" &  Combined_order != "Midpoint_RWA" &  Combined_order != "Final_Washout")

# Run model
ANCOM_feces_v_beef_output_AMR = ancombc2(data = feces_v_beef_start_end_AMR.ps, assay_name = "counts", tax_level = "group",
                                     fix_formula = "Combined_order",
                                     p_adj_method = "holm", 
                                     group = "Combined_order", prv_cut = 0.0, lib_cut = 0, s0_perc = 0.0, struc_zero = TRUE, neg_lb = TRUE,
                                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                                     global = FALSE, pairwise = FALSE,
                                     dunnet = FALSE, trend = FALSE)

# output model results
write.csv(ANCOM_feces_v_beef_output_AMR$res, "ANCOMBC_AMR_ASVtoGroup_feces_v_beef_Combined_order_MM.csv")


### Making figure for feces comparisons ####
