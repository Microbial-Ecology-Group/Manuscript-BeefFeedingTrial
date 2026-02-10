# HFT ANCOMBC
set.seed(123)
library(ANCOMBC)
library(microbiome)
library(dplyr)
library(tidyr)

###
##### Start ANCOMBC - Fecal samples, comparing by Combined_order #####
###
fecal_data.ps <- subset_samples(AMR_data.ps, Sample_Type == "Fecal")


sample_data(fecal_data.ps)$Combined_order <- factor(sample_data(fecal_data.ps)$Combined_order, levels = c("Start_CONV","Midpoint_CONV",
                                                                                                                "End_CONV","Start_RWA","Midpoint_RWA",
                                                                                                                "End_RWA","Final_washout"))
sample_data(fecal_data.ps)$Participant_ID <- as.factor(sample_data(fecal_data.ps)$Participant_ID)

# Trim taxa
#filtered_fecal_data_noSNP_group <- preDA(fecal_data_noSNP_group, min.samples=2, min.reads = 1) # 10% of samples (77/154)
# 203 taxa left

# Change the most specific 
colnames(phyloseq::tax_table(fecal_data.ps))[4] <- "Species"


ANCOM_feces_output = ancombc2(data = fecal_data.ps, assay_name = "counts", tax_level = "Species",
                              fix_formula = "Combined_order + Participant_ID + blinded_treatment",
                              p_adj_method = "holm", 
                              group = "Combined_order", struc_zero = TRUE, neg_lb = TRUE,verbose = TRUE,
                              global = TRUE, pairwise = TRUE, n_cl = 2)

# Extract output for pairwise comparisons
res_pair = ANCOM_feces_output$res_pair


#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_Combined_order = res_pair %>%
  dplyr::select(taxon, contains("Combined_order")) 

df_fig_Combined_order = df_Combined_order %>%
  filter(diff_Combined_orderMidpoint_CONV == 1 | diff_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV == 1 
         | diff_Combined_orderStart_RWA_Combined_orderEnd_CONV  == 1 | diff_Combined_orderMidpoint_RWA_Combined_orderStart_RWA  == 1 
         | diff_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA == 1 | diff_Combined_orderFinal_washout_Combined_orderEnd_RWA == 1 ) %>%
  mutate(lfc_Midpoint_CONV = ifelse(diff_Combined_orderMidpoint_CONV == 1,  # This compares the 2nd time point to Start_CONV
                        lfc_Combined_orderMidpoint_CONV, 0),
         lfc_End_CONV = ifelse(diff_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV== 1, 
                        lfc_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV, 0),
         lfc_Start_RWA = ifelse(diff_Combined_orderStart_RWA_Combined_orderEnd_CONV== 1, 
                        lfc_Combined_orderStart_RWA_Combined_orderEnd_CONV, 0),
         lfc_Midpoint_RWA = ifelse(diff_Combined_orderMidpoint_RWA_Combined_orderStart_RWA== 1, 
                        lfc_Combined_orderMidpoint_RWA_Combined_orderStart_RWA, 0),
         lfc_End_RWA = ifelse(diff_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA== 1, 
                        lfc_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA, 0),
         lfc_Final_washout = ifelse(diff_Combined_orderFinal_washout_Combined_orderEnd_RWA== 1, 
                        lfc_Combined_orderFinal_washout_Combined_orderEnd_RWA, 0)
         ) %>%
  transmute(taxon, 
            `Midpoint_CONV vs. Start_CONV` = round(lfc_Midpoint_CONV, 2),
            `End_CONV vs. Midpoint_CONV` = round(lfc_End_CONV, 2),
            `Start_RWA vs. End_CONV` = round(lfc_Start_RWA, 2),
            `Midpoint_RWA vs. Start_RWA` = round(lfc_Midpoint_RWA, 2),
            `End_RWA vs. Midpoint_RWA` = round(lfc_End_RWA, 2),
            `Final_wahsout vs. End_RWA` = round(lfc_Final_washout, 2)) %>%
  pivot_longer(cols = `Midpoint_CONV vs. Start_CONV`:`End_CONV vs. Midpoint_CONV`:`Start_RWA vs. End_CONV`:`Midpoint_RWA vs. Start_RWA`:`End_RWA vs. Midpoint_RWA`:`Final_wahsout vs. End_RWA`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)


### Heatmap plot ####
# Calculate LFC value ranges
lo = floor(min(df_fig_Combined_order$value))
up = ceiling(max(df_fig_Combined_order$value))
mid = (lo + up)/2



fig_pair = df_fig_Combined_order %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(lo, up),
                       na.value = "white", name = NULL) +
  geom_text(aes(group, taxon, label = ifelse(value == 0, "", value)), size = 6) + # don't show the 0s , #, color = color
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold changes compared to NEG - Day 7") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_pair







###
##### Experimental part of the code, not totally done #####
###
Arm1_TrtA_data_noSNP <- subset_samples(data_noSNP, arm == "1" & blinded_treatment == "A")
sample_data(Arm1_TrtA_data_noSNP)$pool_timepoint <- factor(sample_data(Arm1_TrtA_data_noSNP)$pool_timepoint, levels = c("A","B","C"))

Arm1_TrtA_data_noSNP_group <- tax_glom(Arm1_TrtA_data_noSNP, taxrank = "group")
Arm1_TrtA_data_noSNP_group <- prune_taxa(taxa_sums(Arm1_TrtA_data_noSNP_group) > 0, Arm1_TrtA_data_noSNP_group)
taxa_names(Arm1_TrtA_data_noSNP_group) <- as.data.frame(tax_table(Arm1_TrtA_data_noSNP_group))$Group


output = ancombc2(data = Arm1_TrtA_data_noSNP_group, assay_name = "counts", tax_level = NULL,
                  fix_formula = "pool_timepoint + Participant_ID",
                  rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "pool_timepoint", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 1, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 10, verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 10),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 1),
                  trend_control = list(contrast =
                                         list(matrix(c(1, 0, -1, 1),
                                                     nrow = 2,
                                                     byrow = TRUE)),
                                       node = list(2),
                                       solver = "ECOS",
                                       B = 1))
res_prim = output$res
tab_sens = output$pseudo_sens_tab
res_pair = output$res_pair

res_global = output$res_global
res_dunn = output$res_dunn
res_trend = output$res_trend

#### Calculate sensitivity scores for pool_timepoint ####
## Rename pairwise comparisons
df_bmi = res_pair %>%
  dplyr::select(taxon, contains("timepoint")) 
df_fig_bmi = df_bmi %>%
  filter(diff_pool_timepointB == 1 | diff_pool_timepointC == 1) %>%
  mutate(lfc_B = ifelse(diff_pool_timepointB == 1, 
                        lfc_pool_timepointB, 0),
         lfc_C = ifelse(diff_pool_timepointC== 1, 
                        lfc_pool_timepointC, 0)) %>%
  transmute(taxon, 
            `B vs. A` = round(lfc_B, 2),
            `C vs. A` = round(lfc_C, 2)) %>%
  pivot_longer(cols = `B vs. A`:`C vs. A`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Plot heatmap
# This one doesn't have a mid point since there were only negative fold changes
lo = floor(min(df_fig_timepoint$value))
up = ceiling(max(df_fig_timepoint$value))
mid = (lo + up)/2
fig_timepoint = df_fig_timepoint %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "white", 
                       na.value = "white", limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to the start of the trial") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_timepoint

# Plot heatmap
lo = floor(min(df_fig_timepoint$value))
up = ceiling(max(df_fig_timepoint$value))
mid = (lo + up)/2
fig_timepoint = df_fig_timepoint %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_timepoint





###
#### To run with longitudinal analysis, need global=FALSE so the model can be fit separately for each time point
###


# Remove missing values
fecal_data_noSNP_group <- na.omit(fecal_data_noSNP_group)

# Run ANCOMBC2 with longitudinal analysis
ANCOM_feces_output <- ancombc2(data = fecal_data_noSNP_group, assay_name = "counts", tax_level = NULL,
                               fix_formula = "Combined_order + blinded_treatment + participant_id",
                               rand_formula = NULL,
                               p_adj_method = "holm", pseudo = 0, pseudo_sens = FALSE,
                               group = "Combined_order", struc_zero = TRUE, neg_lb = TRUE, verbose = TRUE,
                               global = TRUE, pairwise = TRUE)

summary(ANCOM_feces_output)

# Extract results from ANCOM output
ANCOM_feces_results <- ANCOM_feces_output$res

# Filter results for significant comparisons
ANCOM_feces_sig <- ANCOM_feces_results %>% filter("p_(Intercept)" <= 0.05)


# extract the ANCOM output table, if empty it means there were no significant differences
ANCOM_table <- as.data.frame(ANCOM_feces_output$res)



# Get pairwise comparison results
pairwise_res <- ANCOM_feces_output$res_pair

# Convert to data.frame
pairwise_df <- as.data.frame(pairwise_res)

# Create barplot
ggplot2::ggplot(pairwise_df, ggplot2::aes(x = factor(taxon), y = Estimate, fill = Comparison)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge") +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "Feature", y = "Estimate") +
  ggplot2::ggtitle("Pairwise comparisons for Combined_order")

library(reshape2)
library(dplyr)
library(tidyverse)

res_pair <- ANCOM_feces_output$res_pair
## Rename pairwise comparisons
df_Combined_order = res_pair %>%
  dplyr::select(taxon, contains("Combined_order")) 
df_fig_Combined_order = df_Combined_order %>%
  filter(diff_Combined_orderMidpoint_CONV == 1 | diff_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV == 1 
         | diff_Combined_orderStart_RWA_Combined_orderEnd_CONV  == 1 | diff_Combined_orderMidpoint_RWA_Combined_orderStart_RWA  == 1 
         | diff_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA == 1 | diff_Combined_orderFinal_washout_Combined_orderEnd_RWA == 1 ) %>%
  mutate(lfc_Midpoint_CONV = ifelse(diff_Combined_orderMidpoint_CONV == 1,  # This compares the 2nd time point to Start_CONV
                                    lfc_Combined_orderMidpoint_CONV, 0),
         lfc_End_CONV = ifelse(diff_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV== 1, 
                               lfc_Combined_orderEnd_CONV_Combined_orderMidpoint_CONV, 0),
         lfc_Start_RWA = ifelse(diff_Combined_orderStart_RWA_Combined_orderEnd_CONV== 1, 
                                lfc_Combined_orderStart_RWA_Combined_orderEnd_CONV, 0),
         lfc_Midpoint_RWA = ifelse(diff_Combined_orderMidpoint_RWA_Combined_orderStart_RWA== 1, 
                                   lfc_Combined_orderMidpoint_RWA_Combined_orderStart_RWA, 0),
         lfc_End_RWA = ifelse(diff_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA== 1, 
                              lfc_Combined_orderEnd_RWA_Combined_orderMidpoint_RWA, 0),
         lfc_Final_washout = ifelse(diff_Combined_orderFinal_washout_Combined_orderEnd_RWA== 1, 
                                    lfc_Combined_orderFinal_washout_Combined_orderEnd_RWA, 0)
  ) %>%
  transmute(taxon, 
            `Midpoint_CONV vs. Start_CONV` = round(lfc_Midpoint_CONV, 2),
            `End_CONV vs. Midpoint_CONV` = round(lfc_End_CONV, 2),
            `Start_RWA vs. End_CONV` = round(lfc_Start_RWA, 2),
            `Midpoint_RWA vs. Start_RWA` = round(lfc_Midpoint_RWA, 2),
            `End_RWA vs. Midpoint_RWA` = round(lfc_End_RWA, 2),
            `Final_wahsout vs. End_RWA` = round(lfc_Final_washout, 2)) %>%
  pivot_longer(cols = `Midpoint_CONV vs. Start_CONV`:`End_CONV vs. Midpoint_CONV`:`Start_RWA vs. End_CONV`:`Midpoint_RWA vs. Start_RWA`:`End_RWA vs. Midpoint_RWA`:`Final_wahsout vs. End_RWA`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

# Plot heatmap
# This one doesn't have a mid point since there were only negative fold changes
lo = floor(min(df_fig_Combined_order$value))
up = ceiling(max(df_fig_Combined_order$value))
mid = (lo + up)/2
fig_timepoint = df_fig_Combined_order %>%
  ggplot2::ggplot(ggplot2::aes(x = group, y = taxon, fill = value)) + 
  ggplot2::geom_tile(color = "black") +
  ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  ggplot2::geom_text(ggplot2::aes(group, taxon, label = value), color = "black", size = 4) +
  ggplot2::labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
fig_timepoint


ANCOM_feces_melted <- melt(as.data.frame(ANCOM_feces_results), id.vars=c("taxon", "lfc_participant_id"))

ANCOM_feces_melted <- ANCOM_feces_melted[order(ANCOM_feces_melted$taxon, ANCOM_feces_melted$variable),]

colnames(ANCOM_feces_melted) <- c("taxon", "participant_id", "variable", "value")

# Pivot data into long format
ANCOM_feces_long <- ANCOM_feces_sig %>% 
  pivot_longer(cols = c(starts_with("diff_"), starts_with("lfc_")), 
               names_to = "metric", 
               values_to = "value") %>% 
  separate(metric, into = c("metric", "group"), sep = "_", remove = FALSE) %>% 
  select(-c(metric, "se_", "W_", "p_", "q_")) %>% 
  group_by(participant_id, Combined_order, group) %>% 
  summarize(value = mean(value), .groups = "drop")

# Plot results
ggplot(ANCOM_feces_long, aes(x = Combined_order, y = value, color = group)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~participant_id, nrow = 2) +
  theme_minimal()