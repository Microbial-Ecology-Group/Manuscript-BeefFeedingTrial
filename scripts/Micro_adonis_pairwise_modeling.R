library(vegan)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

# Distance and metadata you already have:
# Fecal_data_micro.ps.dist
# Fecal_data_micro.ps.df

# Make sure timepoint is a factor
Fecal_data_micro.ps.df$timepoint <- factor(Fecal_data_micro.ps.df$timepoint)

# All unique timepoints
tp_levels <- levels(Fecal_data_micro.ps.df$timepoint)

# All pairwise combinations of timepoints
tp_pairs <- combn(tp_levels, 2, simplify = FALSE)

pairwise_timepoint <- map_dfr(tp_pairs, function(tp_pair) {
  # Subset metadata to the two timepoints
  idx  <- Fecal_data_micro.ps.df$timepoint %in% tp_pair
  meta_sub <- droplevels(Fecal_data_micro.ps.df[idx, ])
  
  # Subset distance matrix to those samples
  dist_mat  <- as.matrix(Fecal_data_micro.ps.dist)
  dist_sub  <- dist_mat[idx, idx]
  dist_sub  <- as.dist(dist_sub)
  
  # PERMANOVA with same structure (trt + timepoint + treatment_order),
  # but only 2 levels of timepoint now; keep strata for repeated measures
  ad <- adonis2(
    dist_sub ~ trt + timepoint + treatment_order,
    data   = meta_sub,
    strata = meta_sub$participant_id,
    permutations = 9999,
    by = "margin"
  )
  
  # Extract the row for timepoint
  row_tp <- ad["timepoint", , drop = FALSE]
  
  tibble(
    timepoint1 = tp_pair[1],
    timepoint2 = tp_pair[2],
    Df         = row_tp$Df,
    SumOfSqs   = row_tp$SumOfSqs,
    R2         = row_tp$R2,
    F          = row_tp$F,
    p          = row_tp$`Pr(>F)`
  )
})

# Optional: adjust for multiple testing
pairwise_timepoint <- pairwise_timepoint %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

pairwise_timepoint

###
####
### By treatment ####
####
###


trt_levels <- levels(as.factor(Fecal_data_micro.ps.df$trt))

pairwise_timepoint_by_trt <- map_dfr(trt_levels, function(tt) {
  # subset to one treatment
  df_tt <- Fecal_data_micro.ps.df %>% filter(trt == tt)
  idx   <- rownames(Fecal_data_micro.ps.df) %in% rownames(df_tt)
  
  dist_mat <- as.matrix(Fecal_data_micro.ps.dist)
  dist_tt  <- as.dist(dist_mat[idx, idx])
  
  # timepoint pairs within this treatment
  tp_levels <- levels(droplevels(df_tt$timepoint))
  tp_pairs  <- combn(tp_levels, 2, simplify = FALSE)
  
  map_dfr(tp_pairs, function(tp_pair) {
    sel <- df_tt$timepoint %in% tp_pair
    meta_sub <- droplevels(df_tt[sel, ])
    
    dist_mat_sub <- as.matrix(dist_tt)
    dist_sub     <- as.dist(dist_mat_sub[sel, sel])
    
    ad <- adonis2(
      dist_sub ~ timepoint + treatment_order,
      data   = meta_sub,
      strata = meta_sub$participant_id,
      permutations = 9999,
      by = "margin"
    )
    
    row_tp <- ad["timepoint", , drop = FALSE]
    
    tibble(
      trt        = tt,
      timepoint1 = tp_pair[1],
      timepoint2 = tp_pair[2],
      Df         = row_tp$Df,
      SumOfSqs   = row_tp$SumOfSqs,
      R2         = row_tp$R2,
      F          = row_tp$F,
      p          = row_tp$`Pr(>F)`
    )
  })
})

pairwise_timepoint_by_trt <- pairwise_timepoint_by_trt %>%
  group_by(trt) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup()

pairwise_timepoint_by_trt