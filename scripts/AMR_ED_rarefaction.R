# Plot rarefaction results
## load tidyverse / dplyr ------------
library(dplyr)   # or: library(tidyverse)
library(tidyr)
library(ggplot2)
library(scales) 
library(purrr)      # map_dfr()

## read the two tables ---------------
rarefaction_counts    <- read.table(
  "AMR/HFT_rarefaction_counts.txt",
  sep   = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)


rarefaction_metadata <- read.table('AMR/HFT_Metadata_2025.txt', header=T, sep='\t', row.names = NULL, quote = "")


## inspect the key columns -----------
names(rarefaction_counts)
names(rarefaction_metadata)

## --- left-join metadata onto counts ---
merged_df <- rarefaction_counts %>% 
  left_join(rarefaction_metadata,
            by = c("Sample" = "MatchID"))

str(merged_df)
## check result -----------------------
head(merged_df)

## Checking why we have rarefaction counts for some samples that are not in the metadata
# Step 1: Count how many unique Samples have NA in NonHost_reads
samples_with_NA <- merged_df %>%
  filter(is.na(NonHost_reads)) %>%
  summarise(Samples = list(unique(Sample)),
            num_samples = n_distinct(Sample))

samples_with_NA$Samples[[1]]

# Step 2: Filter out those rows (keep only rows where NonHost_reads is NOT NA)
filtered_df <- merged_df %>%
  filter(!is.na(NonHost_reads)) %>%
  mutate(NonHost_reads = as.numeric(NonHost_reads))

##
### Now start new calculations ####
##
# Calculate read numbers
#filtered_df <- filtered_df %>% 
#  mutate(Rarefied_reads = SubsamplePercent * NonHost_reads)

filtered_df <- filtered_df %>% 
  mutate(Rarefied_reads = as.numeric(SubsamplePercent) / 100 * NonHost_reads)


# set order
unique(filtered_df$Sample_Type)
filtered_df$Sample_Type <- factor(filtered_df$Sample_Type, levels = c("Fecal", "Beef Sample","Rinsate Pool"))

# add 0's
zero_rows <- filtered_df %>% 
  distinct(Sample, Level, .keep_all = TRUE) %>%  # keep one row per group
  mutate(
    SubsamplePercent = 0,     # optional – set to 0 if you keep this column
    Rarefied_reads   = 0,
    Count     = 0
  )
filtered_df <- bind_rows(filtered_df, zero_rows) %>% 
  arrange(Sample, Level, Rarefied_reads)



# plot the figure at the group level
p <- ggplot(filtered_df %>%                         # start from your latest table
              filter(Level == "group"), # Removed the 0 counts
            aes(x = Rarefied_reads,
                y = Count,
                group  = Sample,                        # one line per sample
                colour = Sample_Type)) +                    # colour by Sample_Type status
  geom_line(size = 1) +
  #geom_smooth(method = "loess",        # <─ smoother instead of raw line
  #            se     = FALSE,          #   hide grey confidence band
  #            span   = 0.5,            #   adjust for more / less smoothing
  #            size   = .6) +
  geom_point(size = 1.5) +
  scale_x_continuous(labels = label_comma()) +
  labs(x = "Rarefied nonhost reads",
       y = "Observed Richness",
       colour = "Sample_Type") +
  scale_colour_brewer(palette = "Set1") +          # pick any palette you like
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),            # ↩ remove gridlines
    panel.grid.minor = element_blank()
  )

p

# Save the plot as png
# ggsave("Example_data/NAHMS_rarefaction_by_Sample_Type_forced0.png",           # file name (extension = PNG)
#        plot   = p,                        # plot object to save
#        width  = 7,                        # inches  (adjust as needed)
#        height = 5,                        # inches
#        dpi    = 300)                      # resolution
# 



#
##
## Calculating Percent_of_total_richness  ####
##
#

# Get per-group total richness at SubsamplePercent == 100
richness_at_100 <- filtered_df %>%
  filter(SubsamplePercent == 100) %>%
  mutate(total_richness_observed = Count) %>%
  select(Sample, Level, total_richness_observed )

# Join it back to the full dataset
filtered_df <- filtered_df %>%
  left_join(richness_at_100, by = c("Sample", "Level")) %>%
  mutate(Percent_of_total_richness = Count / total_richness_observed)

#
##
## Code with smooth lines intersecting at 0,0
##
#

## 1 ───────────────────────────────
##   split data > 0  and build the loess curves
smth_df <- filtered_df %>% 
  filter(Level == "group", Rarefied_reads > 0) %>%      # keep >0 for the fit
  group_by(Sample) %>% 
  arrange(Rarefied_reads) %>% 
  nest() %>%                                            # one row per sample
  mutate(                                              # fit loess per sample
    fit  = map(data, ~ loess(Percent_of_total_richness ~ Rarefied_reads,
                             data = .x, span = .5)),
    grid = map(data, ~ tibble(
      Rarefied_reads = seq(min(.x$Rarefied_reads),
                           max(.x$Rarefied_reads),
                           length.out = 200)))
  ) %>% 
  mutate(pred = map2(grid, fit,
                     ~ mutate(.x,
                              Percent_of_total_richness =
                                predict(.y, newdata = .x)))) %>% 
  select(Sample, pred) %>% 
  unnest(pred) %>%                                      # smoothed lines
  left_join(filtered_df %>% distinct(Sample, Sample_Type, Count),  # add Sample_Type colours
            by = "Sample", relationship = "many-to-many" )

## 2 ───────────────────────────────
##   one straight segment  (0,0)  →  first smoothed point  per sample
link_df <- smth_df %>% 
  group_by(Sample) %>% 
  slice_min(order_by = Rarefied_reads, n = 1) %>%       # first point of curve
  mutate(x0 = 0, y0 = 0)                                # add origin

## 3 ───────────────────────────────
##   plot

### Zoomed in plot ####
p3 <- ggplot() +
  geom_segment(data = link_df,
               aes(x = x0, y = y0,
                   xend = Rarefied_reads,
                   yend = Percent_of_total_richness,
                   colour = Sample_Type, group = Sample),
               linewidth = .3) +
  geom_line(data = smth_df,
            aes(x = Rarefied_reads,
                y = Percent_of_total_richness,
                colour = Sample_Type, group = Sample),
            linewidth = .3) +
  #coord_cartesian(xlim = c(0, 50e6)) +
  #scale_y_continuous(breaks = seq(0, 1, .1)) +
  #scale_x_continuous(
  #  breaks       = seq(0, 50e6, 5e6),   # major ticks every 5 000 000
  #  minor_breaks = seq(0, 50e6, 1e6)
  #) +
  labs(x = "Rarefied Nonhost Reads",
       y = "Proportion of Total Observed Richness",
       colour = "Sample_Type") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.length     = unit(4, "pt"),          # length of all ticks
        axis.ticks.length.x   = unit(4, "pt"),          # (override if wanted)
        axis.ticks.x          = element_line(),         # major ticks (default)
        axis.ticks.x.bottom   = element_line())         # make sure bottom axis drawn


p3

ggsave("Figures/rarefaction_by_Sample_Type_percent_richness_clean_lines.png",           # file name (extension = PNG)
       plot   = p3,                        # plot object to save
       width  = 8,                        # inches  (adjust as needed)
       height = 6,                        # inches
       dpi    = 300)  

#### Rarefied counts ####

# p3_counts <- ggplot() +
#   geom_segment(data = link_df,
#                aes(x = x0, y = y0,
#                    xend = Rarefied_reads,
#                    yend = Count_unique,
#                    colour = Sample_Type, group = Sample),
#                linewidth = .5) +
#   geom_line(data = smth_df,
#             aes(x = Rarefied_reads,
#                 y = Count_unique,
#                 colour = Sample_Type, group = Sample),
#             linewidth = .5) +
#   coord_cartesian(xlim = c(0, 50e6)) +
#   scale_x_continuous(
#     breaks       = seq(0, 50e6, 2e6),   # major ticks every 5 000 000
#     minor_breaks = seq(0, 50e6, 1e6)
#   ) +
#   labs(x = "Rarefied nonhost reads (0 – 50,000,000)",
#        y = "Proportion of total observed richness",
#        colour = "Sample_Type") +
#   scale_colour_brewer(palette = "Set1") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.length     = unit(4, "pt"),          # length of all ticks
#         axis.ticks.length.x   = unit(4, "pt"),          # (override if wanted)
#         axis.ticks.x          = element_line(),         # major ticks (default)
#         axis.ticks.x.bottom   = element_line())         # make sure bottom axis drawn
# 
# 
# p3_counts
# 
# ggsave("Example_data/NAHMS_rarefaction_by_Sample_Type_percent_richness_30M_clean_lines.png",           # file name (extension = PNG)
#        plot   = p3_counts,                        # plot object to save
#        width  = 7,                        # inches  (adjust as needed)
#        height = 5,                        # inches
#        dpi    = 300)  




