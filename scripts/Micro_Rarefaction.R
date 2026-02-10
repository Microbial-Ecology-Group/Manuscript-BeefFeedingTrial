#Rarefaction

#USE DIFFERENT PS OBJECT FOR AMR/MICRO/SNV
#AMR- data_micro.ps
#Micro- data_micro.ps
#SNV- 
#Not using SNV for publication, so not necessary per PM 7/16/2024


# Extract the OTU table as a matrix
otu_mat <- as(otu_table(data_micro.ps), "matrix")

# Ensure that samples are rows and OTUs are columns
if (taxa_are_rows(data_micro.ps)) {
  otu_mat <- t(otu_mat)
}

# Subset Fecal Samples
fecal_subset <- subset_samples(data_micro.ps, sample_type == "Fecal")
otu_mat_fecal <- as(otu_table(fecal_subset), "matrix")

# Ensure that samples are rows and OTUs are columns for the subset
if (taxa_are_rows(fecal_subset)) {
  otu_mat_fecal <- t(otu_mat_fecal)
}

# Compute the rarefaction curves for fecal samples
rarecurves_fecal <- rarecurve(otu_mat_fecal, step = 50, sample = min(rowSums(otu_mat_fecal)), label = FALSE)

# Initialize an empty data frame to collect the results for fecal samples
rarecurve_df_fecal <- data.frame(Sample = character(), Reads = integer(), Species = numeric(), stringsAsFactors = FALSE)

# Extract rarefaction curve data for fecal samples
for (i in seq_along(rarecurves_fecal)) {
  curve <- rarecurves_fecal[[i]]
  subsamples <- attr(curve, "Subsample")
  species_counts <- as.numeric(curve)
  sample_name <- if (is.null(names(rarecurves_fecal)[i])) paste("Sample", i) else names(rarecurves_fecal)[i]
  
  # Create a temporary data frame for the current sample
  temp_df <- data.frame(Sample = rep(sample_name, length(subsamples)), Reads = subsamples, Species = species_counts)
  # Combine with the main data frame
  rarecurve_df_fecal <- rbind(rarecurve_df_fecal, temp_df)
}

# Apply a cutoff to filter reads up to 150,000
#Four samples had more than 150,000 reads
cutoff_reads <- 150000
rarecurve_df_fecal_filtered <- rarecurve_df_fecal[rarecurve_df_fecal$Reads <= cutoff_reads, ]

# Create the rarefaction curve plot for fecal samples with cutoff
Micro_fecal <- ggplot(rarecurve_df_fecal_filtered, aes(x = Reads, y = Species, group = Sample, color = Sample)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Fecal Samples",
    x = "Number of Reads",
    y = "Unique ASVs"
  ) +
  theme(legend.position = "none")
Micro_fecal




#Now Beef Samples

# Subset Beef Samples
beef_subset <- subset_samples(data_micro.ps, sample_type == "Beef sample")
otu_mat_beef <- as(otu_table(beef_subset), "matrix")

# Ensure that samples are rows and OTUs are columns for the subset
if (taxa_are_rows(beef_subset)) {
  otu_mat_beef <- t(otu_mat_beef)
}

# Compute the rarefaction curve
rarecurves_beef <- rarecurve(otu_mat_beef, step = 50, sample = min(rowSums(otu_mat_beef)), label = FALSE)

#empty data frame to collect the results for beef samples
rarecurve_df_beef <- data.frame(Sample = character(), Reads = integer(), Species = numeric(), stringsAsFactors = FALSE)

# Extract rarefaction curve data
for (i in seq_along(rarecurves_beef)) {
  curve <- rarecurves_beef[[i]]
  subsamples <- attr(curve, "Subsample")
  species_counts <- as.numeric(curve)
  sample_name <- if (is.null(names(rarecurves_beef)[i])) paste("Sample", i) else names(rarecurves_beef)[i]
  
  # Create a temporary data frame
  temp_df <- data.frame(Sample = rep(sample_name, length(subsamples)), Reads = subsamples, Species = species_counts)
  # Combine with the main data frame
  rarecurve_df_beef <- rbind(rarecurve_df_beef, temp_df)
}

# Apply a cutoff to filter reads up to 2e5
cutoff_reads <- 150000
rarecurve_df_beef_filtered <- rarecurve_df_beef[rarecurve_df_beef$Reads <= cutoff_reads, ]

# Create the rarefaction curve plot for beef samples with cutoff
Micro_beef <- ggplot(rarecurve_df_beef_filtered, aes(x = Reads, y = Species, group = Sample, color = Sample)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Beef Samples",
    x = "Number of Reads",
    y = "Unique ASVs"
  ) +
  theme(legend.position = "none")
Micro_beef

#Now, plot together in a panel
#Combine plots into panel for publication
# Combine the "Micro_fecal" and "Micro_beef" plots side by side without the legend
micro_combined <- plot_grid(
  Micro_fecal + theme(legend.position = "none"), 
  Micro_beef + theme(legend.position = "none"), 
  ncol = 2
)

# Add a title above the combined plots
micro_final <- plot_grid(
  ggdraw() + draw_label("Microbiome Rarefaction Curves", size = 18),
  micro_combined,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Display the final plot
micro_final


#Secondary curve with individual sample labels because the beef plot has some interesting curves
Micro_beef <- ggplot(rarecurve_df_beef, aes(x = Reads, y = Species, group = Sample, color = Sample)) +
  geom_line() +
  geom_text_repel(aes(label = Sample), 
                  size = 3) +
  theme_minimal() +
  labs(title = "Microbiome Rarefaction Curve (Beef Samples)", 
       x = "Number of Reads", 
       y = "Unique Species") +
  theme(legend.position = "none")
Micro_beef
