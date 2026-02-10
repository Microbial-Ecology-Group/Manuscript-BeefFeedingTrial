read_depths <- sample_sums(AMR_data.ps)

# Summary statistics
summary(read_depths)

# More specific calculations
median_read_depth <- median(read_depths) #320988
min_read_depth <- min(read_depths) #4704
max_read_depth <- max(read_depths) #922839
mean_read_depth <- mean(read_depths) #313367
sd_read_depth <- sd(read_depths) #183357

# Print results
cat("Median Read Depth:", median_read_depth, "\n")
cat("Min Read Depth:", min_read_depth, "\n")
cat("Max Read Depth:", max_read_depth, "\n")
cat("Mean Read Depth:", mean_read_depth, "\n")
cat("Standard Deviation of Read Depth:", sd_read_depth, "\n")

# Create a histogram
ggplot(data.frame(Read_Depth = read_depths), aes(x = Read_Depth)) +
  geom_histogram(binwidth = 1000, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Read Depth Distribution", x = "Read Depth", y = "Count") +
  theme_minimal()


#################################################################################
# Subset data by sample type
fecal_ps <- subset_samples(AMR_data.ps, Sample_Type == "Fecal")
beef_ps <- subset_samples(AMR_data.ps, Sample_Type == "Beef sample")

# Extract read depths for each subset
fecal_read_depths <- sample_sums(fecal_ps)
beef_read_depths <- sample_sums(beef_ps)

# Summary statistics for Fecal samples
cat("FECAL SAMPLE METRICS:\n")
cat("Median Read Depth:", median(fecal_read_depths), "\n") #359526.5
cat("Min Read Depth:", min(fecal_read_depths), "\n")#234902
cat("Max Read Depth:", max(fecal_read_depths), "\n") #922839
cat("Mean Read Depth:", mean(fecal_read_depths), "\n") #401861
cat("Standard Deviation of Read Depth:", sd(fecal_read_depths), "\n\n") #137591

# Summary statistics for Beef samples
cat("BEEF SAMPLE METRICS:\n")
cat("Median Read Depth:", median(beef_read_depths), "\n") #92974
cat("Min Read Depth:", min(beef_read_depths), "\n") #4704
cat("Max Read Depth:", max(beef_read_depths), "\n") #201663
cat("Mean Read Depth:", mean(beef_read_depths), "\n") # 93560.66
cat("Standard Deviation of Read Depth:", sd(beef_read_depths), "\n") #46213.26

# Create histograms for Fecal and Beef Read Depths
fecal_plot <- ggplot(data.frame(Read_Depth = fecal_read_depths), aes(x = Read_Depth)) +
  geom_histogram(binwidth = 1000, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Fecal Sample Read Depth Distribution", x = "Read Depth", y = "Count") +
  theme_minimal()

beef_plot <- ggplot(data.frame(Read_Depth = beef_read_depths), aes(x = Read_Depth)) +
  geom_histogram(binwidth = 1000, fill = "red", color = "black", alpha = 0.7) +
  labs(title = "Beef Sample Read Depth Distribution", x = "Read Depth", y = "Count") +
  theme_minimal()

# Display histograms
fecal_plot
beef_plot



#################################
#Deduplicated reads

# Extract OTU table from phyloseq
AMR_abundance <- otu_table(AMR_data.ps)

# Sum deduplicated reads per sample
dedup_reads_per_sample <- rowSums(AMR_abundance)

# Get median deduplicated reads for fecal and beef samples
fecal_dedup_reads <- dedup_reads_per_sample[sample_data(AMR_data.ps)$Sample_Type == "Fecal"]
beef_dedup_reads <- dedup_reads_per_sample[sample_data(AMR_data.ps)$Sample_Type == "Beef sample"]

# Calculate medians
median_fecal_dedup <- median(fecal_dedup_reads)
median_beef_dedup <- median(beef_dedup_reads)

# Print results
cat("Median deduplicated reads (ARGs) in fecal samples:", median_fecal_dedup, "\n")
cat("Median deduplicated reads (ARGs) in beef samples:", median_beef_dedup, "\n")

