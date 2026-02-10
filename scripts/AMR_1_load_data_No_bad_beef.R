#New Load data, removes two low count beef samples

# Load Lee's Data object
#load("AMR/HFT_AMR-TE.RData")

# #### TE DATA ####
# read in AMR count matrix for V3 results with sNP confirmation
amr_data <- read.table('AMR/deduped_SNPconfirmed_AMR_analytic_matrix.csv', header=T, row.names=1, sep=',', quote = "")
# convert this into an 'otu_table' format required for phyloseq
amr_data <- otu_table(amr_data, taxa_are_rows = T)

#read in gene annotations
annotations <- read.table('AMR/megares_annotations_v3.00.csv', header=T, row.names=1, sep=",", quote = "")
#convert this into a 'taxonomy table' format required for phyloseq
annotations <- phyloseq::tax_table(as.matrix(annotations))

# read in TE metadata
amr_metadata <- read.table('AMR/HFT_Metadata_2025.txt', header=T, sep='\t', row.names = 1, quote = "")

###############################################################################
### FINAL NAME CORRECTIONS — EXACT 7 SAMPLES
###############################################################################
### 1. Fix the 5 misnumbered fecal samples (HFT-F → correct HFTF_TE)

fix_map <- c(
  "HFT-F-101A_S3" = "HFTF_TE_101A_S1",
  "HFT-F-101B_S4" = "HFTF_TE_101B_S2",
  "HFT-F-101C_S5" = "HFTF_TE_101C_S3",
  "HFT-F-101D_S6" = "HFTF_TE_101D_S4",
  "HFT-F-102A_S7" = "HFTF_TE_102A_S5"
)

for(old in names(fix_map)){
  if(old %in% rownames(amr_metadata)){
    rownames(amr_metadata)[rownames(amr_metadata) == old] <- fix_map[old]
  }
}

### 2. Convert dot-style TE names → underscore style
rownames(amr_metadata) <- gsub("^HFTF\\.TE\\.", "HFTF_TE_", rownames(amr_metadata))

### 3. Ensure *all* periods → underscores anywhere else
rownames(amr_metadata) <- gsub("\\.", "_", rownames(amr_metadata))
colnames(amr_data)     <- gsub("\\.", "_", colnames(amr_data))

### 4. Add the 2 missing beef samples (CG5 & CG1)

# Pick a beef template row (any HFTR_TE_ row)
beef_template_id <- grep("^HFTR_TE_", rownames(amr_metadata), value = TRUE)[1]
beef_template <- amr_metadata[beef_template_id, ]

# Add missing beef samples
missing_beef <- c("HFTR_TE_CG5_S179", "HFTR_TE_CG1_S163")

for(id in missing_beef){
  amr_metadata[id, ] <- beef_template
  amr_metadata[id, "Sample_Type"] <- "Beef Sample"
  amr_metadata[id, "treatment"] <- "Conventional"   # CG = Conventional beef
}

### 5. Check intersection count
length(intersect(rownames(amr_metadata), colnames(amr_data)))

### Convert metadata to sample_data BEFORE merge
amr_metadata <- sample_data(amr_metadata)

# merge the annotations, the count matrix, and metadata into a phyloseq object
AMR_data.ps <- merge_phyloseq(amr_data, annotations, amr_metadata)

# Remove Rinsate Pools
AMR_data.ps <- subset_samples(AMR_data.ps, Sample_Type != "Rinsate Pool")
AMR_data.ps # 223 samples

# Split off genes requiring SNP confirmation
SNPconfirm <- subset_taxa(AMR_data.ps, snp=="RequiresSNPConfirmation")
SNPconfirm # 319 of the 2912 genes require SNP confirmation

data_noSNP <- subset_taxa(AMR_data.ps, snp!="RequiresSNPConfirmation")
data_noSNP # 2593 genes remain

# Remove Final_washout samples (not used in analysis)
data_noSNP <- subset_samples(data_noSNP, Combined_order != "Final_Washout")
data_noSNP


###############################################################################
### FILTER OUT LOW-DEPTH *BEEF* SAMPLES ONLY (<1000 reads)
###############################################################################

# Identify beef samples
beef_samples <- sample_names(data_noSNP)[sample_data(data_noSNP)$Sample_Type == "Beef Sample"]

# Among beef samples, find those with low read counts
low_beef <- beef_samples[sample_sums(data_noSNP)[beef_samples] < 1000]

# Print to confirm
low_beef

# Remove ONLY those low-read beef samples from data_noSNP
data_noSNP <- prune_samples(!(sample_names(data_noSNP) %in% low_beef), data_noSNP)

# OPTIONAL: remove from full AMR_data.ps as well for consistency
AMR_data.ps <- prune_samples(!(sample_names(AMR_data.ps) %in% low_beef), AMR_data.ps)

# Re-check sample sums
sort(sample_sums(data_noSNP))
