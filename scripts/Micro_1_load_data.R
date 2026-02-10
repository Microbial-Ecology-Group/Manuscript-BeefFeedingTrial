# import data
qiimedata <- import_biom(
  "16S/table-with-taxonomy-fixed.json",
  "16S/tree.nwk",
  "16S/dna-sequences.fasta"
)



#qiimedata <- import_biom("16S/table-with-taxonomy.json", "16S/tree.nwk", "16S/dna-sequences.fasta")
map_file <- import_qiime_sample_data("16S/HFT_Metadata_9_27_2023.txt") 

sample_names(qiimedata)[1:20]  # look at first 20
rownames(sample_data(map_file))[1:20]

#They dont match, so lets fix them
# Normalize sample names in qiimedata
sample_names(qiimedata) <- gsub("_L001$", "", sample_names(qiimedata))

# Now try merging
data <- merge_phyloseq(qiimedata, map_file)


## DATA EXPLORATION
data # we have 235 samples, 4464 taxa

# extract taxonomy table
tax.data <- data.frame(phyloseq::tax_table(data), stringsAsFactors = FALSE)

# your current data has 7 ranks
colnames(tax.data) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

# fix SILVA-style prefixes (d__, p__, etc)
silva_clean <- function(df) {
  out <- df
  out$Kingdom <- str_replace(out$Kingdom, "d__", "")
  out$Phylum  <- str_replace(out$Phylum,  "p__", "")
  out$Class   <- str_replace(out$Class,   "c__", "")
  out$Order   <- str_replace(out$Order,   "o__", "")
  out$Family  <- str_replace(out$Family,  "f__", "")
  out$Genus   <- str_replace(out$Genus,   "g__", "")
  out$Species <- str_replace(out$Species, "s__", "")
  out
}

tax.data.names <- silva_clean(tax.data)

# convert NAs to blank
tax.data.names[is.na(tax.data.names)] <- ""

# fill unclassified levels *without changing column lengths*
for (i in 1:nrow(tax.data.names)) {
  
  if (tax.data.names$Phylum[i] == "") {
    fill <- paste0("unclassified ", tax.data.names$Kingdom[i])
    tax.data.names[i, 2:7] <- fill
    
  } else if (tax.data.names$Class[i] == "") {
    fill <- paste0("unclassified ", tax.data.names$Phylum[i])
    tax.data.names[i, 3:7] <- fill
    
  } else if (tax.data.names$Order[i] == "") {
    fill <- paste0("unclassified ", tax.data.names$Class[i])
    tax.data.names[i, 4:7] <- fill
    
  } else if (tax.data.names$Family[i] == "") {
    fill <- paste0("unclassified ", tax.data.names$Order[i])
    tax.data.names[i, 5:7] <- fill
    
  } else if (tax.data.names$Genus[i] == "") {
    tax.data.names$Genus[i] <- paste0("unclassified ", tax.data.names$Family[i])
    
  } else if (tax.data.names$Species[i] == "") {
    tax.data.names$Species[i] <- paste0("unclassified ", tax.data.names$Genus[i])
  }
}



uncultured_vals <- c("uncultured", "unculture")

tax.data.names <- tax.data.names %>%
  mutate(
    # make sure we are working with characters
    across(c(Phylum, Class, Order, Family, Genus), as.character),
    
    # choose the first non-uncultured, non-NA label going up the tree
    donor_label = case_when(
      !is.na(Family) & !(Family %in% uncultured_vals) ~ Family,
      !is.na(Order)  & !(Order  %in% uncultured_vals) ~ Order,
      !is.na(Class)  & !(Class  %in% uncultured_vals) ~ Class,
      !is.na(Phylum) & !(Phylum %in% uncultured_vals) ~ Phylum,
      TRUE ~ NA_character_
    ),
    
    # update Genus only when it's uncultured
    Genus = if_else(
      Genus %in% uncultured_vals & !is.na(donor_label),
      paste0(donor_label, "_uncultured"),
      Genus
    )
  ) %>%
  select(-donor_label)  

tax.data.names %>%
  filter(Genus %in% uncultured_vals | grepl("_uncultured$", Genus)) %>%
  select(Phylum, Class, Order, Family, Genus) %>%
  distinct()


head(tax.data.names) # great, no more NAs and no more k__
tax_table(data) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object

#tax_table(data) <- tax_fix(tax_table(data))

tail(tax_table(data), 20) # sweet, lookin good!




##### Clean up sample types ####

# lets split up the controls and samples
controls <- phyloseq::subset_samples(data, Combined_order=="Blank")
controls # 19 samples, 4464 taxa

# Keep only the samples (non-blanks)
samples <- phyloseq::subset_samples(data, Combined_order!="Blank")
any(sample_sums(samples)==0)
sum(taxa_sums(samples)==0) # 95 taxa gone
samples <- prune_taxa(taxa_sums(samples) > 0, samples)
samples # 216 samples, 4369  taxa left

## # lets look at the number of reads per sample and the distribution
sample_sum_df <- data.frame(sum = sample_sums(samples))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) # pretty good

# some QC checks
min(sample_sums(samples)) # 4095
max(sample_sums(samples)) # 518298
sort(sample_sums(samples)) # all good


samples <- prune_samples(sample_sums(samples) > 100, samples)
samples #

#data_micro.ps <- subset_taxa(samples, Domain!="Eukaryota")
#Above line not needed since we switched classifiers. Rename data object below
data_micro.ps <- samples 
data_micro.ps  # 216  samples, 4369   taxa left

# Check sample sums
max(sample_sums(data_micro.ps))
sort(sample_sums(data_micro.ps))


sample_data(data_micro.ps)$Sample_Type

#saveRDS(data_micro.ps, file = "C:/Users/mcclu/Desktop/data_micro_ps.rds")

