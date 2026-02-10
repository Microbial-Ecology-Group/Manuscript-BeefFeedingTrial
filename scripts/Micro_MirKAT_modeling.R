# Packages
library(phyloseq)
library(vegan)
library(MiRKAT)        # for D2K to convert distance -> kernel
# remotes::install_github("hk1785/GLMM-MiRKAT")
library(GLMMMiRKAT)


Fecal_data_micro.ps <- subset_samples(
  data_micro.ps,
  sample_type == "Fecal" & !(timepoint %in% c("Final_Washout")) # ,"Start"
)

# (optional) drop zero-sum taxa left behind
Fecal_data_micro.ps <- prune_taxa(taxa_sums(Fecal_data_micro.ps) > 0, Fecal_data_micro.ps)


# 1) Get metadata (make sure these are factors with the intended baselines)
df <- data.frame(sample_data(Fecal_data_micro.ps))
df$trt             <- factor(df$trt)               # outcome (e.g., 0/1 or two-level factor)
df$timepoint  <- factor(df$timepoint)
df$Timepoint <- factor(df$Timepoint)
df$treatment_order <- factor(df$treatment_order)
df$participant_id  <- factor(df$participant_id)


# (If trt is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$trt) - 1  # baseline = first level in df$trt

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(Fecal_data_micro.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel

# 3) Covariates (no intercept column for X)
X <- model.matrix(~ timepoint, data = df)[, -1, drop = FALSE]

# 4) GLMM-MiRKAT (random intercepts via id)
set.seed(1)
fit <- MiRKAT::GLMMMiRKAT(
  y      = y,
  X      = X,
  Ks     = K,                        # or list(K1, K2, ...) if multiple kernels
  id     = df$participant_id,        # random intercept for subject
  model  = "binomial",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE,                    # random intercept only
  nperm  = 5000                      # (not n.perm)
)

fit$p_values      # kernel-specific p-value(s)
# 0.3007399 6Dec2025, only timepoint
#0.1655669 03Dec2025


### Treatment order ####

# (If trt is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$treatment_order) - 1  # baseline = first level in df$trt

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(Fecal_data_micro.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel

# 3) Covariates (no intercept column for X)
X <- model.matrix(~ Timepoint * treatment_order, data = df)[, -1, drop = FALSE]

# 4) GLMM-MiRKAT (random intercepts via id)
set.seed(1)
trt_order_fit <- fit <- MiRKAT::GLMMMiRKAT(
  y      = y,
  X      = X,
  Ks     = K,                        # or list(K1, K2, ...) if multiple kernels
  id     = df$participant_id,        # random intercept for subject
  model  = "binomial",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE,                    # random intercept only
  nperm  = 5000                      # (not n.perm)
)

trt_order_fit$p_values      # kernel-specific p-value(s)
#0.1655669 03Dec2025


## Beef sample microbiome ####
Beef_AMR_data.ps <- subset_samples(data_micro.ps, sample_type == "Meat Rinsate" )

# (optional) drop zero-sum taxa left behind
Beef_AMR_data.ps <- prune_taxa(taxa_sums(Beef_AMR_data.ps) > 0, Beef_AMR_data.ps)

#Beef_AMR_data.ps <- subset_samples(Beef_AMR_data.ps,!(Participant_ID %in% c("113", "115", "117", "118", "135")))
Beef_AMR_data.ps <- prune_taxa(taxa_sums(Beef_AMR_data.ps) > 0, Beef_AMR_data.ps)

# 1) Get metadata (make sure these are factors with the intended baselines)
df <- data.frame(sample_data(Beef_AMR_data.ps))
df$trt             <- factor(df$trt)               # outcome (e.g., 0/1 or two-level factor)
df$timepoint  <- factor(df$timepoint)
df$Timepoint <- factor(df$Timepoint)
df$treatment_order <- factor(df$treatment_order)
df$participant_id  <- factor(df$participant_id)



# (If trt_id is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$trt) - 1  # baseline = first level in df$trt_id

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(Beef_AMR_data.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel



# 4) GLMM-MiRKAT (random intercepts via id)
set.seed(1)
fit <- MiRKAT::MiRKAT(
  y        = y,
  X        = NULL,
  Ks       = K,
  out_type = "D"   # "D" for dichotomous/binary
)

fit$p_values      # kernel-specific p-value(s)  
# 0.790242 with "Start" sampling times, Treatment as outcome 

#0.4621076 Nov 22 MM
#fit$omnibus_p     # omnibus p-value if you passed multiple Ks

