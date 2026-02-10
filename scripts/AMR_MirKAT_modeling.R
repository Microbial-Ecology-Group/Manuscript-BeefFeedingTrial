# Packages
library(phyloseq)
library(vegan)
library(MiRKAT)        # for D2K to convert distance -> kernel
#remotes::install_github("hk1785/GLMM-MiRKAT")
library(GLMMMiRKAT)


Fecal_AMR_data.ps <- subset_samples(AMR_data.ps, Sample_Type == "Fecal" & !(col_time %in% c("Washout")) # ,"Start"
)

# (optional) drop zero-sum taxa left behind
Fecal_AMR_data.ps <- prune_taxa(taxa_sums(Fecal_AMR_data.ps) > 0, Fecal_AMR_data.ps)

#Fecal_AMR_data.ps <- subset_samples(Fecal_AMR_data.ps,!(Participant_ID %in% c("113", "115", "117", "118", "135")))
Fecal_AMR_data.ps <- prune_taxa(taxa_sums(Fecal_AMR_data.ps) > 0, Fecal_AMR_data.ps)

# 1) Get metadata (make sure these are factors with the intended baselines)
df <- data.frame(sample_data(Fecal_AMR_data.ps))
df$trt_id             <- factor(df$trt_id)               # outcome (e.g., 0/1 or two-level factor)
df$col_time  <- factor(df$col_time)
df$treatment_order <- factor(df$treatment_order)
df$Participant_ID  <- factor(df$Participant_ID)


# (If trt_id is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$trt_id) - 1  # baseline = first level in df$trt_id

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(Fecal_AMR_data.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel

#Fix rowname mismatch 
#rownames(K) <- rownames(as.matrix(D))
#colnames(K) <- rownames(as.matrix(D))
#df <- df[rownames(K), , drop = FALSE]
#y  <- y[match(rownames(K), names(y))]

# 3) Covariates (no intercept column for X)
X <- model.matrix(~ col_time * treatment_order, data = df)[, -1, drop = FALSE]


# 4) GLMM-MiRKAT (random intercepts via id)
set.seed(1)
fit <- MiRKAT::GLMMMiRKAT(
  y      = y,
  X      = X,
  Ks     = K,                        # or list(K1, K2, ...) if multiple kernels
  id     = df$Participant_ID,        # random intercept for subject
  model  = "binomial",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE #,                    # random intercept only
  #nperm  = 5000                      # (not n.perm)
)

fit$p_values      # kernel-specific p-value(s)  
# 0.790242 with "Start" sampling times, Treatment as outcome 

#0.4621076 Nov 22 MM
#fit$omnibus_p     # omnibus p-value if you passed multiple Ks


### Just treatment as outcome, controlling for Time ####

# 1) Get metadata (make sure these are factors with the intended baselines)
df <- data.frame(sample_data(Fecal_AMR_data.ps))
df$trt_id             <- factor(df$trt_id)               # outcome (e.g., 0/1 or two-level factor)

df$Timepoint  <- factor(df$Timepoint)
df$col_time  <- factor(df$col_time)
df$treatment_order <- factor(df$treatment_order)
df$Participant_ID  <- factor(df$Participant_ID)


# (If trt_id is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$trt_id) - 1  # baseline = first level in df$trt_id

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(Fecal_AMR_data.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel

#Fix rowname mismatch 
#rownames(K) <- rownames(as.matrix(D))
#colnames(K) <- rownames(as.matrix(D))
#df <- df[rownames(K), , drop = FALSE]
#y  <- y[match(rownames(K), names(y))]

# 3) Covariates (no intercept column for X)
X <- model.matrix(~ col_time , data = df)[, -1, drop = FALSE]


# 4) GLMM-MiRKAT (random intercepts via id)
set.seed(1)
fit <- MiRKAT::GLMMMiRKAT(
  y      = y,
  X      = X,
  Ks     = K,                        # or list(K1, K2, ...) if multiple kernels
  id     = df$Participant_ID,        # random intercept for subject
  model  = "binomial",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE #,                    # random intercept only
  #nperm  = 5000                      # (not n.perm)
)

fit$p_values  
# 0.41


## Doing by treatment order ####

# (If trt_id is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$treatment_order) - 1  # baseline = first level in df$trt_id

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(Fecal_AMR_data.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel


# 3) Covariates (no intercept column for X)
X <- model.matrix(~ Timepoint * trt_id , data = df)[, -1, drop = FALSE]


# 4) GLMM-MiRKAT (random intercepts via id)
set.seed(1)
fit <- MiRKAT::GLMMMiRKAT(
  y      = y,
  X      = X,
  Ks     = K,                        # or list(K1, K2, ...) if multiple kernels
  id     = df$Participant_ID,        # random intercept for subject
  model  = "binomial",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE #,                    # random intercept only
  #nperm  = 5000                      # (not n.perm)
)

fit$p_values      # kernel-specific p-value(s)  
# 0.1263747

# 0.295141 p value by treatment order, including "start" time
#0.4621076 Nov 22 MM
#fit$omnibus_p     # omnibus p-value if you passed multiple Ks




#######################Run Mirkat without Random Effects #########################
#Added MM 22 Nov 2025

No_Rand_effects.ps <- subset_samples(AMR_data.ps, Sample_Type == "Fecal" & !(col_time %in% c("Washout","Start")) # ,"Start"
)

# (optional) drop zero-sum taxa left behind
No_Rand_effects.ps <- prune_taxa(taxa_sums(No_Rand_effects.ps) > 0, No_Rand_effects.ps)


# 1) Get metadata (make sure these are factors with the intended baselines)
df <- data.frame(sample_data(No_Rand_effects.ps))
df$trt_id             <- factor(df$trt_id)               # outcome (e.g., 0/1 or two-level factor)
df$col_time  <- factor(df$col_time)
df$treatment_order <- factor(df$treatment_order)
df$Participant_ID  <- factor(df$Participant_ID)


# (If trt_id is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$trt_id) - 1  # baseline = first level in df$trt_id

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(No_Rand_effects.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel


# 3) Covariates (no intercept column for X)
X <- model.matrix(~ col_time + treatment_order, data = df)[, -1, drop = FALSE]


#Run Mirkat without Participant ID but WITH ALL participants
set.seed(1)
fit <- MiRKAT(
  y  = y,     # binary outcome (0/1)
  X  = X,     # covariates
  Ks = K      # kernel (or list of kernels)
)

fit$p_values
# 0.9611829


### Beef samples ###
Beef_AMR_data.ps <- subset_samples(AMR_data.ps, Sample_Type == "Beef Sample" & !(col_time %in% c("Washout")))
Beef_AMR_data.ps <- subset_samples(Beef_AMR_data.ps, specific_sample_type != "Sirloin")



# (optional) drop zero-sum taxa left behind
Beef_AMR_data.ps <- prune_taxa(taxa_sums(Beef_AMR_data.ps) > 0, Beef_AMR_data.ps)


# 1) Get metadata (make sure these are factors with the intended baselines)
df <- data.frame(sample_data(Beef_AMR_data.ps))
df$trt_id             <- factor(df$trt_id)               # outcome (e.g., 0/1 or two-level factor)
df$col_time  <- factor(df$col_time)
df$treatment_order <- factor(df$treatment_order)
df$Participant_ID  <- factor(df$Participant_ID)


# (If trt_id is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$trt_id) - 1  # baseline = first level in df$trt_id

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


### Just treatment as outcome, controlling for Time ####

# 1) Get metadata (make sure these are factors with the intended baselines)
df <- data.frame(sample_data(Beef_AMR_data.ps))
df$trt_id             <- factor(df$trt_id)               # outcome (e.g., 0/1 or two-level factor)
df$col_time  <- factor(df$col_time)
df$treatment_order <- factor(df$treatment_order)
df$Participant_ID  <- factor(df$Participant_ID)


# (If trt_id is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$trt_id) - 1  # baseline = first level in df$trt_id

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(Beef_AMR_data.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel

#Fix rowname mismatch 
#rownames(K) <- rownames(as.matrix(D))
#colnames(K) <- rownames(as.matrix(D))
#df <- df[rownames(K), , drop = FALSE]
#y  <- y[match(rownames(K), names(y))]

# 3) Covariates (no intercept column for X)
X <- model.matrix(~ col_time , data = df)[, -1, drop = FALSE]


# 4) GLMM-MiRKAT (random intercepts via id)
set.seed(1)
fit <- MiRKAT::GLMMMiRKAT(
  y      = y,
  X      = X,
  Ks     = K,                        # or list(K1, K2, ...) if multiple kernels
  id     = df$Participant_ID,        # random intercept for subject
  model  = "binomial",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE #,                    # random intercept only
  #nperm  = 5000                      # (not n.perm)
)

fit$p_values  



## Doing by treatment order ####

# (If trt_id is binary, set y to 0/1; if it has >2 levels, see note below)
y <- as.numeric(df$treatment_order) - 1  # baseline = first level in df$trt_id

# 2) Build a distance and convert to a kernel
D <- phyloseq::distance(Beef_AMR_data.ps, method = "bray")
labs <- rownames(as.matrix(D))
df   <- df[labs, , drop = FALSE]             # align to distance order
K    <- MiRKAT::D2K(as.matrix(D))            # distance -> kernel


# 3) Covariates (no intercept column for X)
X <- model.matrix(~ col_time * trt_id , data = df)[, -1, drop = FALSE]


# 4) GLMM-MiRKAT (random intercepts via id)
set.seed(1)
fit <- MiRKAT::GLMMMiRKAT(
  y      = y,
  X      = X,
  Ks     = K,                        # or list(K1, K2, ...) if multiple kernels
  id     = df$Participant_ID,        # random intercept for subject
  model  = "binomial",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE #,                    # random intercept only
  #nperm  = 5000                      # (not n.perm)
)

fit$p_values      # kernel-specific p-value(s)  

# 0.295141 p value by treatment order, including "start" time
#0.4621076 Nov 22 MM
#fit$omnibus_p     # omnibus p-value if you passed multiple Ks


