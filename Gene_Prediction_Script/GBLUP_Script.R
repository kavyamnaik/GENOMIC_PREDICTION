###############################################################
# Genomic Prediction Pipeline
# GBLUP Model - To Predict GEBVs using SNP Markers
# Author: Dr. Kavyashree N M
###############################################################

# --- 1. Load Libraries ---
# ---- Load libraries ----
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(tibble)
library(emmeans)

# ---- Define file paths ----
data_raw_path <- "Gene_Prediction_Data/Raw/Geno_Pheno_Raw_Data.xlsx"
pheno_raw <- read_excel(data_raw_path, sheet ="Phenotype_Raw") #Read phenotypic raw data
geno_raw  <- read_excel(data_raw_path, sheet = "Genotype_Raw")

str(pheno_raw) #To check phenotypic raw data structure
str(geno_raw) #To check genotypic raw data structure

# ==========================================================
# Phenotypic Data Quality Assessment
# ==========================================================
dim(pheno_raw) #To view dimensions (no. of rows & columns)
colnames(pheno_raw) #To view Column names

#Convert design variables to factors (One trait at a time)
pheno_yield <- pheno_raw %>%
  filter(!is.na(Yield)) %>%           # For yield trait
  mutate(
    Genotype    = as.factor(Genotype),
    Environment = as.factor(Environment),
    Rep = as.factor(Rep)
  ) %>%
  select(Genotype, Environment, Rep, Yield)

# Sanity checks (Missing yield- allowed but missing genotype/environment-not allowed)
stopifnot(!anyNA(pheno_yield$Genotype))
stopifnot(!anyNA(pheno_yield$Environment))

cat("Phenotype summary\n")
print(table(pheno_yield$Genotype))
print(table(pheno_yield$Environment))

# Check the structure of data/experimental design
cat("\nNumber of observations per Environment:\n")
table(pheno_yield$Environment)

cat("\nNumber of observations per Rep:\n")
table(pheno_yield$Rep)

cat("\nReplications within each Environment:\n")
with(pheno_yield, table(Environment, Rep))

# =============================================================================
# Environment code harmonization (handling inconsistent location labels)
# =============================================================================

#SAME LOCATION, DIFFERENT CODING HAS BEEN GIVEN. E1 ≈ BLR_Heat, E2 ≈ HYD_Normal
pheno_yield <- pheno_yield %>%
  mutate(Environment = recode(Environment,
                              "E1" = "BLR_Heat",
                              "E2" = "HYD_Normal"))

pheno_yield <- pheno_yield %>%
  mutate(
    Location = case_when(
      Environment %in% c("E1", "BLR_Heat") ~ "Bengaluru",
      Environment %in% c("E2", "HYD_Normal") ~ "Hyderabad",
      TRUE ~ NA_character_
    )) #henceforth we will use locations as environment to avoid confusion, as we renamed E1 and E2

table(pheno_yield$Environment, pheno_yield$Location)
table(pheno_yield$Location, pheno_yield$Rep) #to check replication pattern in each location


# =============================================================================
# Compute BLUEs (Yield ~ Genotype + Location + (1 | Location:Rep))
# =============================================================================
blue_yield <- lmer(
  Yield ~ Genotype + Location + (1 | Location:Rep),
  data = pheno_yield
)

summary(blue_yield)

# Extract genotype BLUEs (This will give adjusted means corrected for Location + Rep.)
blue_means <- emmeans(blue_yield, "Genotype")
blue_means

pheno_for_gblup <- as.data.frame(blue_means) %>%
  select(Genotype, emmean) %>%
  rename(BLUE_Yield = emmean) #Covert into dataframe

pheno_for_gblup #final phenotypic input for gblup

# =============================================================================
# Preparation of genotypic data matrix
# =============================================================================
# Convert genotype data to matrix
geno_mat <- geno_raw %>%
  column_to_rownames("GenotypeID") %>%
  as.matrix()

str(geno_mat)
dim(geno_mat)

#Check genotype names in both phenotypic & genotypic data
pheno_for_gblup$Genotype # Genotypes in phenotype (BLUEs)
rownames(geno_mat) # Genotypes in genotype matrix

#Find common genotypes
common_genotypes <- intersect(
  pheno_for_gblup$Genotype,
  rownames(geno_mat)
)

common_genotypes #To check common genotypes in both phenotypic & genotypic input set
length(common_genotypes) #To check if both phenotypic and genotypic inputs are aligning perfectly, to avoid mismatches

# =============================================================================
# Align genotype matrix to phenotyped genotypes
# =============================================================================
#Subset genotype matrix to use only those genotypes with BLUES and SNP Marker data
geno_use <- geno_mat[common_genotypes, ] 

dim(geno_use) #tells about the number of training individuals and predicting markers in our data
rownames(geno_use)

#Handle missing SNP data by Mean imputation per marker
geno_imp <- apply(geno_use, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

anyNA(geno_imp) #to verify replacement of NA values in marker data
dim(geno_imp)

#Center the marker matrix
Z <- scale(geno_imp, center = TRUE, scale = FALSE) # Center markers by subtracting marker means

round(colMeans(Z), 5)[1:10] #verify centering

# =============================================================================
# Genomic Relationship Matrix (GRM)- VanRaden method
# =============================================================================
# Calculate allele frequencies
p <- colMeans(geno_imp) / 2

# Denominator for VanRaden GRM
denominator <- 2 * sum(p * (1 - p)) #p is the allele frequency at each SNP and 2p(1−p) is the expected genetic variance per marker

# Genomic Relationship Matrix (GRM)
G <- (Z %*% t(Z)) / denominator #G is the estimated genomic similarity between genotypes

dim(G)
round(diag(G), 3)

all(rownames(G) == pheno_for_gblup$Genotype) #To check whether phenotype order align with GRM order

# =============================================================================
# Genomic Prediction-GBLUP
# =============================================================================

#Prepare phenotype vector for GBLUP
y <- pheno_for_gblup$BLUE_Yield #Create phenotype vector y
names(y) <- pheno_for_gblup$Genotype
y

#Fit the GBLUP model using GRM
library(sommer)
pheno_gblup <- data.frame(
  Genotype = names(y),
  Yield = as.numeric(y)
) #Prepare data frame for sommer

str(pheno_gblup)

write.csv(
  pheno_for_gblup,
  "Gene_Prediction_Results/Phenotype_BLUEs_Yield.csv",
  row.names = FALSE
)


gblup_model <- mmer(
  Yield ~ 1,
  random = ~ vs(Genotype, Gu = G),
  rcov   = ~ units,
  data   = pheno_gblup
) #Fit GBLUP Model

# =============================================================================
# Extract GEBVs (Genomic Estimated Breeding Values)
# =============================================================================

# Extract genomic BLUPs (GEBVs)
gblup_effects <- gblup_model$U[["u:Genotype"]]$Yield

#Create clean GEBV table
gblup_results <- data.frame(
  Genotype = names(gblup_effects),
  GEBV = as.numeric(gblup_effects)
)

gblup_results

names(gblup_model$U)

# =============================================================================
# Final ranking of genotypes based on GEBVs
# =============================================================================
gblup_results <- gblup_results %>%
  arrange(desc(GEBV))

gblup_results

write.csv(
  gblup_results,
  "Gene_Prediction_Results/GEBV_Ranking_Yield.csv",
  row.names = FALSE
)

