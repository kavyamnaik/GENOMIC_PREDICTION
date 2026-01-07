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

##################ERROR FIXING##############
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
# Compute BLUEs (Genotypes-> fixed, environment->random )
# =============================================================================
blue_yield <- lmer(
  Yield ~ Genotype + (1 | Environment),
  data = pheno_yield
)
summary(blue_yield)

# Extract genotype BLUPs (these will be used as 'phenotypes' for GBLUP)
blup_gen <- ranef(mod_yield)$Genotype
blup_gen <- data.frame(
  Genotype = rownames(blup_gen),
  GEBV_pheno = blup_gen[,"(Intercept)"],
  row.names = NULL
)