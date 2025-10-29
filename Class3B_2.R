#######################################################################
#### 0. Install and Load Required Packages ####
#######################################################################

# Bioconductor provides R packages for analyzing omics data (genomics, transcriptomics, proteomics etc).

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))
BiocManager::install("GEOquery", force = TRUE)
BiocManager::install("affy", force = TRUE)
BiocManager::install("arrayQualityMetrics", force = TRUE)

# Install CRAN packages for data manipulation
install.packages("dplyr")

# Load Required Libraries
library(GEOquery)             # Download GEO datasets (series matrix, raw CEL files)
library(affy)                 # Pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(dplyr)                # Data manipulation

# -------------------------------------
#### Download Series Matrix Files ####
# -------------------------------------
gse_data <- getGEO("GSE55267", GSEMatrix = TRUE)

# Extract expression data matrix (genes/probes × samples)
expression_data <- exprs(gse_data$GSE55267_series_matrix.txt.gz)


# Extract feature (probe annotation) data
feature_data <-  fData(gse_data$GSE55267_series_matrix.txt.gz)


# Extract phenotype (sample metadata) data
phenotype_data <-  pData(gse_data$GSE55267_series_matrix.txt.gz)

# Check missing values in sample annotation
sum(is.na(phenotype_data$source_name_ch1)) 

# --------------------------------------
#### Download Raw Data (CEL files) ####
# --------------------------------------
# Fetch GEO supplementry files
getGEOSuppFiles("GSE55267", baseDir = "Raw_Data1", makeDirectory = TRUE)

# Untar CEL files if compressed as .tar
untar("Raw_Data1/GSE55267_RAW.tar", exdir = "Raw_Data1/CEL_Files")

# Read CEL files into R as an AffyBatch object
raw_data <- ReadAffy(celfile.path = "Raw_Data1/CEL_Files")

raw_data

# ---------------------------------------------------
#### Quality Control (QC) Before Pre-processing ####
# ---------------------------------------------------
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)

# -------------------------------------------------------
#### RMA (Robust Multi-array Average) Normalization ####
# -------------------------------------------------------
normalized_data <- rma(raw_data)

# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)


processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data)   

# ---------------------------------------------------------------------------
#### Filter Low-Variance Transcripts (“soft” intensity based filtering) ####
# ---------------------------------------------------------------------------
row_median <- rowMedians(as.matrix(processed_data))

# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

# Set a threshold 
threshold <- 3.5 
abline(v = threshold, col = "purple", lwd = 2) 

# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 

colnames(filtered_data) <- rownames(phenotype_data)

processed_data <- filtered_data 

# -----------------------------------
#### Phenotype Data Preparation ####
# -----------------------------------
class(phenotype_data$source_name_ch1) 
colnames(phenotype_data)

unique(phenotype_data$source_name_ch1)
unique(phenotype_data$characteristics_ch1.3)


# Define experimental groups (normal vs cancer)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("Pre-treated fresh frozen FL biopsy", "Purified from normal reactive tonsil"),
                 label = c("Cancer", "Normal"))

class(groups)
levels(groups)
# ---------------------------------------------------------------------------
#### Save workspace ####
# ---------------------------------------------------------------------------

save.image(file = "GSE30727.RData")

# ---------------------------------------------------------------------------
#### Loading workspace ####
# ---------------------------------------------------------------------------
load("GSE30727.RData")
