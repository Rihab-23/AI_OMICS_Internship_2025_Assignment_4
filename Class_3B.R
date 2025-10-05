# ===============================================================
#                    Task_BENSALEK Rihab
# ---------------------------------------------------------
# ---------------------------------------------------------
#     Module II: Introduction to Genomics Data Analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
#         Assignment 4:   Microarray Data Analysis
# ===============================================================

# -------------------------------
# 0. Install and Load Packages
# -------------------------------
# Install Bioconductor packages
BiocManager::install(c("lumi", "arrayQualityMetrics", "illuminaio"), ask = FALSE)
BiocManager::install("illuminaio", force = TRUE)
BiocManager::install("lumi", force = TRUE)
BiocManager::install("arrayQualityMetrics", force = TRUE)

# Install CRAN packages for data manipulation
install.packages("dplyr")

# Load Required Libraries
library(Biobase)              # for ExpressionSet / AnnotatedDataFrame
library(lumi)                 # Illumina expression arrays
library(arrayQualityMetrics)  # QC reports
library(illuminaio)           # readIDAT function
library(dplyr)                # data manipulation

# -------------------------------
# 1. Set Paths
# -------------------------------
idat_dir <- "C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_II/E-MTAB-8148_(1)"
metadata_file1 <- "C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_II/MAGE-TAB_Files/E-MTAB-8148.sdrf.txt"

# -------------------------------
# 2. List IDAT Files
# -------------------------------
idat_files <- list.files(idat_dir, pattern = "[iI][dD][aA][tT]$", full.names = TRUE)
cat("Total IDAT files found:", length(idat_files), "\n")

# -------------------------------
# 3. Read IDAT Files in Batches
# -------------------------------
batch_size <- 50
intensity_matrix <- NULL
probe_names <- NULL
successful_files <- character(0)

for(i in seq(1, length(idat_files), by = batch_size)){
  batch_files <- idat_files[i:min(i+batch_size-1, length(idat_files))]
  cat("Reading batch", i, "to", min(i+batch_size-1, length(idat_files)), "\n")
  batch_list <- lapply(batch_files, function(f){
    ok <- try(readIDAT(f), silent = TRUE)
    if(inherits(ok, "try-error")) {
      warning("Failed to read: ", f)
      return(NULL)
    }
    ok
  })
  # remove NULLs (failed reads) and keep names of successful files
  valid_idx <- !sapply(batch_list, is.null)
  if(!any(valid_idx)) next
  batch_list <- batch_list[valid_idx]
  successful_files <- c(successful_files, batch_files[valid_idx])
  
  if(is.null(probe_names)) probe_names <- batch_list[[1]]$Manifest$Name
  batch_matrix <- sapply(batch_list, function(x){
    # adapt depending on readIDAT structure; try these options
    if(!is.null(x$Quants$Mean)) return(x$Quants$Mean)
    if(!is.null(x$Intensity)) return(x$Intensity)
    stop("Cannot find intensity vector in readIDAT output")
  })
  # ensure matrix orientation: probes × samples
  if(is.vector(batch_matrix)) batch_matrix <- matrix(batch_matrix, ncol = 1)
  if(is.null(intensity_matrix)) intensity_matrix <- batch_matrix else intensity_matrix <- cbind(intensity_matrix, batch_matrix)
  
  rm(batch_list, batch_matrix); gc()
}

# sanity checks
if(is.null(intensity_matrix)) stop("No IDATs were successfully read.")
rownames(intensity_matrix) <- probe_names
colnames(intensity_matrix) <- gsub("\\.idat$","", basename(successful_files))

cat("Final intensity matrix dims (probes x samples):", dim(intensity_matrix), "\n")
cat("Successful files read:", length(successful_files), "\n")

# 2) Align phenotype metadata to the actual samples we have
# Try to find column in phenotype_data that holds the raw filename (common names: Array.Data.File, Source.Name)
# We'll try Array.Data.File first, otherwise look for a column that contains sample names
md <- phenotype_data

# Detect candidate column that contains the file names (loosely)
cands <- c("Array.Data.File","Array.Data.File_name","Array.Data.File","Array.Data.File.URI","Source.Name")
found <- intersect(cands, colnames(md))
if(length(found) == 0){
  # fallback: look for any column which contains the first sample name (partial match)
  maybe <- sapply(md, function(col) any(grepl(gsub("\\.idat$","", basename(successful_files)[1]), col, ignore.case = TRUE)))
  if(any(maybe)){
    found <- names(maybe)[which(maybe)[1]]
    message("Auto-detected phenotype column: ", found)
  } else {
    message("No obvious phenotype column found. Showing head of metadata columns:")
    print(colnames(md))
    stop("Please identify which metadata column matches the sample names from IDAT filenames.")
  }
}

file_col <- found[1]
cat("Using phenotype column:", file_col, "\n")

# create matching vector: remove .idat suffixes from phenotype column values for matching
phen_names <- gsub("\\.idat$","", as.character(md[[file_col]]))
sample_names_actual <- colnames(intensity_matrix)

ord <- match(sample_names_actual, phen_names)
if(any(is.na(ord))){
  missing_samples <- sample_names_actual[is.na(ord)]
  stop("The following samples are present in intensity matrix but not found in phenotype metadata column '", file_col, "':\n", paste(missing_samples, collapse = ", "))
}

rownames(intensity_matrix) <- probe_names
colnames(intensity_matrix) <- gsub("\\.idat$","", basename(successful_files))

cat("Intensity matrix dims (probes x samples):", dim(intensity_matrix), "\n")
cat("Successful files read:", length(successful_files), "\n")

# -------------------------------
# 4. Load Metadata
# -------------------------------
phenotype_data <- read.delim(metadata_file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Optional: check column names and candidate columns for disease info
colnames(phenotype_data)
candidate_cols <- grep("Characteristics|Factor Value|Source Type", colnames(phenotype_data), value = TRUE)

# -------------------------------
# 5. Align Phenotype Metadata with Samples
# -------------------------------
# Match sample names between metadata and intensity matrix
phen_names <- gsub("\\.idat$","", phenotype_data$Array.Data.File)
ord <- match(colnames(intensity_matrix), phen_names)
md2 <- phenotype_data[ord, , drop = FALSE]
rownames(md2) <- colnames(intensity_matrix)

# -------------------------------
# 6. Build ExpressionSet
# -------------------------------
pheno_adf <- AnnotatedDataFrame(md2)
eset <- ExpressionSet(assayData = as.matrix(intensity_matrix), phenoData = pheno_adf)

# Quick checks to make sure data loaded successfully & correctly 
cat("ExpressionSet dimensions (probes x samples):", dim(exprs(eset)), "\n") # dim(exprs(eset)) → prints the number of probes × number of samples, tells how many features and samples are there
cat("pData rows:", nrow(pData(eset)), 
    "should equal ncol(exprs(eset)):", 
    ncol(exprs(eset)), "\n") # nrow(pData(eset)) vs ncol(exprs(eset)) → verifies that the number of rows in the phenotype data matches the number of samples in the expression matrix.

# -------------------------------
# 7. QC Before Normalization
# -------------------------------
arrayQualityMetrics(expressionset = eset,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)
# OUTPUT: HTML report with flagged arrays (outliers)

# -------------------------------
# 8. Normalization
# -------------------------------
exprs(eset) <- normalizeBetweenArrays(exprs(eset), method = "quantile")
eset_norm <- eset

# QC after normalization
arrayQualityMetrics(expressionset = eset_norm,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)
# OUTPUT: HTML report post-normalization

# -------------------------------
# 9. Filter Low-Intensity Probes
# -------------------------------
processed_data <- exprs(eset_norm)
row_median <- apply(processed_data, 1, median)

hist(row_median, breaks = 100, freq = FALSE, main = "Median Probe Intensity")
threshold <- 120  # chosen based on histogram inspection
abline(v = threshold, col = "blue", lwd = 2)

filtered_data <- processed_data[row_median > threshold, ]
cat("Probes remaining after filtering:", nrow(filtered_data), "\n")
summary(row_median) # check the distribution of probe medians

#Probes (features) present before filtering" -> nrow(exprs(eset_norm))
#Transcripts remained after filtering low-intensity probes -> nrow(filtered_data)

# -------------------------------
# 10. Define Experimental Groups
# -------------------------------
# Identify which metadata column contains disease info
groups <- factor(phenotype_data$Characteristics.disease.,
                 levels = c("Normal mucosa adjacent", "colorectal cancer"),
                 labels = c("normal", "colorectal cancer"))

# Assign groups to ExpressionSet
pData(eset_norm)$group <- groups
table(pData(eset_norm)$group)
# OUTPUT: counts of normal vs cancer samples

# -------------------------------
# 11. Summary for Assignment
# -------------------------------
nrow(exprs(eset_norm))          # Probes before filtering
nrow(filtered_data)             # Probes after filtering
length(colnames(eset_norm))     # Total samples
table(pData(eset_norm)$group)   # Sample counts by group

# -------------------------------
# 12.Saving the ExpressionSet after normalization
# -------------------------------
save(eset_norm, file = "eset_norm.RData")

# -------------------------------
# 13.Saving the normalized expression matrix & filtered data
# -------------------------------
write.csv(exprs(eset_norm), file = "normalized_expression_matrix.csv")
write.csv(filtered_data, file = "filtered_expression_matrix.csv")

# -------------------------------
# 14.Saving the entire R workspace
# -------------------------------
save.image(file = "BENSALEK_Rihab_Class3B_Assignment4.RData")

# -------------------------------
# 15.Loading Workspace
# -------------------------------
load("BENSALEK_Rihab_Class3B_Assignment4.RData")