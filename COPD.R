########################
# Pre-processing starts
########################
matrixes <- list.files(path = "data_required/", full.names = T, pattern = ".gz")


# Read and merge all files by "ID_REF"
merged_data <- lapply(matrixes, function(file) {
  # Read file
  df <- read.table(file, header = TRUE, sep = "\t", comment.char = "!", fill = TRUE, quote = "", stringsAsFactors = FALSE)
  # Keep only relevant columns
  df
})

# View merged result
# head(merged_data)


# Rename the first column in each data frame to "ID" for uniform merging
merged_data <- lapply(merged_data, function(df) {
  colnames(df)[1] <- "ID"
  return(df)
})

# Now merge all data frames by "ID" using Reduce
final_merged_data <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), merged_data)

# View merged data
# head(final_merged_data)
dim(final_merged_data)


# Clean up column names
clean_names <- function(names_vec) {
  names_vec <- gsub("^X\\.", "", names_vec)         # Remove leading "X."
  names_vec <- gsub("\\.\\.x$|\\.\\.y$|\\.$", "", names_vec)  # Remove trailing "..x", "..y", or "."
  return(names_vec)
}

# Apply to your data
colnames(final_merged_data) <- clean_names(colnames(final_merged_data))
names(final_merged_data)

# Check for duplicates
# dup_cols <- colnames(final_merged_data)[duplicated(colnames(final_merged_data))]
# print(dup_cols)

# Remove duplicated columns, keeping the first occurrence
# final_merged_data <- final_merged_data[, !duplicated(colnames(final_merged_data))]

# Remove double quotes from 'ID' column
final_merged_data$ID <- gsub('"', '', final_merged_data$ID)

# soft file for annotation
soft <- read.delim("data_required/GPL570-55999.txt")
soft <- as.data.frame(soft[,c(1,11)])


# Merge gene symbols from soft into final_merged_data based on ID
final_annotated_data <- merge(soft, final_merged_data, by = "ID")


# Remove rows with duplicated Gene.Symbols, keeping the first occurrence
final_annotated_data_unique <- final_annotated_data[!duplicated(final_annotated_data$Gene.Symbol) & 
                                                      !is.na(final_annotated_data$Gene.Symbol), ]


# Assign 'ID' column as rownames
rownames(final_annotated_data_unique) <- final_annotated_data_unique$Gene.Symbol

# Drop the 'ID' column
final_annotated_data_unique$ID <- NULL
final_annotated_data_unique$Gene.Symbol <- NULL

dim(final_annotated_data_unique)
# write.csv(final_annotated_data_unique,"output/final_annotated_data_unique.csv")



# Load phenotype data
# Load required libraries
library(GEOquery)
library(Biobase)
library(dplyr)
library(stringr)

# Function to clean and set row names
clean_and_set_rownames <- function(df) {
  df$ID <- str_replace_all(df$ID, '"', '')  # Remove quotation marks
  rownames(df) <- df$ID
  df$ID <- NULL  # Optionally remove the ID column after setting row names
  return(df)
}

# List of GEO datasets
gse_ids <- c("GSE11784", "GSE13896", "GSE11906", "GSE37768", "GSE130928")

# Function to annotate samples based on dataset
annotate_samples <- function(gse_id) {
  gse <- getGEO(gse_id, GSEMatrix = TRUE)
  if (length(gse) > 1) idx <- grep("GPL", attr(gse, "names")) else idx <- 1
  gse_data <- gse[[idx]]
  pheno_data <- pData(gse_data)
  
  # Annotate phenotype
  pheno_data$phenotype <- case_when(
    grepl("non-smoker", pheno_data$title, ignore.case = TRUE) |
      grepl("tissue_NS", pheno_data$title, ignore.case = TRUE) |
      grepl("\\bNS\\b", pheno_data$title, ignore.case = TRUE) ~ "Healthy Non-Smoker",
    
    grepl("COPD", pheno_data$title, ignore.case = TRUE) ~ "Smoker with COPD",
    
    grepl("smoker", pheno_data$title, ignore.case = TRUE) & 
      !grepl("COPD", pheno_data$title, ignore.case = TRUE) ~ "Smoker",
    
    grepl("SYMs", pheno_data$title, ignore.case = TRUE) |
      grepl("\\bS\\b", pheno_data$title, ignore.case = TRUE) |
      grepl("tissue_S", pheno_data$title, ignore.case = TRUE) |
      grepl("LODs", pheno_data$title, ignore.case = TRUE) ~ "Smoker",
    
    TRUE ~ "Unknown"
  )
  
  # Add GSE_ID column
  pheno_data$GSE_ID <- gse_id
  
  return(pheno_data[, c("GSE_ID", "title", "phenotype")])
}

# Annotate samples for all datasets
annotated_samples_list <- lapply(gse_ids, annotate_samples)

# Combine all annotations into a single data frame
all_annotated_samples <- bind_rows(annotated_samples_list)


# View result
# View(all_annotated_samples)

# List of Unknown samples
Unknown <- data.frame(all_annotated_samples[all_annotated_samples$phenotype == "Unknown", ])
# write.csv(Unknown, "output/Unknown samples removed.csv", row.names = F)


# Remove unknown phenotypes
# all_annotated_samples_clean <- all_annotated_samples %>%
#   filter(phenotype != "Unknown")

all_annotated_samples_withIDs <- data.frame("GSM_IDs" = row.names(all_annotated_samples), all_annotated_samples)


# Remove '...numbers' from GSM IDs
all_annotated_samples_withIDs$GSM_IDs <- gsub("\\.\\.\\..*", "", all_annotated_samples_withIDs$GSM_IDs)
# write.csv(all_annotated_samples_withIDs, "output/all_annotated_samples_withIDs.csv", row.names = T)

# View the combined annotations
# View(all_annotated_samples_withIDs)


# Remove rows with duplicates
all_annotated_samples_withIDs_unique <- all_annotated_samples_withIDs[!duplicated(all_annotated_samples_withIDs$GSM_IDs) & 
                                                         !is.na(all_annotated_samples_withIDs$GSM_IDs), ]

# Sort the dataframe by phenotype
all_annotated_samples_withIDs_unique_sorted <- all_annotated_samples_withIDs_unique[order(all_annotated_samples_withIDs_unique$phenotype), ]
row.names(all_annotated_samples_withIDs_unique_sorted) <- seq_len(nrow(all_annotated_samples_withIDs_unique_sorted))

# Final phenotype dataframe
phenotype <- all_annotated_samples_withIDs_unique_sorted[, c(1, 2, 4)]

# remove Unknown samples
phenotype <- phenotype[phenotype$phenotype != "Unknown", ]

# write.csv(phenotype, "output/phenotype_final.csv", row.names = F) # nolint





# 1. Keep only columns in final_annotated_data_unique that are in phenotype$GSM_IDs
common_samples <- intersect(phenotype$GSM_IDs, colnames(final_annotated_data_unique))
filtered_expression <- final_annotated_data_unique[, common_samples]


# 2. Sort phenotype by phenotype column (e.g., alphabetical order) (already done before but just to be sure)
phenotype_sorted <- phenotype[order(phenotype$phenotype), ]
table(phenotype_sorted$phenotype)

# 3. Reorder the columns of expression data to match sorted phenotype
filtered_expression <- filtered_expression[, phenotype_sorted$GSM_IDs]
# write.csv(filtered_expression,"filtered_expression.csv")

# Subset GSM IDs by phenotype group
hn_ids <- phenotype_sorted$GSM_IDs[phenotype_sorted$phenotype == "Healthy Non-Smoker"]
smoker_ids <- phenotype_sorted$GSM_IDs[phenotype_sorted$phenotype == "Smoker"]
copd_ids <- phenotype_sorted$GSM_IDs[phenotype_sorted$phenotype == "Smoker with COPD"]

# Create subset expression matrices
expr_hn <- filtered_expression[, hn_ids]
expr_smoker <- filtered_expression[, smoker_ids]
expr_copd <- filtered_expression[, copd_ids]
# Healthy vs COPD expression matrix
expr_hnVScopd <- filtered_expression[,c(hn_ids,copd_ids)]
# write.csv(expr_hnVScopd,"output/expr_hnVScopd.csv")
# Smoker vs COPD expression matrix
expr_smokerVScopd <- filtered_expression[,c(smoker_ids,copd_ids)]
# write.csv(expr_smokerVScopd,"output/expr_smokerVScopd.csv")



# Subset phenotype information for the selected samples
pheno_hnVScopd <- phenotype_sorted[phenotype_sorted$GSM_IDs %in% colnames(expr_hnVScopd), ]

# Ensure phenotype vector matches the column order
pheno_hnVScopd <- pheno_hnVScopd[match(colnames(expr_hnVScopd), pheno_hnVScopd$GSM_IDs), ]
#####################
# Pre-processing end
#####################


#############################
# Batch effect removal starts
#############################
# Create a sample group mapping
sample_groups <- data.frame(
  SampleID = colnames(expr_hnVScopd),
  Group = c(rep("HN", 161), rep("COPD", 92))  # Adjust based on your actual group sizes
)
# Extract gene symbols
gene_symbols <- sapply(strsplit(rownames(expr_hnVScopd), " /// "), `[`, 1)

# Identify valid (non-NA and non-duplicated) entries
valid <- !is.na(gene_symbols) & !duplicated(gene_symbols)

# Subset the expression matrix and assign new rownames
expr_hnVScopd <- expr_hnVScopd[valid, ]
rownames(expr_hnVScopd) <- gene_symbols[valid]


# log normalization
expr_hnVScopd_log <- log2(expr_hnVScopd + 1)

# Load the limma package for normalization
library(limma)
expr_log2_norm <- normalizeBetweenArrays(expr_hnVScopd_log, method = "quantile")



# Batch effect removal
library(sva)

# Batch correction
batch_corrected <- ComBat(dat = as.matrix(expr_log2_norm), batch = pheno_hnVScopd$phenotype, par.prior = TRUE)
write.csv(batch_corrected, "output/batch_corrected.csv")
head(batch_corrected[1:5, 1:5])
batch_corrected <- read.csv("output/batch_corrected.csv", row.names = 1)


#############################
# Batch effect removal done
#############################
