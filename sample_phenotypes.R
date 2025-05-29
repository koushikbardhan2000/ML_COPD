# Load required libraries
library(GEOquery)
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
View(all_annotated_samples)

# List of Unknown samples
Unknown <- data.frame(all_annotated_samples[all_annotated_samples$phenotype == "Unknown", ])
write.csv(Unknown, "output/Unknown samples removed.csv", row.names = F)


# Remove unknown phenotypes
# all_annotated_samples_clean <- all_annotated_samples %>%
#   filter(phenotype != "Unknown")

all_annotated_samples_withIDs <- data.frame("GSM_IDs" = row.names(all_annotated_samples), all_annotated_samples)


# Remove '...numbers' from GSM IDs
all_annotated_samples_withIDs$GSM_IDs <- gsub("\\.\\.\\..*", "", all_annotated_samples_withIDs$GSM_IDs)
write.csv(all_annotated_samples_withIDs, "output/all_annotated_samples_withIDs.csv", row.names = T)

# View the combined annotations
View(all_annotated_samples_withIDs)


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

write.csv(phenotype, "output/phenotype_final.csv", row.names = F) # nolint