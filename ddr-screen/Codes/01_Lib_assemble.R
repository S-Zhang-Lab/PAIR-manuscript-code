library(stringi)
library(tidyr)
library(dplyr)
library(DNABarcodes)
library(stringr)

# generate random 14 bit DNA barcode with the desired length and Hamming distance 5 for guaranteed 2 bit error correction
initial_barcodes_14 <- create.dnabarcodes(14, dist = 5, cores=24)
initial_barcodes <- initial_barcodes_14

# Function to check for continuous identical nucleotides
has_continuous_nucleotides <- function(barcode, n = 3) {
  for (i in 1:(nchar(barcode) - n + 1)) {
    if (substring(barcode, i, i + n - 1) == paste0(rep(substring(barcode, i, i), n), collapse = "")) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# Function to check for balanced nucleotides with a tolerance level
is_balanced <- function(barcode, tolerance = 2) {
  counts <- table(strsplit(barcode, NULL)[[1]])
  # Ensure that counts of A, T, C, G are approximately equal within the given tolerance
  return(all(abs(counts - mean(counts)) <= tolerance))
}

# Filter barcodes for continuous identical nucleotides
filtered_barcodes <- initial_barcodes[!sapply(initial_barcodes, has_continuous_nucleotides)]

# Further filter for balanced nucleotides
balanced_barcodes <- filtered_barcodes[sapply(filtered_barcodes, is_balanced)]

# Check if we have enough barcodes after filtering
if (length(balanced_barcodes) < 11025) {
  stop("Not enough barcodes available after filtering. Consider relaxing some constraints.")
}

# Randomly pick 11025 barcodes
set.seed(123)  # Set seed for reproducibility
selected_barcodes <- sample(balanced_barcodes, 11025, replace = FALSE)

# Analyze the selected barcodes using analyse.barcodes
barcode_analysis <- analyse.barcodes(selected_barcodes)

# Output the analysis results and selected barcodes
print(barcode_analysis)

# view selected barcodes
selected_barcodes[1:50]


# ==  Step 2 Assemble crRNA design table ==
# load crRNA design table
crRNA_table <- read.csv("Misc/crRNA_table.csv", sep="\t")
df <- crRNA_table
df$num <- seq(1:5)
head(df)

# Convert sequences to uppercase
df$CRISPRa <- toupper(df$CRISPRa)
df$CasRx <- toupper(df$CasRx)

# Convert to long format
df_long <- pivot_longer(df, cols = c(CRISPRa, CasRx), names_to = "Type", values_to = "Sequence")
df_long$name <- paste(df_long$ID, df_long$Type, df_long$num, sep = "_")
head(df_long)

# Define oligo elements
e1 <- "AGGGCCTATTTCCCATGATTcgtctcacaccg"
e2 <- "N20"
e3 <- "gttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgcttCAAGTAAACCCCTACCAACTGGTCGGGGTTTGAAAC"
e4 <- "n23"
e5 <- "TTTTTTT"
e6 <- "BC14"
e7 <- "ctacagagacgcacttgtacttcagcggtc"

# Extract sequences based on Type
CRISPRa_sequences <- df_long %>% filter(Type == "CRISPRa")
CasRx_sequences <- df_long %>% filter(Type == "CasRx")

# Create all possible combinations of indices
combinations <- expand.grid(CRISPRa_index = 1:nrow(CRISPRa_sequences), 
                            CasRx_index = 1:nrow(CasRx_sequences))

# Add sequence and name columns based on indices
combinations <- combinations %>%
  mutate(CRISPRa = CRISPRa_sequences$Sequence[CRISPRa_index],
         CasRx = CasRx_sequences$Sequence[CasRx_index],
         CRISPRa_name = CRISPRa_sequences$name[CRISPRa_index],
         CasRx_name = CasRx_sequences$name[CasRx_index])
head(combinations)

combinations_all <- combinations %>%
  mutate(Final_oligos = paste0(e1, CRISPRa, e3, CasRx, e5, selected_barcodes, e7))
combinations_all[11020:11025,5:7]

write.csv(combinations_all, file = "./Misc/Final_oligos_allcombinations.csv", row.names = FALSE)

# Extract rows where both CRISPRa and CasRx columns start with "NT"
filtered_combinations <- combinations_all[grepl("^NT", combinations_all$CRISPRa_name) & grepl("^NT", combinations_all$CasRx_name), ]
head(filtered_combinations)

# Extract the gene names from the "name" column for comparison
combinations_all <- combinations_all %>%
  mutate(CRISPRa_ID = sub("_CRISPRa.*", "", CRISPRa_name),
         CasRx_ID = sub("_CasRx.*", "", CasRx_name))
head(combinations_all)

# Filter out combinations where the sequences are from the same gene name (ID)
combinations_all <- combinations_all %>% filter(CRISPRa_ID != CasRx_ID)
head(combinations_all)

# Create the final sequences
combinations_all$CRISPRa_ID <- NULL
combinations_all$CasRx_ID <- NULL

Final <- rbind(combinations_all, filtered_combinations)

head(Final)
tail(Final)

# Extract barcode sequences
Final$BC14 <- str_extract(Final$Final_oligos, "(?<=TTTTTTT).{14}(?=ctacagagacgca)")

# View results
print(Final)

# Save the final table to a CSV file with comma separation
write.csv(Final, file = "./Misc/Final_oligos.csv", row.names = FALSE)

# add oligo names
Final_oligos <- read.csv("./Misc/Final_oligos.csv")

Final_oligos$ID <- paste("PAIR", seq(1:10525), sep = "_")
write.csv(Final_oligos, file = "./Misc/Final_oligos_withID.csv", row.names = FALSE)
