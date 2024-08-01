library(stringi)
library(tidyr)
library(dplyr)

# generate random 10bit DNA barcode
BC <- read.csv("Misc/BC.txt", sep="\t")

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

# Define e1, e3, e5 and e7
e1 <- "GAGGGCCTATTTCCCATGATTcgtctcacacc"
e3 <- "gttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgcttCAAGTAAACCCCTACCAACTGGTCGGGGTTTGAAAC"
e5 <- "TTTTTTT"
e6 <- BC
e7 <- "ctacagagacgcacttgtacttcagcggtca"

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

# Extract the gene names from the "name" column for comparison
combinations <- combinations %>%
  mutate(CRISPRa_ID = sub("_CRISPRa.*", "", CRISPRa_name),
         CasRx_ID = sub("_CasRx.*", "", CasRx_name))
head(combinations)

# Filter out combinations where the sequences are from the same gene name (ID)
combinations <- combinations %>% filter(CRISPRa_ID != CasRx_ID)
head(combinations)

# Create the final sequences
combinations <- combinations %>%
  mutate(Final_oligos = paste0(e1, CRISPRa, e3, CasRx, e5))

# Create the final table
Final <- data.frame(seq = seq(1:nrow(combinations)), combinations)

# Rename the sequence column to "seq_ID"
colnames(Final)[1] <- "seq_ID"

# Display the final data frame
print(head(Final))
dim(Final)

# Save the final table to a CSV file with comma separation
write.csv(Final, file = "./Misc/Final_oligos.csv", row.names = FALSE)

