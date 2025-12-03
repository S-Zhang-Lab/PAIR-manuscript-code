rm(list = ls())
set.seed(7)
library(data.table) 
library(Seurat)
library(dplyr)
library(fgsea)
library(qs)

# Import raw data and metadata----
scRNA_seq <- qread("./data/seu_prep.qs")
meta_data <- scRNA_seq@meta.data

# Convert loge to log2
scRNA_seq <- log2(exp(scRNA_seq[["RNA"]]$data)) %>% as.matrix()

all(colnames(scRNA_seq) == rownames(meta_data)) # TRUE

treatments <- c("CTRL", "RNP")

treatment <- "RNP"

chosen_idx <- meta_data$treatment == treatment
meta_data <- meta_data[chosen_idx, ]
scRNA_seq <- scRNA_seq[, chosen_idx]
all(colnames(scRNA_seq) == rownames(meta_data)) # TRUE

NT_CRISPRa_NT_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NT_CRISPRa_NT_CasRx"]

NBN_CRISPRa_NT_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NBN_CRISPRa_NT_CasRx"]
MRE11_CRISPRa_NT_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "MRE11_CRISPRa_NT_CasRx"]

NT_CRISPRa_53BP1_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NT_CRISPRa_53BP1_CasRx"]
NT_CRISPRa_DNPKcs_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NT_CRISPRa_DNPKcs_CasRx"]
NT_CRISPRa_KU70_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NT_CRISPRa_KU70_CasRx"]
NT_CRISPRa_PQ_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NT_CRISPRa_PQ_CasRx"]

NBN_CRISPRa_53BP1_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NBN_CRISPRa_53BP1_CasRx"]
NBN_CRISPRa_DNPKcs_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NBN_CRISPRa_DNPKcs_CasRx"]
NBN_CRISPRa_KU70_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NBN_CRISPRa_KU70_CasRx"]
NBN_CRISPRa_PQ_CasRx_exp <- scRNA_seq[, meta_data$assigned_tag == "NBN_CRISPRa_PQ_CasRx"]

# Calculate log2 fold changes
NBN_NT_over_NT_NT <- sort(rowMeans(NBN_CRISPRa_NT_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                          decreasing = TRUE)
MRE11_NT_over_NT_NT <- sort(rowMeans(MRE11_CRISPRa_NT_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                          decreasing = TRUE)
NT_53BP1_over_NT_NT <- sort(rowMeans(NT_CRISPRa_53BP1_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                             decreasing = TRUE)
NT_DNPKcs_over_NT_NT <- sort(rowMeans(NT_CRISPRa_DNPKcs_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                            decreasing = TRUE)
NT_KU70_over_NT_NT <- sort(rowMeans(NT_CRISPRa_KU70_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                             decreasing = TRUE)
NT_PQ_over_NT_NT <- sort(rowMeans(NT_CRISPRa_PQ_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                           decreasing = TRUE)

NBN_DNPKcs_over_NT_NT <- sort(rowMeans(NBN_CRISPRa_DNPKcs_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                             decreasing = TRUE)
NBN_KU70_over_NT_NT <- sort(rowMeans(NBN_CRISPRa_KU70_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                           decreasing = TRUE)
NBN_PQ_over_NT_NT <- sort(rowMeans(NBN_CRISPRa_PQ_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                         decreasing = TRUE)
NBN_53BP1_over_NT_NT <- sort(rowMeans(NBN_CRISPRa_53BP1_CasRx_exp) - rowMeans(NT_CRISPRa_NT_CasRx_exp), 
                          decreasing = TRUE)

pwy <- "c2"
signature_all <- readRDS(paste0("./data/pathways/", pwy, ".rds"))

base_excision_repair <- c(
  "WP_BASE_EXCISION_REPAIR",
  "REACTOME_POLB_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR",
  "REACTOME_PCNA_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR",
  "REACTOME_BASE_EXCISION_REPAIR",
  "REACTOME_BASE_EXCISION_REPAIR_AP_SITE_FORMATION",
  "REACTOME_DISEASES_OF_BASE_EXCISION_REPAIR",
  "KEGG_BASE_EXCISION_REPAIR"
)

dsb_repair <- c(
  "WP_DNA_IRDOUBLE_STRAND_BREAKS_DSBS_AND_CELLULAR_RESPONSE_VIA_ATM",
  "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
  "REACTOME_SENSING_OF_DNA_DOUBLE_STRAND_BREAKS",
  "REACTOME_DNA_DOUBLE_STRAND_BREAK_RESPONSE",
  "REACTOME_PROCESSING_OF_DNA_DOUBLE_STRAND_BREAK_ENDS"
)

recombinational_repair <- c(
  "WP_HOMOLOGOUS_RECOMBINATION",
  "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR",
  "REACTOME_MEIOTIC_RECOMBINATION",
  "REACTOME_INHIBITION_OF_DNA_RECOMBINATION_AT_TELOMERE",
  "KEGG_HOMOLOGOUS_RECOMBINATION"
)

cell_cycle_checkpoint <- c(
  "REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_MACHINERY_AT_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS",
  "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
  "REACTOME_G2_M_DNA_REPLICATION_CHECKPOINT",
  "REACTOME_G2_M_CHECKPOINTS",
  "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS",
  "REACTOME_MITOTIC_SPINDLE_CHECKPOINT",
  "REACTOME_CELL_CYCLE_CHECKPOINTS",
  "REACTOME_THE_ROLE_OF_GTSE1_IN_G2_M_PROGRESSION_AFTER_G2_CHECKPOINT"
)

signature <- signature_all[cell_cycle_checkpoint]

# Run fGSEA----
fgseaRes_NBN_NT_over_NT_NT <- fgsea(signature, NBN_NT_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NBN_NT_over_NT_NT <- fgseaRes_NBN_NT_over_NT_NT[fgseaRes_NBN_NT_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NBN_NT_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NBN_NT_over_NT_NT_", pwy, ".csv"))

fgseaRes_MRE11_NT_over_NT_NT <- fgsea(signature, MRE11_NT_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_MRE11_NT_over_NT_NT <- fgseaRes_MRE11_NT_over_NT_NT[fgseaRes_MRE11_NT_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_MRE11_NT_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/MRE11_NT_over_NT_NT_", pwy, ".csv"))

fgseaRes_NT_53BP1_over_NT_NT <- fgsea(signature, NT_53BP1_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NT_53BP1_over_NT_NT <- fgseaRes_NT_53BP1_over_NT_NT[fgseaRes_NT_53BP1_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NT_53BP1_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NT_53BP1_over_NT_NT_", pwy, ".csv"))

fgseaRes_NT_DNPKcs_over_NT_NT <- fgsea(signature, NT_DNPKcs_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NT_DNPKcs_over_NT_NT <- fgseaRes_NT_DNPKcs_over_NT_NT[fgseaRes_NT_DNPKcs_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NT_DNPKcs_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NT_DNPKcs_over_NT_NT_", pwy, ".csv"))

fgseaRes_NT_KU70_over_NT_NT <- fgsea(signature, NT_KU70_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NT_KU70_over_NT_NT <- fgseaRes_NT_KU70_over_NT_NT[fgseaRes_NT_KU70_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NT_KU70_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NT_KU70_over_NT_NT_", pwy, ".csv"))

fgseaRes_NT_PQ_over_NT_NT <- fgsea(signature, NT_PQ_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NT_PQ_over_NT_NT <- fgseaRes_NT_PQ_over_NT_NT[fgseaRes_NT_PQ_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NT_PQ_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NT_PQ_over_NT_NT_", pwy, ".csv"))

fgseaRes_NBN_DNPKcs_over_NT_NT <- fgsea(signature, NBN_DNPKcs_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NBN_DNPKcs_over_NT_NT <- fgseaRes_NBN_DNPKcs_over_NT_NT[fgseaRes_NBN_DNPKcs_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NBN_DNPKcs_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NBN_DNPKcs_over_NT_NT_", pwy, ".csv"))

fgseaRes_NBN_KU70_over_NT_NT <- fgsea(signature, NBN_KU70_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NBN_KU70_over_NT_NT <- fgseaRes_NBN_KU70_over_NT_NT[fgseaRes_NBN_KU70_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NBN_KU70_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NBN_KU70_over_NT_NT_", pwy, ".csv"))

fgseaRes_NBN_PQ_over_NT_NT <- fgsea(signature, NBN_PQ_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NBN_PQ_over_NT_NT <- fgseaRes_NBN_PQ_over_NT_NT[fgseaRes_NBN_PQ_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NBN_PQ_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NBN_PQ_over_NT_NT_", pwy, ".csv"))

fgseaRes_NBN_53BP1_over_NT_NT <- fgsea(signature, NBN_53BP1_over_NT_NT, minSize=15, maxSize=500)
# fgseaRes_NBN_53BP1_over_NT_NT <- fgseaRes_NBN_53BP1_over_NT_NT[fgseaRes_NBN_53BP1_over_NT_NT$padj < 0.05, ]
fwrite(fgseaRes_NBN_53BP1_over_NT_NT, paste0("./res/2025_1120/pwy_results", "_", treatment, "/NBN_53BP1_over_NT_NT_", pwy, ".csv"))


