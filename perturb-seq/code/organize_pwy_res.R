rm(list = ls())
set.seed(7)
library(data.table) 
library(Seurat)
library(dplyr)
library(fgsea)
library(qs)
library(VennDiagram)


venn_regions_2 <- function(A, B) {
  # extract the names of the variables passed to the function
  arg_names <- as.list(match.call())[-1]
  nameA <- deparse(arg_names[[1]])
  nameB <- deparse(arg_names[[2]])
  
  # compute regions
  regions <- list(
    A_only = setdiff(A, B),
    B_only = setdiff(B, A),
    A_B    = intersect(A, B)
  )
  
  # rename based on input variable names
  names(regions) <- c(
    paste0(nameA, "_only"),
    paste0(nameB, "_only"),
    paste0(nameA, "_", nameB)
  )
  
  return(regions)
}


venn_regions <- function(A, B, C) {
  # extract the names of the variables passed to the function
  arg_names <- as.list(match.call())[-1]
  nameA <- deparse(arg_names[[1]])
  nameB <- deparse(arg_names[[2]])
  nameC <- deparse(arg_names[[3]])
  
  # construct region names dynamically
  names_list <- list(
    setdiff(A, union(B, C)),
    setdiff(B, union(A, C)),
    setdiff(C, union(A, B)),
    
    setdiff(intersect(A, B), C),
    setdiff(intersect(A, C), B),
    setdiff(intersect(B, C), A),
    
    Reduce(intersect, list(A, B, C))
  )
  
  names(names_list) <- c(
    paste0(nameA, "_only"),
    paste0(nameB, "_only"),
    paste0(nameC, "_only"),
    
    paste0(nameA, "_", nameB),
    paste0(nameA, "_", nameC),
    paste0(nameB, "_", nameC),
    
    paste0(nameA, "_", nameB, "_", nameC)
  )
  
  return(names_list)
}


treatments <- c("CTRL", "RNP")
treatment <- "RNP"


# ==== NBN_DNPKcs ====
NBN_NT <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NBN_NT_over_NT_NT_c2.csv"))
NBN_NT <- NBN_NT$pathway

NT_DNPKcs <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NT_DNPKcs_over_NT_NT_c2.csv"))
NT_DNPKcs <- NT_DNPKcs$pathway
NBN_DNPKcs <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NBN_DNPKcs_over_NT_NT_c2.csv"))
NBN_DNPKcs <- NBN_DNPKcs$pathway

venn_list <- list(
  NBN_NT = NBN_NT,
  NT_DNPKcs = NT_DNPKcs,
  NBN_DNPKcs = NBN_DNPKcs
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,               # draw to R instead of a file
  fill = c("#A0CBE8", "#F28E2B", "#76B7B2"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  lwd = 2,
  cat.pos = c(-20, 20, 0),
  cat.dist = c(0.05, 0.05, 0.05)
)

png(paste0("./res/2025_1113/organize_pwy_", treatment, "/venn_NBN_DNPKcs.png"),
    width = 7, height = 7, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(venn.plot)
dev.off()

regions <- venn_regions(NBN_NT, NT_DNPKcs, NBN_DNPKcs)
# Determine max length among all regions
max_len <- max(lengths(regions))

# Pad each region with NA to equalize lengths
regions_padded <- lapply(regions, function(x) {
  c(x, rep(NA, max_len - length(x)))
})

# Convert to data frame
regions_df <- as.data.frame(regions_padded, stringsAsFactors = FALSE)
write.csv(regions_df,
          file = paste0("./res/2025_1113/organize_pwy_", treatment, "/venn_regions_NBN_DNPKcs.csv"),
          row.names = FALSE,
          na = "")


# ==== NBN_KU70 ====
NT_KU70 <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NT_KU70_over_NT_NT_c2.csv"))
NT_KU70 <- NT_KU70$pathway
NBN_KU70 <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NBN_KU70_over_NT_NT_c2.csv"))
NBN_KU70 <- NBN_KU70$pathway


venn_list <- list(
  NBN_NT = NBN_NT,
  NT_KU70 = NT_KU70,
  NBN_KU70 = NBN_KU70
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,               # draw to R instead of a file
  fill = c("#A0CBE8", "#F28E2B", "#76B7B2"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  lwd = 2,
  cat.pos = c(-20, 20, 0),
  cat.dist = c(0.05, 0.05, 0.05)
)

png(paste0("./res/2025_1113/organize_pwy_", treatment, "/venn_NBN_KU70.png"),
    width = 7, height = 7, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(venn.plot)
dev.off()

regions <- venn_regions(NBN_NT, NT_KU70, NBN_KU70)
# Determine max length among all regions
max_len <- max(lengths(regions))

# Pad each region with NA to equalize lengths
regions_padded <- lapply(regions, function(x) {
  c(x, rep(NA, max_len - length(x)))
})

# Convert to data frame
regions_df <- as.data.frame(regions_padded, stringsAsFactors = FALSE)
write.csv(regions_df,
          file = paste0("./res/2025_1113/organize_pwy_", treatment, "/venn_regions_NBN_KU70.csv"),
          row.names = FALSE,
          na = "")


# ==== NBN_PQ ====
NT_PQ <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NT_PQ_over_NT_NT_c2.csv"))
NT_PQ <- NT_PQ$pathway
NBN_PQ <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NBN_PQ_over_NT_NT_c2.csv"))
NBN_PQ <- NBN_PQ$pathway

venn_list <- list(
  NBN_NT = NBN_NT,
  NT_PQ = NT_PQ,
  NBN_PQ = NBN_PQ
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,               # draw to R instead of a file
  fill = c("#A0CBE8", "#F28E2B", "#76B7B2"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  lwd = 2,
  cat.pos = c(-20, 20, 0),
  cat.dist = c(0.05, 0.05, 0.05)
)

png(paste0("./res/2025_1113/organize_pwy_", treatment, "/venn_NBN_PQ.png"),
    width = 7, height = 7, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(venn.plot)
dev.off()

regions <- venn_regions(NBN_NT, NT_PQ, NBN_PQ)
# Determine max length among all regions
max_len <- max(lengths(regions))

# Pad each region with NA to equalize lengths
regions_padded <- lapply(regions, function(x) {
  c(x, rep(NA, max_len - length(x)))
})

# Convert to data frame
regions_df <- as.data.frame(regions_padded, stringsAsFactors = FALSE)
write.csv(regions_df,
          file = paste0("./res/2025_1113/organize_pwy_", treatment, "/venn_regions_NBN_PQ.csv"),
          row.names = FALSE,
          na = "")


# ==== NBN_53BP1 ====
NT_53BP1 <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NT_53BP1_over_NT_NT_c2.csv"))
NT_53BP1 <- NT_53BP1$pathway
NBN_53BP1 <- read.csv(paste0("./res/2025_1113/pwy_results_", treatment, "/NBN_53BP1_over_NT_NT_c2.csv"))
NBN_53BP1 <- NBN_53BP1$pathway

venn_list <- list(
  NBN_NT = NBN_NT,
  NT_53BP1 = NT_53BP1,
  NBN_53BP1 = NBN_53BP1
)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,               # draw to R instead of a file
  fill = c("#A0CBE8", "#F28E2B", "#76B7B2"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  lwd = 2,
  cat.pos = c(-20, 20, 0),
  cat.dist = c(0.05, 0.05, 0.05)
)

png(paste0("./res/2025_1113/organize_pwy_", treatment, "/venn_NBN_53BP1.png"),
    width = 7, height = 7, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(venn.plot)
dev.off()


regions <- venn_regions(NBN_NT, NT_53BP1, NBN_53BP1)
# Determine max length among all regions
max_len <- max(lengths(regions))

# Pad each region with NA to equalize lengths
regions_padded <- lapply(regions, function(x) {
  c(x, rep(NA, max_len - length(x)))
})

# Convert to data frame
regions_df <- as.data.frame(regions_padded, stringsAsFactors = FALSE)
write.csv(regions_df,
          file = paste0("./res/2025_1113/organize_pwy_", treatment, "/venn_regions_NBN_53BP1.csv"),
          row.names = FALSE,
          na = "")
