rm(list = ls())
set.seed(7)

library(dplyr)
library(msigdb)
library(GSEABase)
msigdb.mm <- getMsigdb(org = 'hs', id = 'SYM', version = '7.4')
msigdb.mm <- appendKEGG(msigdb.mm) # KEGG is included for mouse
listCollections(msigdb.mm)

pwy_collections <- c("h", "c2")

for (i in 1:length(pwy_collections)) {
  pwy <- pwy_collections[i]
  
  chosen_pathways <- subsetCollection(msigdb.mm, collection =  c(pwy))
  
  listSubCollections(chosen_pathways)
  
  # 4. Create pathway list for hallmark gene sets----
  signature <- c()
  # names(signature) <- names(chosen_pathways)
  for (i in 1:length(chosen_pathways)) {
    geneID_name_i <- chosen_pathways[[i]]@setName
    signature[[geneID_name_i]] <- chosen_pathways[[i]]@geneIds
  }
  save_path <- paste0("./data/pathways/", pwy, ".rds")
  saveRDS(signature, save_path)
}
