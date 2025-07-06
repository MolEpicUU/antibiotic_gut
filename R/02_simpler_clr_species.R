# Project Antibiotic use and gut microbiota 

library(data.table)

# This script will perform a CLR transformation on the gut microbiota relative abundances

  # CLR transformation function:
  clr_fun <- function(spmat){
  
  zeros <- lapply(1:ncol(spmat), function(i) which(spmat[,i]==0))
  
  min_non_zero <- min(spmat[spmat>0])
  
  tempmat <- spmat + min_non_zero
  
  mgsclr <- as.matrix(as.data.frame(compositions::clr(tempmat)))

  # Replacing the inicial zeros with the minimal value for every species.
  for(i in 1:ncol(mgsclr)){
    if(length(zeros[[i]])>0){
    min <- min(mgsclr[,i][ -zeros[[i]] ])
    mgsclr[,i][ zeros[[i]]  ] <- min
    }
  }
  
  data.frame(Subject = rownames(mgsclr), mgsclr)
  }
  

  # SIMPLER --------------------------------------------------------------------
  message("SIMPLER")
  simpler_mgs <- fread('/proj/simp2023007/data/microbiome/processed/simpler_metagenomics_mgs_relative_abundances_v1.0.tsv', data.table = F)
  simpler <- fread("/proj/nobackup/simp2023007/users/baldanzi/work/simpler_working_dataset.csv")
  
  simpler_mgs <- simpler_mgs[simpler_mgs$SIMPKEY %in% simpler$SIMPKEY, ]
  
  rownames(simpler_mgs) <- simpler_mgs$SIMPKEY
  simpler_mgs <- as.matrix(simpler_mgs[,-1])
  
  allzeros <- colnames(simpler_mgs)[which(colSums(simpler_mgs)==0)]
  simpler_mgs <- simpler_mgs[, -which(colnames(simpler_mgs) %in% allzeros)]

  simpler_clr <- clr_fun(simpler_mgs)
  setnames(simpler_clr, "Subject", "SIMPKEY")
  fwrite(simpler_clr, file = '/proj/nobackup/simp2023007/users/baldanzi/work/simpler_clr_species.tsv')
  
  
  message("End")
  
  