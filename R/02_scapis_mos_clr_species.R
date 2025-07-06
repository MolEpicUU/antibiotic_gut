# Project Antibiotic use and gut microbiota 

library(data.table)

# This script will perform a CLR transformation on the gut microbiota relative abundances

  # CLR transformation function:
  clr_fun <- function(spmat){
  
  zeros <- lapply(1:ncol(spmat), function(i) which(spmat[,i]==0))
  
  min_non_zero <- min(spmat[spmat>0])
  
  tempmat <- spmat + min_non_zero
  
  message("\nCLR Transformation started")
  
  mgsclr <- as.matrix(as.data.frame(compositions::clr(tempmat)))

  # Replacing the inicial zeros with the minimal value for every species.
  for(i in 1:ncol(mgsclr)){
    
    min <- min(mgsclr[,i][ -zeros[[i]] ])
    mgsclr[,i][ zeros[[i]]  ] <- min
  }
    data.frame(Subject = rownames(mgsclr), mgsclr)
  
  }
  
  # SCAPIS ---------------------------------------------------------------------
  message("SCAPIS")
  scapis_mgs <- fread('data/scapis/raw/scapis_metagenomics_mgs_relative_abundances_v1.0.tsv', data.table=F)

  rownames(scapis_mgs) <- scapis_mgs$scapis_id
  scapis_mgs$scapis_id <- NULL
  scapis_mgs <- as.matrix(scapis_mgs)

  scapis_clr <- clr_fun(scapis_mgs)
  fwrite(scapis_clr, file = 'atb_gut/work/scapis_clr_species.tsv')


  # MOS ------------------------------------------------------------------------
  message("MOS")
  mos <- fread("atb_gut/work/mos_working_dataset.tsv", na.strings = c("NA",NA,""))
  microb <- fread('MOS/Metagenomics/mos_metagenomics_mgs_relative_abundances_v1.0.tsv', data.table = F)
  
  microb <- microb[microb$lopnrMOS %in% mos$lopnrMOS, ]

  rownames(microb) <- microb$lopnrMOS
  
  allzeros <- names(microb[,-1])[colSums(microb[,-1])==0]
  microb <- microb[, -which(colnames(microb) %in% allzeros)]
  
  
  
  microb_clr <- clr_fun(as.matrix(microb[,-1]))
  
  message("\nCLR Transformation finished")
  setnames(microb_clr, "Subject", "lopnrMOS")
  fwrite(microb_clr, file = 'atb_gut/work/mos_clr_species.tsv')
  
  message("\nEnd")
  
  