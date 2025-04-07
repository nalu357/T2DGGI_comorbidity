library(tidyverse)
library(mrclust)

################################################
#------------------- Functions -----------------
################################################
get.data <- function(df, data.dir){
  data <- TwoSampleMR::format_data(type = data.dir,
                                   dat = as.data.frame(df),
                                   snp_col = "rs_id",
                                   beta_col = "beta",
                                   se_col = "standard_error",
                                   effect_allele_col = "effect_allele",
                                   other_allele_col = "other_allele",
                                   eaf_col = "effect_allele_frequency",
                                   pval_col = "p_value",
                                   chr_col = "chr",
                                   pos_col = "pos",
                                   ncase_col = "ncases",
                                   samplesize_col = "n")
  return(data)
}

wrap.MRClust <- function(exposure.file, outcome.file, exposure.name, outcome.name, ivs, savepath, ancest){
  res.dt <- data.table::data.table()
  if(length(ivs)==0) return(res.dt)
  
  #----------- Get exposure -----------
  exposure <- get.data(data.table::fread(exposure.file), data.dir="exposure")
  exposure <- exposure[exposure$SNP %in% ivs,]
  if(nrow(exposure)==0) return(res.dt)
  exposure$exposure <- exposure.name
  print("Exposure Read")
  print(colnames(exposure))
  
  
  #----------- Get outcome -----------
  outcome <- get.data(data.table::fread(outcome.file), data.dir="outcome")
  outcome$outcome <- outcome.name
  print("Outcome Read")
  print(colnames(outcome))

  #----------- Harmonize data -----------
  data <- tryCatch(
    TwoSampleMR::harmonise_data(exposure_dat=exposure, outcome_dat=outcome),
    error=function(e) e
  )
  
  if(inherits(data, "error")) {
    print("Harmonization failed.")
    return(res.dt) 
  }
  
  if(nrow(data)==0) return(res.dt)
    
  if(nrow(subset(data, mr_keep==TRUE))==0) {
    print("No SNP left after data harmonization, not able to run MR Clust :(")
    return(res.dt)
  }

  data$id.exposure <- data$exposure
  data$id.outcome <- data$outcome
  data <- subset(data, mr_keep==TRUE)

  #----------- Run MRClust  -----------
  snp_names = data$SNP
  bx = data$beta.exposure
  bxse = data$se.exposure
  by = data$beta.outcome
  byse = data$se.outcome
  ratio_est = by/bx
  ratio_est_se = byse/abs(bx)
  
  res_em = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx, by = by, bxse = bxse, byse = byse, obs_names = snp_names)
  res80 = mrclust::pr_clust(dta = res_em$results$best, prob = 0.8, min_obs =  4)

  keep80 = which(snp_names %in% res80$observation)
  bx80 = bx[keep80]
  bxse80 = bxse[keep80]
  by80 = by[keep80]
  byse80 = byse[keep80]
  snp_names80 = snp_names[keep80]	
  
  Plot.FullProb = res_em$plots$two_stage +
    ggplot2::xlim(0, max(abs(bx) + 2*bxse)) +
    ggplot2::xlab(paste0("Genetic Association with ",exposure.name)) +
    ggplot2::ylab(paste0("Genetic Association with ",outcome.name)) +
    ggplot2::ggtitle("")
	
  Plot.Prob80 = two_stage_plot(res = res80, bx = bx80, by = by80, bxse = bxse80, byse = byse80, obs_names = snp_names80) + 
    ggplot2::xlim(0, max(abs(bx80) + 2*bxse80)) +
	ggplot2::xlab(paste0("Genetic association with ",exposure.name)) + 
	ggplot2::ylab(paste0("Genetic association with ",outcome.name)) + 
	ggplot2::ggtitle("")	

  #----------- Save Datasets  -----------
  write_csv(res_em$results$all, paste0(savepath,"AllSNP_",exposure.name,"_",ancest,"_Exp.csv"))
  write_csv(res_em$results$best, paste0(savepath,"BestSNP_",exposure.name,"_",ancest,"_Exp.csv"))
  write_csv(res80, paste0(savepath,"PR80_",exposure.name,"_",ancest,"_Exp.csv"))
  ggsave(filename = paste0(savepath,"Full_",exposure.name,"_",ancest,"_Exp.png"), plot = Plot.FullProb, width = 7, height = 7)
  ggsave(filename = paste0(savepath,"Prob80_",exposure.name,"_",ancest,"_Exp.png"), plot = Plot.Prob80, width = 7, height = 7)
  
}

RunMRClust <- function(datasets.dt){
for(pop in unique(datasets.dt$ANCESTRY)){
    for(t in datasets.dt[ANCESTRY==pop, TRAIT]){
      trait <- datasets.dt[ANCESTRY==pop & TRAIT==t]
	  #Read T2D IVs
	  t2d.all.ivs <- data.table::fread(paste0(trait$PATH, "t2d_", pop, "_ivs.txt"))
	  t2d.all.ivs[, final.ivs:=ifelse(clumped.iv=="", clumped.ivs, clumped.iv), by=seq_len(nrow(t2d.all.ivs))]
      
      forward <- wrap.MRClust(exposure.file=paste0(trait$PATH, trait$TRAIT, "_", pop, "_forMR.txt"),  
                         outcome.file=paste0(trait$PATH, "t2d_", pop, "_forMR.txt"), 
                         exposure.name=trait$TRAIT,
                         outcome.name="t2d",
                         ivs=data.table::fread(paste0(trait$PATH, trait$TRAIT, "_", pop, "_ivs.txt"))$clumped.ivs,
						 savepath = paste0(trait$PATH),
             ancest = pop)
      print("Forward Finished")
	  t2d.ivs <- t2d.all.ivs$final.ivs
	  backward <- wrap.MRClust(exposure.file=paste0(trait$PATH, "t2d_", pop, "_forMR.txt"),  
                             outcome.file=paste0(trait$PATH, trait$TRAIT, "_", pop, "_forMR.txt"),
                             exposure.name="t2d",
                             outcome.name=trait$TRAIT,
                             ivs=t2d.ivs,
							 savepath = paste0(trait$PATH),
               ancest = pop)
	  print("Backward Finished")
							 }					 
}
}

################################################
#------------------- Testing -------------------
################################################

#Important Variable
setwd("/mnt/data1/home/davis/T2DGGI_MR/Comorbid/")
info.datasets <- data.table::as.data.table(readxl::read_excel("remaining_datasets.xlsx"))

#Test with AD
RunMRClust(info.datasets[15])

#Test with FDR Significant
RunMRClust(info.datasets)