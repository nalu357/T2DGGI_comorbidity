.libPaths(c(.libPaths(), "/home/itg/ana.arruda/R/x86_64-pc-linux-gnu-library/4.1",
            "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/Rpackages"))
library(dplyr)
library(data.table)

################################################
#------------------- Functions -----------------
################################################
local.clump <- function(data, pop="EUR", pval=5e-08){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ieugwasr::ld_clump(dplyr::tibble(rsid=data$SNP, pval=data$pval.exposure), #id=data$id.exposure),
                                   plink_bin=genetics.binaRies::get_plink_binary(),
                                   bfile=paste0("/lustre/groups/itg/shared/referenceData/1kG/EAS-SAS-AMR-EUR-AFR_1kg_v3/", pop),
                                   clump_p=pval)
  return(data[data$SNP %in% clumped.dt$rsid, ])
}

read.adiposity <- function(trait, project_folder){
  if (trait %in% c("bmi", "whr")){
    dt <- data.table::fread(paste0(project_folder, "data/", trait, ".giant-ukbb.meta-analysis.combined.23May2018.txt.gz"))
    return(dt[, .(chr=CHR, pos=POS, rs_id=sub(":.*", "", SNP), effect_allele=Tested_Allele, other_allele=Other_Allele, 
                  effect_allele_frequency=Freq_Tested_Allele, beta=BETA,standard_error=SE, p_value=P, n=N)])  
  }
  else if (trait %in% c("body_fat_percentage", "whole_body_fat_mass")){
    dt <- data.table::fread(paste0(project_folder, "data/", trait, "_neale_rsid.txt.gz"))
    return(dt[, .(chr=CHR, pos=POS, rs_id=SNP, effect_allele=minor_allele, effect_allele_frequency=minor_AF, 
                  other_allele=major_allele, beta, standard_error=se, p_value=pval, n=n_complete_samples)])  
  }
  else if (trait %in% c("ldl", "tgs", "hdl", "apoa", "apob")){
    dt <- data.table::fread(paste0(project_folder, "data/", trait, "_neale_rsid.txt.gz"))
    return(dt[, .(chr=CHR, pos=POS, rs_id=rsID, effect_allele=minor_allele, effect_allele_frequency=minor_AF, 
                  other_allele=major_allele, beta, standard_error=se, p_value, n=n_complete_samples)])  
  }
  else {
    return(print("Trait not recognized.")) 
  }
}

read.magic <- function(trait, project_folder){
  dt <- data.table::fread(paste0(project_folder, "data/MAGIC1000G_", trait, "_EUR.tsv.gz"))
  return(dt[, .(chr=chromosome, pos=base_pair_location, rs_id=variant, effect_allele, other_allele, 
                effect_allele_frequency, beta, standard_error, p_value, n=sample_size)])
}

get.filename <- function(trait, project_folder){
  if(trait %in% c("2hGlu", "FI", "FG", "HbA1c")){
    dt <- read.magic(trait, project_folder)
  }
  else if(trait %in% c("bmi", "whr", "whole_body_fat_mass", "body_fat_percentage", "ldl", "tgs", "hdl", "apoa", "apob")){
    dt <- read.adiposity(trait, project_folder)
  }
  else print("ERROR: Subset not known")
  
  return(dt)
}

get.t2d.data <- function(info.t2d, path_T2DGGI, ancestry_T2DGGI){
  info.t2d.ancestry <- subset(info.t2d, ancestry == ancestry_T2DGGI)
  t2d <- data.table::fread(paste0(path_T2DGGI, ancestry_T2DGGI, "_MetalFixed_LDSC-CORR_Results1TBL_rsid.gz"))
  t2d <- t2d[, .(chr.b38=as.integer(sub("chr", "", CHR_HG38)),
                 pos.b38=POS_HG38,
                 chr.b37=as.integer(sub("chr", "", CHR_HG19)),
                 pos.b37=POS_HG19, 
                 rs_id=RSID,
                 p_value=as.numeric(`P-value`), 
                 beta=Effect,
                 standard_error=StdErr, 
                 effect_allele_frequency=Freq1, 
                 n=info.t2d.ancestry$n.cases+info.t2d.ancestry$n.controls, 
                 ncases=info.t2d.ancestry$n.cases, 
                 effect_allele=toupper(Allele1), 
                 other_allele=toupper(Allele2))]
  return(t2d)
}

get_trait_data <- function(trait, tmp_data, type="binary"){
  dt <- data.table::fread(trait$FILE, tmpdir=tmp_data)
  if(type=="binary"){
    dt <- dt[, .(chr=ifelse(is.na(trait$chr.col), NA, as.integer(sub("chr", "", get(trait$chr.col)))),
                 pos=ifelse(is.na(trait$pos.col), NA, get(trait$pos.col)),
                 p_value=get(trait$pval.col), 
                 beta=get(trait$beta.col),
                 standard_error=get(trait$se.col), 
                 rs_id=sub(":.*", "", get(trait$rsid.col)), 
                 effect_allele=toupper(get(trait$ea.col)), 
                 other_allele=toupper(get(trait$nea.col)),
                 effect_allele_frequency=ifelse(is.na(trait$eaf.col), NA, get(trait$eaf.col)), 
                 n=ifelse(is.na(trait$n.col), trait$N, get(trait$n.col)),
                 ncases=ifelse(is.na(trait$ncases.col), trait$NCASES, get(trait$ncases.col)))]
  }
  else {
    dt <- dt[, .(chr=ifelse(is.na(trait$chr.col), NA, as.integer(sub("chr", "", get(trait$chr.col)))),
                 pos=ifelse(is.na(trait$pos.col), NA, get(trait$pos.col)),
                 p_value=get(trait$pval.col), 
                 beta=get(trait$beta.col),
                 standard_error=get(trait$se.col), 
                 rs_id=sub(":.*", "", get(trait$rsid.col)), 
                 effect_allele=toupper(get(trait$ea.col)), 
                 other_allele=toupper(get(trait$nea.col)),
                 effect_allele_frequency=ifelse(is.na(trait$eaf.col), NA, get(trait$eaf.col)), 
                 n=ifelse(is.na(trait$n.col), trait$N, get(trait$n.col)))]
  }
  
  #Keep only SNVs with rsID
  dt <- dt[effect_allele %in% c("A", "T", "C", "G")]
  dt <- dt[other_allele %in% c("A", "T", "C", "G")]
  dt<- subset(dt, !is.na(rs_id))
  dt<- subset(dt, !grepl(",", rs_id))
  dt<- dt[rs_id!="."]
  dt<- dt[rs_id!=""]
  
  return(dt)
}

get.data <- function(file, data.type, data.dir){
  if(data.type=="quant"){
    data <- TwoSampleMR::format_data(type = data.dir,
                                     dat = file,
                                     snp_col = "rs_id",
                                     beta_col = "beta",
                                     se_col = "standard_error",
                                     effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",
                                     eaf_col = "effect_allele_frequency",
                                     pval_col = "p_value",
                                     chr_col="chr",
                                     pos_col="pos",
                                     samplesize_col="n")
  }
  else if(data.type=="binary"){
    data <- TwoSampleMR::format_data(type = data.dir,
                                     dat = file,
                                     snp_col = "rs_id",
                                     beta_col = "beta",
                                     se_col = "standard_error",
                                     effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",
                                     eaf_col = "effect_allele_frequency",
                                     pval_col = "p_value",
                                     chr_col="chr",
                                     pos_col="pos",
                                     ncase_col="ncases",
                                     samplesize_col="n")
  }
  else {
    print("Data type not recognized (options: quant or binary).")
  }
  
  return(data)
}

run.mvmregger <- function(data, phen){
  #----------- Run MV MR-Egger from MVMR package -----------
  mv_exp_beta <- as.data.frame(data$exposure_beta)
  mv_exp_se <- as.data.frame(data$exposure_se)
  for(i in 1:length(colnames(mv_exp_beta))){
    assign(paste0("exposure", i), data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure)
    colnames(mv_exp_beta)[1] <- paste0(data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure, "_beta")
    colnames(mv_exp_se)[1] <- paste0(data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure, "_se")
  }
  
  mv_out_beta <- as.data.frame(data$outcome_beta)
  mv_out_se <- as.data.frame(data$outcome_se)
  colnames(mv_out_beta) <- c(paste(phen, "beta", sep="_"))
  colnames(mv_out_se) <- c(paste(phen, "se", sep="_"))
  
  mv_exp_beta$SNP <- row.names(mv_exp_beta)
  mv_exp_se$SNP <- row.names(mv_exp_se)
  
  mv_exp <- merge(mv_exp_beta, mv_exp_se, by="SNP")
  mv_exp <- cbind(mv_exp, mv_out_beta)
  mv_exp <- cbind(mv_exp, mv_out_se)
  
  F.data <- MVMR::format_mvmr(BXGs = mv_exp[,c(2,3)],
                              BYG = mv_exp[,6],
                              seBXGs = mv_exp[,c(4,5)],
                              seBYG = mv_exp[,7],
                              RSID = mv_exp[,1])
  
  sres <- MVMR::strength_mvmr(r_input = F.data, gencov = 0)
  pres <- MVMR::pleiotropy_mvmr(r_input = F.data, gencov = 0)
  res_mvmr <- MVMR::ivw_mvmr(r_input = F.data)
  res_mvmr2 <- MVMR::qhet_mvmr(F.data, matrix(c(1,0,0,1), nrow=2, ncol=2), CI = F, iterations = 100)
  
  dt.egger <- merge(data.table::as.data.table(res_mvmr) %>% .[, .(exposure=c(exposure1, exposure2), effect=Estimate, se=`Std. Error`, tvalue=`t value`, `Pr(>|t|)`)],
                    data.table::as.data.table(res_mvmr2) %>% .[, .(exposure=c(exposure1, exposure2), effect.qhet=`Effect Estimates`)],
                    by="exposure") %>% 
    .[, `:=`(outcome=phen,Fstat=paste(sres$exposure1, sres$exposure2, sep="_"), Qstat=pres$Qstat, Qpval=pres$Qpval)]
  return(dt.egger)
}

run.mvmr <- function(data, phen, filename, egger){
  if(nrow(data$exposure_beta)<2) {
    print("After harmonizing data, not enough SNPs left :(")
    dt <- data.table::data.table()
    return(dt)
  }
  else {
    print("Multiple SNPs in common between exposure and outcome, running MVMR")
    
    # Assumption 1: check strength of IVs with F-statistic and remove the ones with Fstat<10
    # data$fstat_per_snp <- (data$beta.exposure)^2/(data$se.exposure)^2
    # data <- data[data$fstat_per_snp>10,]
    
    #----------- Run MV TwoSampleMR ivw -----------
    res <- TwoSampleMR::mv_multiple(data, pval_threshold = 5e-01, plots=TRUE)
    
    #------------ Save output and plots
    dt <- data.table::as.data.table(res$result[,-c(1,3)])
    
    if(egger & nrow(data$expname)==2){
      #----------- Run MV MR-Egger from MVMR package -----------
      mv_exp_beta <- as.data.frame(data$exposure_beta)
      mv_exp_se <- as.data.frame(data$exposure_se)
      for(i in 1:length(colnames(mv_exp_beta))){
        assign(paste0("exposure", i), data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure)
        colnames(mv_exp_beta)[1] <- paste0(data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure, "_beta")
        colnames(mv_exp_se)[1] <- paste0(data$expname[data$expname$id.exposure==colnames(mv_exp_beta)[i],]$exposure, "_se")
      }
      
      mv_out_beta <- as.data.frame(data$outcome_beta)
      mv_out_se <- as.data.frame(data$outcome_se)
      colnames(mv_out_beta) <- c(paste(phen, "beta", sep="_"))
      colnames(mv_out_se) <- c(paste(phen, "se", sep="_"))
      
      mv_exp_beta$SNP <- row.names(mv_exp_beta)
      mv_exp_se$SNP <- row.names(mv_exp_se)
      
      mv_exp <- merge(mv_exp_beta, mv_exp_se, by="SNP")
      mv_exp <- cbind(mv_exp, mv_out_beta)
      mv_exp <- cbind(mv_exp, mv_out_se)
      
      F.data <- MVMR::format_mvmr(BXGs = mv_exp[,c(2,3)],
                                  BYG = mv_exp[,6],
                                  seBXGs = mv_exp[,c(4,5)],
                                  seBYG = mv_exp[,7],
                                  RSID = mv_exp[,1])
      sres <- MVMR::strength_mvmr(r_input = F.data, gencov = 0)
      pres <- MVMR::pleiotropy_mvmr(r_input = F.data, gencov = 0)
      res_qhet <- MVMR::qhet_mvmr(F.data, matrix(c(1,0,0,1), nrow=2, ncol=2), CI = F, iterations = 100)
      dt$effect.qhet <- res_qhet$`Effect Estimates`
      dt$Fstat <- c(sres$exposure1, sres$exposure2)
      dt$Qstat <- pres$Qstat
      dt$Qpval <- pres$Qpval
    }
  }
  return(dt)
}

get_odds_ratios_mvmr <- function(mr_res){
  mr_res$or <- exp(mr_res$b)
  mr_res$or_lci95 <- exp(mr_res$b - 1.96 * mr_res$se)
  mr_res$or_uci95 <- exp(mr_res$b + 1.96 * mr_res$se)
  return(mr_res)
}

wrap.MVMR <- function(exposure.lst, outcome.lst, outcome.data, exposure1.data, exposure2.data, output.file, ivs, custom.clump.p=5e-08, egger=TRUE){
  # exposure.lst = list of lists of vectors
  dt <- data.table::data.table()
  exposure1.data <- data.frame(exposure1.data)
  exposure1 <- get.data(exposure1.data, data.dir="exposure", data.type=exposure.lst[[1]][3])
  exposure1$exposure <- paste(exposure.lst[[1]][1], exposure.lst[[1]][2], sep=".")
  
  for (exp in exposure.lst[-1]){
    #----------- Read and clump exposure -----------
    exposure2.data <- data.frame(exposure2.data)
    exposure2 <- get.data(exposure2.data, data.dir="exposure", data.type=exp[3])
    exposure2$exposure <- paste(exp[1], exp[2], sep=".")
    
    clump.exposure <- dplyr::bind_rows(exposure1, exposure2)
  }
  
  for (out in outcome.lst) {
    #----------- Get outcome -----------
    outcome.data <- data.frame(outcome.data)
    outcome <- get.data(outcome.data, data.dir="outcome", data.type=out[3])
    outcome$outcome <- paste(out[1], out[2], sep=".")
    
    #----------- Harmonize data -----------
    data <- TwoSampleMR::mv_harmonise_data(exposure_dat=clump.exposure, outcome_dat=outcome)
    
    #----------- Run TwoSampleMR for unfiltered data -----------
    dt <- rbind(dt, run.mvmr(data, filename=output.file, phen=out[1], egger=egger), fill=TRUE)
  }
  
  dt <- dt[!is.na(exposure)]
  dt <- get_odds_ratios_mvmr(dt)
  
  data.table::fwrite(dt, paste0("mvmr/", output.file, ".csv"))
  return(dt)
}

########################################################################################
#---------------------------- Get IVs for cardiometabolic traits ----------------------#
########################################################################################
cm.datasets.dt <- data.table::as.data.table(readxl::read_excel("cluster_traits.xlsx"))

ivs.dt<- data.table::data.table()
for (t in c("TG", "HDL", "SAT", "VAT", "GFAT", "PI", "FI", "FG","RG", "Glucose", "bmi", "whr", "body_fat_percentage")){
  for(pop in unique(cm.datasets.dt[TRAIT==t, ANCESTRY])){
    #----------- Read and clump exposure -----------
    trait <- cm.datasets.dt[ANCESTRY==pop & TRAIT==t]
    exposure <- get.data(get_trait_data(trait, tmp_data=tmp_dir), data.dir="exposure", data.type="quant")
    exposure$exposure <- paste(t, pop, sep=".")
    clump.exp <- local.clump(data=exposure, pop=pop, pval=5e-8)
    
    #Keep only F-stat>10
    exposure$snp_fstat=(exposure$beta.exposure)^2/(exposure$se.exposure)^2
    
    #Only keep the variants that are clumped
    exposure$iv=ifelse(exposure$SNP %in% clump.exp$SNP & exposure$snp_fstat>=10, TRUE, FALSE)
    ivs <- exposure[exposure$iv==TRUE,]$SNP
    
    ivs.dt <- rbind(ivs.dt, data.table::data.table(trait=paste(t, pop, sep="_"), snp=ivs))
  }
}
data.table::fwrite(ivs.dt, "mvmr/all_cluster_traits_ivs.csv")

########################################################################################
#------------------------------------ Variables ---------------------------------------#
########################################################################################
pathMR="/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/"
path_T2DGGI = "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/datasets/T2DGGI_bothbuilds/"
tmp_dir = "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/tmpdata/"
setwd(pathMR)

cm.ivs <- data.table::fread("/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/mvmr/all_cluster_traits_ivs.csv")

info.T2DGGI <- data.frame(ancestry = c("AFR", "EAS", "EUR", "AMR", "SAS"),
                          n.cases = c(50251, 88109, 242283, 29375, 16832),
                          n.controls = c(103909, 339395, 1569734, 59368, 33767))

datasets.dt <- data.table::as.data.table(readxl::read_excel("info_datasets.xlsx"))
cm.datasets.dt <- data.table::as.data.table(readxl::read_excel("cluster_traits.xlsx"))
results <- unique(data.table::fread("ResultsMR/new/all_analysis_results_w_meta_analysis.tsv"))
fdr.sig.results <- results[p_value_fdr<0.05]

########################################################################################
#------------------------------------ Running -----------------------------------------#
########################################################################################
cluster.traits <- data.table::data.table(cluster=c("Obesity", 
                                                   "Beta cell +PI", 
                                                   "Residual glycaemic", 
                                                   "Beta cell -PI", 
                                                   "Lipodystrophy", 
                                                   "Metabolic syndrome", 
                                                   "Body fat",
                                                   "Full"),
                                         traits=c(list(c("bmi", "whr", "body_fat_percentage", "HDL")), 
                                                  list(c("FI", "FG", "PI", "Glucose", "bmi")),  # "FIadjBMI", "FGadjBMI" 
                                                  list(c("FI", "FG", "PI", "Glucose", "bmi")),  # "FIadjBMI", "FGadjBMI" 
                                                  list(c("FI", "FG", "PI", "Glucose", "bmi")),  # "FIadjBMI", "FGadjBMI" 
                                                  list(c("TG", "FI", "whr", "HDL", "body_fat_percentage", "GFAT", "bmi")),  # "FIadjBMI", 
                                                  list(c("FG", "whr", "TG", "FI", "HDL", "VAT", "SAT", "Glucose", "bmi")),  #"FIadjBMI", "FGadjBMI" 
                                                  list(c("body_fat_percentage", "bmi", "SAT", "VAT")),
                                                  list(c("bmi", "whr", "body_fat_percentage", "HDL", "FI", "FG", "PI", "Glucose", "VAT", "SAT", "GFAT", "TG"))))
univariate.cm.mr <- data.table::fread("mvmr/processed_univariate_mr_results_cm_traits.csv")
sig.univariate.cm.mr <- univariate.cm.mr[p_value_fdr<0.05]

res.dt <- data.table::data.table()
for(pop in unique(fdr.sig.results$ancestry)){
  if(pop=="Meta_Analysis") next
  
  t2d <- get.t2d.data(info.T2DGGI, path_T2DGGI, pop)
  
  for(t in unique(fdr.sig.results[ancestry==pop, outcome])){
    cur.results <- fdr.sig.results[ancestry==pop & outcome==t]
    trait.dt <- datasets.dt[ANCESTRY==pop & TRAIT==t]
    
    t2d.all.ivs <- data.table::fread(paste0(trait.dt$PATH, "t2d_", pop, "_ivs.txt"))
    if("clumped.iv" %in% colnames(t2d.all.ivs)) t2d.all.ivs[, final.ivs:=ifelse(clumped.iv=="", clumped.ivs, clumped.iv), by=seq_len(nrow(t2d.all.ivs))]
    else t2d.all.ivs[, final.ivs:=clumped.ivs]
    
    cur.trait <- get_trait_data(trait.dt, tmp_data=tmp_dir)
    
    for(i in 1:nrow(cur.results)){
      c <- cur.results[i, cluster]
      if(c!="Full") {
        t2d_ivs <- t2d.all.ivs[`Cluster assignment`==c, final.ivs] 
      }
      else {
        t2d_ivs <- t2d.all.ivs$final.ivs
      }
      
      if(length(sig.univariate.cm.mr[outcome==t & ancestry==pop, unique(exposure)])==0) next
      for(exp2 in sig.univariate.cm.mr[outcome==t & ancestry==pop, unique(exposure)]) { 
        if(nrow(cm.ivs[trait==paste(exp2, pop, sep="_")])==0) next
        cur.ivs <- unique(c(t2d_ivs, cm.ivs[trait %in% paste(exp2, pop, sep="_"), snp]))
        cur.ivs <- local.clump(data=data.table::data.table(SNP=cur.ivs, pval.exposure=0), pop=pop, pval=1)
        cur.ivs <- cur.ivs$SNP
        exposure2 <- get_trait_data(cm.datasets.dt[ANCESTRY==pop & TRAIT==exp2], tmp_data=tmp_dir, type="quant")
        exposure2 <- exposure2[rs_id %in% cur.ivs]
        cur.res <- wrap.MVMR(exposure.lst=list(c(paste0("t2d.", gsub(" ", "_", c)), pop, "binary"), c(exp2, pop, "quant")),
                             outcome.lst=list(c(t, pop, "binary")),
                             exposure1.data=t2d[rs_id %in% cur.ivs],
                             exposure2.data=exposure2,
                             outcome.data=cur.trait[rs_id %in% cur.ivs],
                             ivs=cur.ivs,
                             output.file=paste0("mvmr_", t, "_", pop, "_", gsub(" ", "_", c), "_", exp2))
        cur.res$cluster <- c
        res.dt <- rbind(res.dt, cur.res)
        
        rm(exposure2)
      }
    }
    rm(cur.trait)
  }
}

data.table::fwrite(res.dt, "mvmr/mvmr_results_new_approach.csv")