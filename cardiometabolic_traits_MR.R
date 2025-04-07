library(dplyr)
library(data.table)

################################################
#------------------- Functions -----------------
################################################
local.clump <- function(data, pop="EUR", pval=5e-08){
  # https://mrcieu.github.io/ieugwasr/articles/local_ld.html
  clumped.dt <- ieugwasr::ld_clump(dplyr::tibble(rsid=data$rs_id, pval=data$p_value),
                                   plink_bin=genetics.binaRies::get_plink_binary(),
                                   bfile=paste0("/lustre/groups/itg/shared/referenceData/1kG/EAS-SAS-AMR-EUR-AFR_1kg_v3/", pop),
                                   clump_p=pval)
  return(data[rs_id %in% clumped.dt$rsid])
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

add_odds_ratios <- function(dt, beta="beta", se="standard_error"){
  dt$or <- exp(dt[, get(beta)])
  dt$or_lci95 <- exp(dt[, get(beta)] - 1.96 * dt[, get(se)])
  dt$or_uci95 <- exp(dt[, get(beta)] + 1.96 * dt[, get(se)])
  return(dt)
}

run.SensitivityAnalysis <- function(data, res, plt=TRUE, het=TRUE, steiger = TRUE){
  fstat_overall <- mean((data$beta.exposure)^2/(data$se.exposure)^2)
  res[, Fstat:=fstat_overall]
  
  if(het){
    # check heterogeneity with Cochran's Q-statistic
    het <- data.table::as.data.table(TwoSampleMR::mr_heterogeneity(data)) 
    res <- merge(res, het[, .(exposure, outcome, method, Q, Q_df, Q_pval)],
                 by=c("exposure", "outcome", "method"), all.x=TRUE)
  }
  if(plt){
    # Assumption 2: check pleiotropy with MR-Egger intercept
    plt <- data.table::as.data.table(TwoSampleMR::mr_pleiotropy_test(data))  #MR Egger intercept for directional pleiotropy
    res[method=="MR Egger", `:=` (egger_intercept=plt$egger_intercept, egger_intercept_se=plt$se, egger_intercept_pval=plt$pval)]
  }
  if(steiger){
    #Run the Steiger filtering
    data.steiger <- tryCatch(
      TwoSampleMR::steiger_filtering(data),
      error=function(e) e
    )
    
    if(inherits(data.steiger, "error")) print("Steiger-filter yielded error")
    else{
      data.steiger <- subset(data.steiger, steiger_dir & steiger_pval<0.05)
      res.steiger <- TwoSampleMR::mr(data.steiger, method = "mr_ivw")
      if(nrow(res.steiger)!=0){
        res.steiger <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res.steiger))
        res.steiger[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
        #Compute heterogeneity
        het <- data.table::as.data.table(TwoSampleMR::mr_heterogeneity(data.steiger, method = "mr_ivw"))
        res.steiger <- merge(res.steiger, het[, .(exposure, outcome, method, Q, Q_df, Q_pval)],
                             by=c("exposure", "outcome", "method"), all.x=TRUE)
        res.steiger$method <- "Steiger Inverse variance weighted"
        res <- rbind(res, res.steiger, fill = TRUE)
      }
    }
  }
  
  return(res)
}

run.TwoSampleMR <- function(data){
  if(nrow(data)==0)  {
    print("No SNP in common between exposure and outcome :(")
    tmp.dt <- data.table::data.table()
    return(tmp.dt)
  }
  else if (nrow(data)==1) {
    fstat <- (data$beta.exposure)^2/(data$se.exposure)^2
    print("Only one SNP in common between exposure and outcome, running Wald ratio method")
    res <- TwoSampleMR::mr(data, method_list="mr_wald_ratio")
    res <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res))
    res[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
    
    tmp.dt <- data.table::as.data.table(res) %>% .[, `:=`(Fstat=fstat,Q=as.numeric(NA),Q_df=as.numeric(NA),Q_pval=as.numeric(NA),plt.egger_intercept=as.numeric(NA),plt.pval=as.numeric(NA),plt.se=as.numeric(NA))]
    return(tmp.dt)
  }
  else {
    print("Multiple SNPs in common between exposure and outcome, running MR")
    
    res <- TwoSampleMR::mr(data, method_list=c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression", "mr_penalised_weighted_median", "mr_sign"))
    res <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res))
    res[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
    
    # Run sensitivity analysis
    res <- run.SensitivityAnalysis(data, res)
    
    return(res)
  }
}

wrap.MR <- function(exposure.file, outcome.file, exposure.name, outcome.name, ivs){
  if(length(ivs)==0) return(NULL)
  
  #----------- Get exposure -----------
  exposure <- get.data(data.table::fread(exposure.file,data.table=FALSE), data.type="quant", data.dir="exposure")
  exposure <- exposure[exposure$SNP %in% ivs,]
  if(nrow(exposure)==0) return(res.dt)
  exposure$exposure <- exposure.name
  
  #----------- Get outcome -----------
  outcome <- get.data(data.table::fread(outcome.file,data.table=FALSE), data.type="binary", data.dir="outcome")
  outcome$outcome <- outcome.name
  
  #----------- Harmonize data -----------
  data <- tryCatch(
    TwoSampleMR::harmonise_data(exposure_dat=exposure, outcome_dat=outcome),
    error=function(e) e
  )
  
  if(inherits(data, "error")) {
    print("Harmonization failed.")
    return(NULL) 
  }
  
  if(nrow(data)==0) return(NULL)
  
  if(nrow(subset(data, mr_keep==TRUE))==0) {
    print("No SNP left after data harmonization, not able to run MR :(")
    return(NULL)
  }
  data$id.exposure <- data$exposure
  data$id.outcome <- data$outcome  
  data <- subset(data, mr_keep==TRUE)
  
  #----------- Run TwoSampleMR with sensitivity analyses -----------
  res.mr <- run.TwoSampleMR(data)
  
  #Compute the I2 statistics
  res.mr$I2 <- (res.mr$Q-res.mr$Q_df)/res.mr$Q
  
  return(res.mr)
}

function.process.cmMR <- function(res.mr, weighted.mode.sensitivity = F){
  #If no look at Radial filtering, remove it from results
  res.mr <- subset(res.mr, method != "Sign concordance test")
  
  #Agreed not to look at weighted mode
  if(!weighted.mode.sensitivity){
    res.mr <- subset(res.mr, method != "Weighted mode")
  }
  #Subset the results to keep the IVW p-value and effect size if at least 3 snps
  res.mr.IVW <- subset(res.mr, method == "Inverse variance weighted")
  #Only keep Wald Ratio if only one or two SNPs
  res.mr.wald <- subset(res.mr, method == "Wald ratio")
  
  res.mr.all <- data.table::data.table()
  
  ###########################Sensitivity analyses
  for(out in unique(res.mr$outcome)){
    for (pop in unique(res.mr[outcome==out, ancestry])){
      #Flag if the MR Egger intercept is significant
      res.mr.MREgger <- subset(res.mr, outcome==out & method == "MR Egger" & ancestry==pop)
      res.mr.MREgger$pval_fdr <- p.adjust(res.mr.MREgger$egger_intercept_pvalue, method = "fdr")
      MREgger.sig <- subset(res.mr.MREgger, pval_fdr < 0.05)$exposure
      #Look at whether heterogeneity is above threshold
      Het.sig <- subset(res.mr.IVW, outcome==out & I2>0.5 & ancestry==pop)$exposure
      #Flag reverse causation (Steiger filtered IVW has a different direction of effect as IVW)
      Steiger.diff <- res.mr[outcome==out & ancestry==pop & method %in% c("Steiger Inverse variance weighted","Inverse variance weighted"), length(unique(sign(beta)))>1, by=exposure] %>% .[V1==TRUE, exposure]
      
      tmp.res <- res.mr[outcome==out & ancestry==pop]
      if(nrow(tmp.res)==0) next
      beta.sensitivity <- data.frame(MREgger = ifelse(length(subset(tmp.res, method == "MR Egger")$beta)!=0, subset(tmp.res, method == "MR Egger")$beta, NA),
                                     WeightedMedian = ifelse(length(subset(tmp.res, method == "Weighted median")$beta)!=0, subset(tmp.res, method == "Weighted median")$beta, NA),
                                     Steiger = ifelse(length(subset(tmp.res, method == "Steiger Inverse variance weighted")$beta)!=0, subset(tmp.res, method == "Steiger Inverse variance weighted")$beta, NA))
      
      Prop.SameDir <- sum(sign(beta.sensitivity)!=sign(subset(res.mr.IVW, outcome==out & ancestry==pop)$beta), na.rm=TRUE)
      
      #Flag the results that don't have same direction
      res.mr.IVW[outcome==out & ancestry==pop, DiffDirection:=ifelse(Prop.SameDir!=0, TRUE, FALSE)]
      
      #Combine the two files
      res.mr.out <- rbind(res.mr.IVW[outcome==out & ancestry==pop], res.mr.wald[outcome==out & ancestry==pop], fill=TRUE)
      
      #Flag the results that have significant het, pleio
      res.mr.out$FlagPleiotropy <- ifelse(res.mr.out$exposure %in% unique(MREgger.sig), T, F)
      res.mr.out$FlagHeterogeneity <- ifelse(res.mr.out$exposure %in% unique(Het.sig), T, F)
      res.mr.out$FlagSteiger <- ifelse(res.mr.out$exposure %in% unique(Steiger.diff), T, F)
      
      res.mr.all <- rbind(res.mr.all, res.mr.out, fill=TRUE)
    }
  }
  
  #Compute the final ajusted p-value
  res.mr.all$p_value_fdr <- p.adjust(res.mr.all$p_value, method = "fdr")
  res.mr.all[, p_value_plot:=ifelse(p_value_fdr<0.05 & DiffDirection==FALSE & FlagHeterogeneity==FALSE & FlagPleiotropy==FALSE, 0, 2)]
  
  return(res.mr.all)                            
}

################################################
#------------------- Variables -----------------
################################################
.libPaths(c("/home/itg/ana.arruda/R/x86_64-pc-linux-gnu-library/4.1",
            "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/Rpackages"))
pathMR="/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/"
tmp_dir = "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/tmpdata/"
setwd(pathMR)

########################################################################################
#---------------------------- Get IVs for new traits ----------------------------------#
########################################################################################
cm.datasets.dt <- data.table::as.data.table(readxl::read_excel("cluster_traits.xlsx"))
cm.datasets.dt[ANCESTRY=="EUR" & TRAIT=="body_fat_percentage", FILE:="/lustre/groups/itg/teams/zeggini/projects/SCZ_T2D/data/body_fat_percentage_neale_rsid.txt.gz"]
cm.datasets.dt[ANCESTRY=="EUR" & TRAIT=="PI", FILE:="/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/data/cluster_traits/proinsulin_rsid.tsv.gz"]
cm.datasets.dt[ANCESTRY=="EUR" & TRAIT=="PI", rsid.col:="rs_id"]
ivs.dt<- data.table::data.table()
for (t in c("TG", "HDL", "SAT", "VAT", "GFAT", "PI", "FI", "FG", "RG", "Glucose", "bmi", "whr", "body_fat_percentage")){
  for(pop in unique(cm.datasets.dt[TRAIT==t, ANCESTRY])){
    #----------- Read and clump exposure -----------
    trait <- cm.datasets.dt[ANCESTRY==pop & TRAIT==t]
    exposure <- get_trait_data(trait, tmp_data=tmp_dir)
    clump.exp <- local.clump(data=exposure, pop=pop, pval=5e-8)
    
    #Keep only F-stat>10
    exposure[, snp_fstat:=(beta)^2/(standard_error)^2]
    
    #Only keep the variants that are clumped
    exposure[, iv:=ifelse(rs_id %in% clump.exp$rs_id & snp_fstat>=10, TRUE, FALSE)]
    
    ivs.dt <- rbind(ivs.dt, data.table::data.table(trait=paste(t, pop, sep="_"), snp=exposure[iv==TRUE, rs_id]))
    data.table::fwrite(exposure[iv==TRUE], paste0("mvmr/", t, "_", pop, "_forMR.txt"))
  }
}
data.table::fwrite(ivs.dt, "mvmr/all_cluster_traits_ivs.csv")

########################################################################################
#---------------------------- Prepare disease files -----------------------------------#
########################################################################################
cm.ivs <- data.table::fread("mvmr/all_cluster_traits_ivs.csv")

datasets.dt <- data.table::as.data.table(readxl::read_excel("info_datasets.xlsx"))
cm.datasets.dt <- data.table::as.data.table(readxl::read_excel("cluster_traits.xlsx"))
results <- unique(data.table::fread("ResultsMR/new/all_analysis_results_w_meta_analysis.tsv"))
fdr.sig.results <- results
pop="EUR"
unique(results$outcome)
t=""

for(pop in unique(fdr.sig.results$ancestry)){
  if(pop=="Meta_Analysis") next
  for(t in unique(fdr.sig.results[ancestry==pop, outcome])){
    trait <- datasets.dt[ANCESTRY==pop & TRAIT==t]
    dt <- get_trait_data(trait, tmp_dir)
    data.table::fwrite(dt[rs_id %in% unique(cm.ivs$snp)], paste0("mvmr/", t, "_", pop, "_cmMR.txt"))
  }
}

########################################################################################
#------------------------------------- Main -------------------------------------------#
########################################################################################
# Read CM traits pre-calculated clumped IVs
cm.ivs <- data.table::fread("mvmr/all_cluster_traits_ivs.csv")

datasets.dt <- data.table::as.data.table(readxl::read_excel("info_datasets.xlsx"))
cm.datasets.dt <- data.table::as.data.table(readxl::read_excel("cluster_traits.xlsx"))
cm.datasets.dt <- cm.datasets.dt[!(grep("FIadjBMI|FGadjBMI|hdl|tgs", TRAIT))]
results <- unique(data.table::fread("ResultsMR/new/all_analysis_results_w_meta_analysis.tsv"))
fdr.sig.results <- results[p_value_fdr<0.05]

cm.ivs[, tmp := sub("_(?=[^_]*$)", "@@", trait, perl = TRUE)]
cm.ivs[, c("cm.trait", "ancestry") := tstrsplit(trait, "@@")]
cm.ivs[, tmp:=NULL]

res.dt <- data.table::data.table()
for(pop in unique(fdr.sig.results$ancestry)){
  if(pop=="Meta_Analysis") next
  for(disease in unique(fdr.sig.results[ancestry==pop, outcome])){
    trait <- datasets.dt[ANCESTRY==pop & TRAIT==disease]
    
    for(t in unique(cm.datasets.dt[ANCESTRY==pop, TRAIT])){
      cur.res <- wrap.MR(exposure.file=paste0("mvmr/", t, "_", pop, "_forMR.txt"),  
                         outcome.file=paste0("mvmr/", disease, "_", pop, "_cmMR.txt"),
                         exposure.name=t,
                         outcome.name=disease,
                         ivs=cm.ivs[trait==paste(t, pop, sep="_"), snp])
      if(!is.null(cur.res)){
        cur.res$ancestry <- pop
        res.dt <- rbind(res.dt, cur.res, fill=TRUE)
      }
    }
  }
}
data.table::setnames(res.dt, old=c("b", "se", "pval", "or", "or_lci95", "or_uci95", "egger_intercept_se", "egger_intercept_pval"),
                     new=c("beta", "standard_error", "p_value", "odds_ratio", "ci_lower", "ci_upper", "egger_intercept_standard_error", "egger_intercept_pvalue"))
data.table::fwrite(res.dt, "mvmr/univariate_mr_results_cm_traits.csv")

# Process results
processed.dt <- function.process.cmMR(res.dt)
processed.dt[, `:=`(egger_intercept_pvalue=NULL, Distortion_pval=NULL, egger_intercept=NULL, egger_intercept_standard_error=NULL)]
processed.dt[p_value_fdr<0.05, .(exposure, outcome, ancestry, beta, p_value_fdr)]
data.table::fwrite(processed.dt, "mvmr/processed_univariate_mr_results_cm_traits.csv")
