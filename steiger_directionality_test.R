# https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#mr-steiger-directionality-test
# https://mrcieu.github.io/TwoSampleMR/reference/get_r_from_lor.html
# https://mrcieu.github.io/TwoSampleMR/reference/mr_steiger.html

################################################################################
#---------------------------------  Variables ----------------------------------
################################################################################
.libPaths(c("/home/itg/ana.arruda/R/x86_64-pc-linux-gnu-library/4.1",
            "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/Rpackages"))
pathMR="/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/"
setwd(pathMR)

################################################################################
#-----------------------------  forward direction ------------------------------
################################################################################
datasets.dt <- data.table::as.data.table(readxl::read_excel("info_datasets.xlsx"))
results <- data.table::fread("ResultsMR/new/all_analysis_results_w_meta_analysis.tsv")
fdr.sig.results <- results[p_value_fdr<0.05 & ancestry!="Meta_Analysis"]
prevalence.dt <- data.table::data.table(disease=unique(fdr.sig.results$outcome),
                                        prev=0.1)
res.dt <- data.table::data.table()
for(i in 1:nrow(fdr.sig.results)){
  pop <- fdr.sig.results[i, ancestry]
  c <- fdr.sig.results[i, cluster]
  t <- fdr.sig.results[i, outcome]
  trait <- datasets.dt[ANCESTRY==pop & TRAIT==t]
  
  t2d.all.ivs <- data.table::fread(paste0(trait$PATH, "t2d_", pop, "_ivs.txt"))
  if("clumped.iv" %in% colnames(t2d.all.ivs)) t2d.all.ivs[, final.ivs:=ifelse(clumped.iv=="", clumped.ivs, clumped.iv), by=seq_len(nrow(t2d.all.ivs))]
  else t2d.all.ivs[, final.ivs:=clumped.ivs]
  exposure <- data.table::fread(paste0(trait$PATH, "t2d_", pop, "_forMR.txt"))
  outcome <- data.table::fread(paste0(trait$PATH, trait$TRAIT, "_", pop, "_forMR.txt"))
  if(any(is.na(outcome$n))) outcome$n=rep(as.integer(datasets.dt[ANCESTRY==pop & TRAIT==t, N]), nrow(outcome))
  if(any(is.na(outcome$ncases))) outcome$ncases=rep(as.integer(datasets.dt[ANCESTRY==pop & TRAIT==t, NCASES]), nrow(outcome))
  if(any(is.na(outcome$effect_allele_frequency))){
    outcome$effect_allele_frequency <- NULL
    outcome <- merge(outcome, exposure[, .(rs_id, effect_allele_frequency, A1=effect_allele)], by="rs_id")
    outcome[, effect_allele_frequency:=ifelse(effect_allele==A1, effect_allele_frequency, 1-effect_allele_frequency)]
    outcome[, A1:=NULL]
  }
    
  if(c=="Full") t2d_ivs <- t2d.all.ivs[final.ivs %in% outcome$rs_id & final.ivs %in% exposure$rs_id, final.ivs]
  else t2d_ivs <- t2d.all.ivs[`Cluster assignment`==c & final.ivs %in% outcome$rs_id & final.ivs %in% exposure$rs_id, final.ivs]
  
  if(sum(is.na(outcome$effect_allele_frequency))>0) {
    outcome[, effect_allele_frequency:=NULL]
    outcome <- merge(outcome, exposure[, .(rs_id, effect_allele_frequency)], by="rs_id")
  }
  outcome <- unique(outcome, by="rs_id")
  
  est.r.exp <- TwoSampleMR::get_r_from_lor(lor=exposure[rs_id %in% t2d_ivs, beta],
                                           af=exposure[rs_id %in% t2d_ivs, effect_allele_frequency],
                                           ncase=exposure[rs_id %in% t2d_ivs, ncases],
                                           ncontrol=exposure[rs_id %in% t2d_ivs, n]-exposure[rs_id %in% t2d_ivs, ncases],
                                           prevalence=0.11)
  est.r.out <- TwoSampleMR::get_r_from_lor(lor=outcome[rs_id %in% t2d_ivs, beta],
                                           af=outcome[rs_id %in% t2d_ivs, effect_allele_frequency],
                                           ncase=outcome[rs_id %in% t2d_ivs, ncases],
                                           ncontrol=outcome[rs_id %in% t2d_ivs, n]-outcome[rs_id %in% t2d_ivs, ncases],
                                           prevalence=prevalence.dt[disease==t, prev])
  res <- TwoSampleMR::mr_steiger(p_exp = exposure[rs_id %in% t2d_ivs, p_value], 
                                 p_out = outcome[rs_id %in% t2d_ivs, p_value], 
                                 n_exp = exposure[rs_id %in% t2d_ivs, n], 
                                 n_out = outcome[rs_id %in% t2d_ivs, n], 
                                 r_exp = est.r.exp,
                                 r_out = est.r.out)
  res.dt <- rbind(res.dt, data.table::data.table(disease=t, cluster=c, ancestry=pop, correct_causal_direction=res$correct_causal_direction))
}

data.table::fwrite(res.dt, "ResultsMR/new/result_directionality_test_ana.tsv")

################################################################################
#-----------------------------  reverse direction ------------------------------
################################################################################
results <- data.table::fread("ResultsMR/new/reverse_results.tsv")
fdr.sig.results <- results[p_value_fdr<0.05 & ancestry!="Meta_Analysis"]
prevalence.dt <- data.table::data.table(disease=unique(fdr.sig.results$exposure),
                                        prev=0.1)

res.dt <- data.table::data.table()
for(i in 1:nrow(fdr.sig.results)){
  pop <- fdr.sig.results[i, ancestry]
  c <- fdr.sig.results[i, cluster]
  t <- fdr.sig.results[i, exposure]
  trait <- datasets.dt[ANCESTRY==pop & TRAIT==t]
  
  ivs <- data.table::fread(paste0(trait$PATH, t, "_", pop, "_ivs.txt"))
  ivs[, final.ivs:=clumped.ivs]
  outcome <- data.table::fread(paste0(trait$PATH, "t2d_", pop, "_forMR.txt"))
  exposure <- data.table::fread(paste0(trait$PATH, trait$TRAIT, "_", pop, "_forMR.txt"))
  if(any(is.na(exposure$n))) exposure$n=rep(as.integer(datasets.dt[ANCESTRY==pop & TRAIT==t, N]), nrow(exposure))
  if(any(is.na(exposure$ncases))) exposure$ncases=rep(as.integer(datasets.dt[ANCESTRY==pop & TRAIT==t, NCASES]), nrow(exposure))
  if(any(is.na(exposure$effect_allele_frequency))){
    exposure$effect_allele_frequency <- NULL
    exposure <- merge(exposure, outcome[, .(rs_id, effect_allele_frequency, A1=effect_allele)], by="rs_id")
    exposure[, effect_allele_frequency:=ifelse(effect_allele==A1, effect_allele_frequency, 1-effect_allele_frequency)]
    exposure[, A1:=NULL]
  }
  
  ivs <- ivs[final.ivs %in% outcome$rs_id & final.ivs %in% exposure$rs_id, final.ivs]
  est.r.exp <- TwoSampleMR::get_r_from_lor(lor=exposure[rs_id %in% ivs, beta],
                                           af=exposure[rs_id %in% ivs, effect_allele_frequency],
                                           ncase=exposure[rs_id %in% ivs, ncases],
                                           ncontrol=exposure[rs_id %in% ivs, n]-exposure[rs_id %in% ivs, ncases],
                                           prevalence=prevalence.dt[disease==t, prev])
  est.r.out <- TwoSampleMR::get_r_from_lor(lor=outcome[rs_id %in% ivs, beta],
                                           af=outcome[rs_id %in% ivs, effect_allele_frequency],
                                           ncase=outcome[rs_id %in% ivs, ncases],
                                           ncontrol=outcome[rs_id %in% ivs, n]-outcome[rs_id %in% ivs, ncases],
                                           prevalence=0.11)
  res <- TwoSampleMR::mr_steiger(p_exp = exposure[rs_id %in% ivs, p_value], 
                                 p_out = outcome[rs_id %in% ivs, p_value], 
                                 n_exp = exposure[rs_id %in% ivs, n], 
                                 n_out = outcome[rs_id %in% ivs, n], 
                                 r_exp = est.r.exp,
                                 r_out = est.r.out)
  res.dt <- rbind(res.dt, data.table::data.table(disease=t, cluster=c, ancestry=pop, correct_causal_direction=res$correct_causal_direction))
}

data.table::fwrite(res.dt, "ResultsMR/new/result_directionality_test_reverse_direction.tsv")




