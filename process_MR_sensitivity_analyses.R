################################################################################
#--------------------------------- FUNCTIONS ----------------------------------#
################################################################################
function.processMR <- function(res.mr, weighted.mode.sensitivity = F){
  #If no look at Radial filtering, remove it from results
  res.mr <- subset(res.mr, method != "Sign concordance test")
  
  #Agreed not to look at weighted mode
  if(!weighted.mode.sensitivity){
    res.mr <- subset(res.mr, method != "Weighted mode")
  }
  #Subset the results to keep the IVW p-value and effect size if at least 3 snps
  res.mr.IVW <- subset(res.mr, method == "Inverse variance weighted" & nsnp > 2)
  #Only keep Wald Ratio if only one or two SNPs
  res.mr.wald <- subset(res.mr, method == "Wald ratio")
  
  res.mr.all <- data.table::data.table()
  
  ###########################Sensitivity analyses
  for(out in unique(res.mr$outcome)){
    #Flag if the MR Egger intercept is significant
    res.mr.MREgger <- subset(res.mr, outcome==out & method == "MR Egger")
    res.mr.MREgger$pval_fdr <- p.adjust(res.mr.MREgger$egger_intercept_pvalue, method = "fdr")
    MREgger.sig <- subset(res.mr.MREgger, pval_fdr < 0.05)$cluster
    #Look at whether heterogeneity is above threshold
    Het.sig <- subset(res.mr.IVW, outcome==out & I2>0.5)$cluster
    #Flag reverse causation (Steiger filtered IVW has a different direction of effect as IVW)
    Steiger.diff <- res.mr[outcome==out & method %in% c("Steiger Inverse variance weighted","Inverse variance weighted"), length(unique(sign(beta)))>1, by=cluster] %>% .[V1==TRUE, cluster]
    
    for (z in unique(res.mr$cluster)){
      tmp.res <- res.mr[cluster==z & outcome==out]
      if(nrow(tmp.res)==0) next
      if(any(grepl(tmp.res$method, pattern = "MR_PRESSO"))){
        #Modify the distortion p-value field to have it in numeric
        subset(tmp.res, method == "MR_PRESSO_IVW_Outliers")[, Distortion_pval:=ifelse(Distortion_pval == "<0.001", 1e-4, as.numeric(Distortion_pval))]
        tmp.res[, Distortion_pval:=as.numeric(Distortion_pval)]
        # res.mr$Distortion_pval <- as.numeric(ifelse(subset(res.mr, cluster == z & method == "MR_PRESSO_IVW_Outliers")$Distortion_pval == "<0.001", 1e-4, subset(res.mr, cluster == z & method == "MR_PRESSO_IVW_Outliers")$Distortion_pval))
        beta.sensitivity <- data.frame(MRPRESSO = ifelse(subset(tmp.res, method == "MR_PRESSO_IVW_Outliers")$Distortion_pval<0.05 & !is.na(subset(tmp.res, method == "MR_PRESSO_IVW_Outliers")$Distortion_pval), subset(tmp.res, method == "MR_PRESSO_IVW_Outliers")$beta, subset(tmp.res, method == "MR_PRESSO_IVW_Raw")$beta),
                                       MREgger = ifelse(length(subset(tmp.res, method == "MR Egger")$beta)!=0, subset(tmp.res, method == "MR Egger")$beta, NA),
                                       WeightedMedian = ifelse(length(subset(tmp.res, method == "Weighted median")$beta)!=0, subset(tmp.res, method == "Weighted median")$beta, NA),
                                       Steiger = ifelse(length(subset(tmp.res, method == "Steiger Inverse variance weighted")$beta)!=0, subset(tmp.res, method == "Steiger Inverse variance weighted")$beta, NA))
      }else{
        beta.sensitivity <- data.frame(MREgger = ifelse(length(subset(tmp.res, method == "MR Egger")$beta)!=0, subset(tmp.res, method == "MR Egger")$beta, NA),
                                       WeightedMedian = ifelse(length(subset(tmp.res, method == "Weighted median")$beta)!=0, subset(tmp.res, method == "Weighted median")$beta, NA),
                                       Steiger = ifelse(length(subset(tmp.res, method == "Steiger Inverse variance weighted")$beta)!=0, subset(tmp.res, method == "Steiger Inverse variance weighted")$beta, NA))
      }
      Prop.SameDir <- sum(sign(beta.sensitivity)!=sign(subset(res.mr.IVW, cluster==z & outcome==out)$beta), na.rm=TRUE)
      
      #Flag the results that don't have same direction
      res.mr.IVW[cluster==z & outcome==out, DiffDirection:=ifelse(Prop.SameDir!=0, TRUE, FALSE)]
    }
    
    #Combine the two files
    res.mr.out <- rbind(res.mr.IVW[outcome==out], res.mr.wald[outcome==out], fill=TRUE)
    
    #Flag the results that have significant het, pleio
    res.mr.out$FlagPleiotropy <- ifelse(res.mr.out$cluster %in% unique(MREgger.sig), T, F)
    res.mr.out$FlagHeterogeneity <- ifelse(res.mr.out$cluster %in% unique(Het.sig), T, F)
    res.mr.out$FlagSteiger <- ifelse(res.mr.out$cluster %in% unique(Steiger.diff), T, F)
    
    res.mr.all <- rbind(res.mr.all, res.mr.out, fill=TRUE)
  }
  
  #Compute the final ajusted p-value
  res.mr.all$p_value_fdr <- p.adjust(res.mr.all$p_value, method = "fdr")
  res.mr.all[, p_value_plot:=ifelse(p_value_fdr<0.05 & DiffDirection==FALSE & FlagHeterogeneity==FALSE & FlagPleiotropy==FALSE, 0, 2)]
  
  return(res.mr.all)                            
}

add_odds_ratios <- function(dt){
  dt$odds_ratio <- exp(dt$beta)
  dt$ci_lower <- exp(dt$beta - 1.96 * dt$standard_error)
  dt$ci_upper <- exp(dt$beta + 1.96 * dt$standard_error)
  return(dt)
}

process.meta.results <- function(meta.dt, cur.analysis="all"){
  meta.dt <- meta.dt[analysis==cur.analysis]
  meta.dt <- meta.dt[, .(outcome, exposure, method, cluster, analysis,
                         ancestry="Meta_Analysis",
                         beta=meta.re.b,
                         standard_error=meta.re.se,
                         p_value=meta.re.p,
                         I2=meta.re.I2/100)]
  meta.dt <- add_odds_ratios(meta.dt) 
  meta.dt[, FlagHeterogeneity:=ifelse(I2>0.5, TRUE, FALSE)]
  return(meta.dt)  
}

################################################################################
#--------------------------------- VARIABLES ----------------------------------#
################################################################################
library(dplyr)
library(data.table)
pathMR="/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/"
setwd(pathMR)

################################################################################
#----------------------------------- MAIN -------------------------------------#
################################################################################
for(analysis in c("one_per_locus", "clumped", "cluster_specific")){
  dt <- data.table()
  all.files <- list.files("ResultsMR/new/")[grepl(paste0("_", analysis, ".txt"), list.files("ResultsMR/new/"))]
  dt <- data.table::rbindlist(lapply(paste0("ResultsMR/new/", all.files), data.table::fread), fill=TRUE)
  dt[, `:=`(plt.egger_intercept=NULL, plt.pval=NULL, plt.se=NULL)]
  
  #------------- Add meta-analysis results
  meta.dt <- data.table()
  for (t in c("asthma", "depression", "scz","OP", "COPD", "RA")){
    meta.dt <- rbind(meta.dt, data.table::fread(paste0("ResultsMR/new/meta_analysis_", t, "_results.txt")), fill=TRUE)
  }
  
  meta.dt <- process.meta.results(meta.dt, cur.analysis=analysis)
  
  #-------------- Forward direction
  processed.dt <- function.processMR(dt[outcome!="t2d"])
  processed.dt[, `:=`(egger_intercept_pvalue=NULL, Distortion_pval=NULL, egger_intercept=NULL, egger_intercept_standard_error=NULL)]
  unique(processed.dt)[p_value_plot==0] %>% .[, .(outcome, cluster, beta, standard_error, p_value)]
  processed.dt[p_value_fdr<0.05]
  
  final.dt <- rbind(processed.dt, meta.dt, fill=TRUE)
  final.dt$p_value_fdr <- p.adjust(final.dt$p_value, method = "fdr")
  final.dt[, y.axis.name:=paste(outcome, ancestry, sep="_")]
  
  data.table::fwrite(final.dt, paste0("ResultsMR/new/", analysis, "_results.tsv"))

  # Output grouping
  fwrite(final.dt[p_value_fdr<0.05, .(outcome, ancestry, cluster, DiffDirection, FlagPleiotropy, FlagHeterogeneity, FlagSteiger)], 
         paste0("ResultsMR/new/", analysis, "_fdr_significant_sensitivity.csv"))
}

############ Compare to main analysis and plot forest plots ###############
one.per.locus <- unique(fread("ResultsMR/new/one_per_locus_results.tsv"))
one.per.locus[, id:=paste(outcome, ancestry, cluster, sep="_")]
one.per.locus[, plot.id:=paste(outcome, ancestry, sep="_")]
one.per.locus <- one.per.locus[method=="Inverse variance weighted"]
clumped <- unique(fread("ResultsMR/new/clumped_results.tsv"))
clumped[, id:=paste(outcome, ancestry, cluster, sep="_")]
clumped[, plot.id:=paste(outcome, ancestry, sep="_")]
clumped <- clumped[method=="Inverse variance weighted"]
cluster.specific <- unique(fread("ResultsMR/new/cluster_specific_results.tsv"))
cluster.specific[, id:=paste(outcome, ancestry, cluster, sep="_")]
cluster.specific[, plot.id:=paste(outcome, ancestry, sep="_")]
cluster.specific <- cluster.specific[method=="Inverse variance weighted"]

all <- fread("ResultsMR/new/all_analysis_results_w_meta_analysis.tsv")
all[, id:=paste(outcome, ancestry, cluster, sep="_")]
all[, plot.id:=paste(outcome, ancestry, sep="_")]
all.sig <- unique(all[p_value_fdr<0.05, plot.id])

dt <- merge(all, one.per.locus[, .(id, beta.one.per.locus=beta, se.one.per.locus=standard_error, pval.fdr.one.per.locus=p_value_fdr)], by="id", all.x=TRUE)
dt <- merge(dt, clumped[, .(id, beta.clumped=beta, pval.fdr.clumped=p_value_fdr)], by="id", all.x=TRUE)
dt <- merge(dt, cluster.specific[, .(id, beta.cluster.specific=beta, pval.fdr.cluster.specific=p_value_fdr)], by="id", all.x=TRUE)
dt[, same.dir.one_per_locus:=ifelse(sign(beta)==sign(beta.one.per.locus), TRUE, FALSE), by=seq_len(nrow(dt))]
dt[, same.dir.clumped:=ifelse(sign(beta)==sign(beta.clumped), TRUE, FALSE), by=seq_len(nrow(dt))]
dt[, same.dir.cluster_specific:=ifelse(sign(beta)==sign(beta.cluster.specific), TRUE, FALSE), by=seq_len(nrow(dt))]

dt[same.dir.one_per_locus==F |same.dir.clumped==F| same.dir.cluster_specific==F]  # 82/352
dt[same.dir.one_per_locus==F | same.dir.cluster_specific==F]  # 61/352
dt[ancestry=="EUR" & (same.dir.one_per_locus==F | same.dir.cluster_specific==F)]  # 26/176
dt[same.dir.cluster_specific==F]  # 44/176
dt[ancestry=="EUR" & same.dir.cluster_specific==F]  # 20/176
dt[ancestry=="EUR" & same.dir.cluster_specific==F, .(y.axis.name, cluster, beta, beta.cluster.specific)]  # 20/176

unique(dt[p_value_fdr<0.05 & (pval.fdr.one.per.locus>0.05 | pval.fdr.cluster.specific>0.05), 
   .(id, p_value_fdr, pval.fdr.one.per.locus, pval.fdr.cluster.specific)], by=c("id", "p_value_fdr"))  # pval.fdr.clumped, 

unique(dt[p_value_fdr>0.05 & (pval.fdr.one.per.locus<0.05 | pval.fdr.cluster.specific<0.05), 
   .(id, p_value_fdr, pval.fdr.one.per.locus, pval.fdr.cluster.specific)])  # pval.fdr.clumped, 

plot.dt <- dt[id %in% all.sig]

plot.dt <- rbind(all[plot.id %in% all.sig, .(analysis="all", beta, standard_error, p_value, p_value_fdr, id, outcome, cluster, ancestry, plot.id)],
                 one.per.locus[plot.id %in% all.sig, .(analysis="one.per.locus", beta, standard_error, p_value, p_value_fdr, id, outcome, cluster, ancestry, plot.id)],
                 cluster.specific[plot.id %in% all.sig, .(analysis="cluster.specific", beta, standard_error, p_value, p_value_fdr, id, outcome, cluster, ancestry, plot.id)])

info.datasets <- as.data.table(readxl::read_excel("info_datasets.xlsx"))
plot.dt <- merge(plot.dt, unique(info.datasets[, .(outcome=TRAIT, disease, CLASS)]), by="outcome", all.x=TRUE)
plot.dt[cluster=="Full", cluster:="All"]
plot.dt[ancestry=="Meta_Analysis", ancestry:="Meta analysis"]
plot.dt[disease=="Anorexia", disease:="Anorexia nervosa"]
plot.dt[, analysis:=ifelse(analysis=="one.per.locus", "One variant per locus", ifelse(analysis=="cluster.specific", "Cluster-specific clumping", "All T2D index variants"))]
plot.dt <- add_odds_ratios(plot.dt)
plot.dt[, twoFC_beta:=log(2)*beta]
plot.dt[, twoFC_standard_error:=log(2)*standard_error]

library(ggplot2)
for(i in all.sig){
  ggforestplot::forestplot(
    df = plot.dt[plot.id==i],
    name = cluster,
    estimate = twoFC_beta,
    se = twoFC_standard_error,
    colour = analysis,
    pvalue = p_value_fdr,
    psignif = 0.05,
    title = paste(unique(plot.dt[plot.id==i, disease]), paste0("(", unique(plot.dt[plot.id==i, ancestry]), ")"), sep="\n"),
    xlab = "OR (95% CI) per 2-fold increase in genetically \ndetermined T2D risk",
    y.lab = NA, 
    logodds = TRUE) +
    ggplot2::scale_x_continuous(breaks = seq(0, max(plot.dt$ci_upper), by = 0.1)) +
    theme(axis.text.y = element_text(size = 14), plot.title = element_text(hjust = 0.5)) +
    labs(color='Approach to select T2D IVs')
  
  ggplot2::ggsave(paste0("paper/SupFigs/", i, "_ivs_selection.png"), dpi=600, height=6, width=7.5)
}
