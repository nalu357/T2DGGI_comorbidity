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
    for (pop in unique(res.mr[outcome==out, ancestry])){
      #Flag if the MR Egger intercept is significant
      res.mr.MREgger <- subset(res.mr, outcome==out & method == "MR Egger" & ancestry==pop)
      res.mr.MREgger$pval_fdr <- p.adjust(res.mr.MREgger$egger_intercept_pvalue, method = "fdr")
      MREgger.sig <- subset(res.mr.MREgger, pval_fdr < 0.05)$cluster
      #Look at whether heterogeneity is above threshold
      Het.sig <- subset(res.mr.IVW, outcome==out & I2>0.5 & ancestry==pop)$cluster
      #Flag reverse causation (Steiger filtered IVW has a different direction of effect as IVW)
      Steiger.diff <- res.mr[outcome==out & ancestry==pop & method %in% c("Steiger Inverse variance weighted","Inverse variance weighted"), length(unique(sign(beta)))>1, by=cluster] %>% .[V1==TRUE, cluster]
      
      for (z in unique(res.mr$cluster)){
        tmp.res <- res.mr[cluster==z & outcome==out & ancestry==pop]
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
        Prop.SameDir <- sum(sign(beta.sensitivity)!=sign(subset(res.mr.IVW, cluster==z & outcome==out & ancestry==pop)$beta), na.rm=TRUE)
        
        #Flag the results that don't have same direction
        res.mr.IVW[cluster==z & outcome==out & ancestry==pop, DiffDirection:=ifelse(Prop.SameDir!=0, TRUE, FALSE)]
      }
      
      #Combine the two files
      res.mr.out <- rbind(res.mr.IVW[outcome==out & ancestry==pop], res.mr.wald[outcome==out & ancestry==pop], fill=TRUE)
      
      #Flag the results that have significant het, pleio
      res.mr.out$FlagPleiotropy <- ifelse(res.mr.out$cluster %in% unique(MREgger.sig), T, F)
      res.mr.out$FlagHeterogeneity <- ifelse(res.mr.out$cluster %in% unique(Het.sig), T, F)
      res.mr.out$FlagSteiger <- ifelse(res.mr.out$cluster %in% unique(Steiger.diff), T, F)
      
      res.mr.all <- rbind(res.mr.all, res.mr.out, fill=TRUE)
    }
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
analysis="all"
dt <- data.table()
all.files <- list.files("ResultsMR/new/")[grepl(paste0("_", analysis, ".txt"), list.files("ResultsMR/new/"))]
dt <- data.table::rbindlist(lapply(paste0("ResultsMR/new/", all.files), data.table::fread), fill=TRUE)
dt[, `:=`(plt.egger_intercept=NULL, plt.pval=NULL, plt.se=NULL)]
dt <- dt[outcome!="alzheimer" & exposure!="alzheimer"]
fwrite(dt, "ResultsMR/new/all_analysis_results.txt")

#################################################
#-------------- Forward direction ---------------
#################################################
#-------------- T2D --> disease direction
processed.dt <- function.processMR(dt[outcome!="t2d"])
processed.dt[, `:=`(egger_intercept_pvalue=NULL, Distortion_pval=NULL, egger_intercept=NULL, egger_intercept_standard_error=NULL)]
unique(processed.dt)[p_value_plot==0] %>% .[, .(outcome, cluster, beta, standard_error, p_value)]
processed.dt[p_value_fdr<0.05, .(outcome, ancestry, cluster, beta, p_value)]

#------------- Add meta-analysis results
meta.dt <- data.table()
for (t in c("asthma", "OP", "COPD", "RA", "depression", "scz")){
  meta.dt <- rbind(meta.dt, data.table::fread(paste0("ResultsMR/new/meta_analysis_", t, "_results.txt")), fill=TRUE)
}

meta.dt <- process.meta.results(meta.dt)

#------------ Add reverse direction to adjust p-value
reverse.dt <- data.table::fread("ResultsMR/new/reverse_results.tsv")
final.dt <- data.table::fread("ResultsMR/new/all_analysis_results_w_meta_analysis.tsv")
final.dt <- rbind(final.dt, reverse.dt, fill=TRUE)
final.dt$p_value_fdr <- p.adjust(final.dt$p_value, method = "fdr")
final.dt[p_value_fdr<0.05, .(outcome, ancestry, cluster, beta, p_value)]
data.table::fwrite(final.dt, "ResultsMR/new/all_analysis_results_w_meta_analysis_reverse.tsv")

final.dt[ancestry!="Meta_Analysis", p_value_plot:=ifelse(p_value_fdr<0.05 & DiffDirection==FALSE & FlagHeterogeneity==FALSE & FlagPleiotropy==FALSE, 0, 2)]
final.dt[ancestry=="Meta_Analysis", p_value_plot:=ifelse(p_value_fdr<0.05 & FlagHeterogeneity==FALSE, 0, 2)]
final.dt[ancestry=="Meta_Analysis" & p_value_fdr<0.05]
final.dt[, y.axis.name:=paste(outcome, ancestry, sep="_")]
final.dt[p_value_fdr<0.05 & DiffDirection==TRUE, .(outcome, ancestry, cluster, beta, p_value)]  # no consistent direction of effect across methods
final.dt[p_value_fdr<0.05 & FlagPleiotropy==TRUE, .(outcome, ancestry, cluster, beta, p_value)]  # evidence of pleiotropy
final.dt[p_value_fdr<0.05 & FlagHeterogeneity==TRUE, .(outcome, ancestry, cluster, beta, p_value)] # evidence of heterogeneity
final.dt[p_value_fdr<0.05 & FlagSteiger==TRUE, .(outcome, ancestry, cluster, beta, p_value)] # evidence of reverse causation
final.dt[p_value_plot==0] # most robust results

# Plot FDR significant results
ggforestplot::forestplot(
  df = final.dt[p_value_fdr<0.05],
  name = y.axis.name,
  estimate = beta,
  se = standard_error,
  pvalue = p_value_plot,
  psignif = 1,
  colour = cluster,
  xlab = "Odds ratio (95% CI)",
  y.lab = NA, 
  title = "All FDR significant results",
  logodds = TRUE)
ggplot2::ggsave("ResultsMR/new/plots/all_fdr_significant.png", dpi=600, height=11, width=7)

# Output grouping
fwrite(processed.dt[p_value_fdr<0.05, .(outcome, ancestry, cluster, DiffDirection, FlagPleiotropy, FlagHeterogeneity, FlagSteiger)], "ResultsMR/new/all_fdr_significant_sensitivity.csv")

# Ancestry-specific results plot
result <- fread("ResultsMR/new/all_analysis_results_w_meta_analysis.tsv")
plot.dt <- result[p_value_fdr<0.05 & ancestry!="Meta_Analysis" & outcome %in% c("OP", "depression", "asthma", "scz", "RA", "COPD")]
ggforestplot::forestplot(
  df = plot.dt,
  name = y.axis.name,
  estimate = beta,
  se = standard_error,
  pvalue = p_value_plot,
  psignif = 1,
  colour = cluster,
  xlab = "Odds ratio (95% CI)",
  y.lab = NA, 
  title = "Ancestry-specific FDR significant results",
  logodds = TRUE)
ggplot2::ggsave("ResultsMR/new/plots/ancestry_specific_fdr_significant.png", dpi=600, height=8, width=7)

data.table::fwrite(final.dt, "ResultsMR/new/all_analysis_results_w_meta_analysis.tsv")

#################################################
#-------------- Reverse direction ---------------
#################################################
reverse.dt <- function.processMR(dt[outcome=="t2d"])
reverse.dt[p_value_fdr<0.05]
reverse.dt[p_value_fdr<0.05 & DiffDirection==TRUE]
reverse.dt[p_value_fdr<0.05 & DiffDirection==TRUE & FlagPleiotropy==TRUE]  # evidence of pleiotropy
reverse.dt[p_value_fdr<0.05 & DiffDirection==TRUE & FlagHeterogeneity==TRUE] # evidence of heterogeneity
reverse.dt[p_value_fdr<0.05 & DiffDirection==TRUE & FlagSteiger==TRUE] # evidence of reverse causation
reverse.dt[, y.axis.name:=paste(exposure, ancestry, sep="_")]
data.table::fwrite(reverse.dt, "ResultsMR/new/reverse_results.tsv")
fwrite(reverse.dt[p_value_fdr<0.05, .(outcome, ancestry, cluster, DiffDirection, FlagPleiotropy, FlagHeterogeneity, FlagSteiger)], "ResultsMR/new/reverse_all_fdr_significant_sensitivity.csv")

#------------- Add meta-analysis results
reverse.meta.dt <- data.table()
for (t in c("asthma", "OP", "COPD", "RA", "depression", "scz")){
  reverse.meta.dt <- rbind(reverse.meta.dt, data.table::fread(paste0("ResultsMR/new/reverse_meta_analysis_", t, "_results.txt")), fill=TRUE)
}

reverse.meta.dt <- process.meta.results(reverse.meta.dt)

reverse.final.dt <- rbind(reverse.dt, reverse.meta.dt, fill=TRUE)
reverse.final.dt$p_value_fdr <- p.adjust(reverse.final.dt$p_value, method = "fdr")
reverse.final.dt[p_value_fdr<0.05, .(exposure, outcome, ancestry, cluster, beta, p_value)]

reverse.final.dt[ancestry!="Meta_Analysis", p_value_plot:=ifelse(p_value_fdr<0.05 & DiffDirection==FALSE & FlagHeterogeneity==FALSE & FlagPleiotropy==FALSE, 0, 2)]
reverse.final.dt[ancestry=="Meta_Analysis", p_value_plot:=ifelse(p_value_fdr<0.05 & FlagHeterogeneity==FALSE, 0, 2)]
reverse.final.dt[ancestry=="Meta_Analysis" & p_value_fdr<0.05]
reverse.final.dt[, y.axis.name:=paste(outcome, ancestry, sep="_")]
reverse.final.dt[p_value_fdr<0.05 & DiffDirection==TRUE]  # no consistent direction of effect across methods
reverse.final.dt[p_value_fdr<0.05 & FlagPleiotropy==TRUE]  # evidence of pleiotropy
reverse.final.dt[p_value_fdr<0.05 & FlagHeterogeneity==TRUE] # evidence of heterogeneity
reverse.final.dt[p_value_fdr<0.05 & FlagSteiger==TRUE] # evidence of reverse causation
reverse.final.dt[p_value_plot==0] # most robust results
data.table::fwrite(reverse.final.dt, "ResultsMR/new/reverse_analysis_results_w_meta_analysis.tsv")
