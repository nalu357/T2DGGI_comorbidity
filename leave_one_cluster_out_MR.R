################################################################################
#--------------------------------- FUNCTIONS ----------------------------------#
################################################################################
inverse.ivs.MR <- function(datasets.dt, path.code, ref.bfile){
  one.per.locus <- data.table::fread("one_ivs_per_locus.txt")
  
  for(pop in unique(datasets.dt$ANCESTRY)){
    ref.bfile.path <- paste0(ref.bfile, pop)
    clumped <- data.table::fread(paste0("t2d_ivs/index_snps_clumped_t2d_", pop, ".txt"))
    
    for(t in datasets.dt[ANCESTRY==pop, TRAIT]){
      trait <- datasets.dt[ANCESTRY==pop & TRAIT==t]
      
      # Read T2D ivs
      t2d.all.ivs <- data.table::fread(paste0(trait$PATH, "t2d_", pop, "_ivs.txt"))
      if("clumped.iv" %in% colnames(t2d.all.ivs)) t2d.all.ivs[, final.ivs:=ifelse(clumped.iv=="", clumped.ivs, clumped.iv), by=seq_len(nrow(t2d.all.ivs))]
      else t2d.all.ivs[, final.ivs:=clumped.ivs]

      t2d.ivs <- t2d.all.ivs$final.ivs
      
      cluster.mr <- data.table::data.table()
      for (c in unique(t2d.all.ivs$`Cluster assignment`)){
        if(c=="Liver/lipid metabolism") next  # do not run for the cluster that has only 2 SNPs
        
        ivs.cluster <- t2d.all.ivs[`Cluster assignment`!=c, final.ivs]
        
        cluster.mr <- rbind(cluster.mr, wrap.MR(exposure.file=paste0(trait$PATH, "t2d_", pop, "_forMR.txt"),
                                                outcome.file=paste0(trait$PATH, trait$TRAIT, "_", pop, "_forMR.txt"),
                                                exposure.name="t2d",
                                                outcome.name=trait$TRAIT,
                                                ivs=ivs.cluster,
                                                cluster.name=c,
                                                ref.bfile.path,
                                                correlated=TRUE), fill=TRUE)
      }
      cluster.mr$ancestry <- pop
      cluster.mr$cluster <- paste0("without ", cluster.mr$cluster)
      if(nrow(cluster.mr)==0) return(cluster.mr)
      
      data.table::setnames(cluster.mr, old=c("b", "se", "pval", "or", "or_lci95", "or_uci95"),
                           new=c("beta", "standard_error", "p_value", "odds_ratio", "ci_lower", "ci_upper"))
      
      data.table::fwrite(cluster.mr[, `:=` (ancestry=pop)], paste0("ResultsMR/new/inverse_ivs_mr_results_", trait$TRAIT, "_", pop, ".txt"), na = "NA", sep = "\t", quote = F)
    }
  }
}

meta.MR <- function(datasets.dt, trait){
  meta.dt <- data.table::data.table()
  for (pop in unique(datasets.dt$ANCESTRY)){
    for (a in analysis){
      res.dt <- data.table::fread(paste0("ResultsMR/new/inverse_ivs_mr_results_", trait, "_", pop, ".txt"))
      meta.dt <- rbind(meta.dt, 
                       res.dt[method %in% c("Wald ratio", "Inverse variance weighted"), .(outcome, exposure, cluster, ancestry, cluster, analysis, method, nsnp, beta, standard_error, p_value)], 
                       fill=TRUE)
    }
  }
  
  res.meta <- data.table::data.table()
  for(c in unique(meta.dt$cluster)){
    for (a in analysis){
      dt <- meta.dt[exposure!=trait & cluster==c & analysis==a]
      # Run meta-analysis
      re.meta <- metafor::rma(yi=beta, sei=standard_error, data=dt, method = "REML", test="z")
      re.meta.knha <- metafor::rma(yi=beta, sei=standard_error, data=dt, method = "REML", test="knha")
      fe.meta <- metafor::rma(yi=beta, sei=standard_error, data=dt, method = "FE")

      # Convert data table from long to wide
      dt <- dt[, exposure:=sub("__[^__]+$", "", exposure)]
      dt <- dt[, .SD, .SDcols = unique(names(dt))]
      dt <- data.table::dcast(dt, outcome + analysis + cluster + exposure + method ~ ancestry, value.var=c("nsnp", "beta", "standard_error", "p_value"), sep=".")

      # Add meta-analysis results
      dt$nstudies <- re.meta$k
      dt$meta.re.b <- re.meta$beta[,1]
      dt$meta.re.se <- re.meta$se
      dt$meta.re.p <- re.meta$pval
      dt$meta.re.I2 <- re.meta$I2
      dt$meta.re.QE <- re.meta$QE
      dt$meta.re.QEp <- re.meta$QEp
      dt$meta.re.knha.b <- re.meta.knha$beta[,1]
      dt$meta.re.knha.se <- re.meta.knha$se
      dt$meta.re.knha.p <- re.meta.knha$pval
      dt$meta.re.knha.I2 <- re.meta.knha$I2
      dt$meta.re.knha.QE <- re.meta.knha$QE
      dt$meta.re.knha.QEp <- re.meta.knha$QEp
      dt$meta.fe.b <- fe.meta$beta[,1]
      dt$meta.fe.se <- fe.meta$se
      dt$meta.fe.p <- fe.meta$pval
      dt$meta.fe.I2 <- fe.meta$I2
      dt$meta.fe.QE <- fe.meta$QE
      dt$meta.fe.QEp <- fe.meta$QEp

      res.meta <- rbind(res.meta, dt, fill=TRUE)
    }
  }
  
  data.table::fwrite(res.meta, paste0("ResultsMR/new/inverse_ivs_meta_analysis_", trait, "_results.txt"))
  return(res.meta)
}

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

process.meta.results <- function(meta.dt){
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
.libPaths(c("/home/itg/ana.arruda/R/x86_64-pc-linux-gnu-library/4.1",
            "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/Rpackages"))
pathMR="/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/"
setwd(pathMR)
path_T2DGGI = "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/datasets/T2DGGI_bothbuilds/"
path_1kG = "/lustre/groups/itg/shared/referenceData/1kG/EAS-SAS-AMR-EUR-AFR_1kg_v3/"
path_code = paste0(pathMR, "scripts/")
tmp_dir = "/lustre/groups/itg/teams/zeggini/projects/T2D-diamante/T2DGGI_MR/project2/tmpdata/"

##Dataframe with the different datasets to use
info.datasets <- data.table::as.data.table(readxl::read_excel("info_datasets.xlsx"))
source(paste0(path_code, "MR_functions.R"))

#Number of cases and controls in each of the ancestries (numbers obtained from Kostas)
info.T2DGGI <- data.frame(ancestry = c("AFR", "EAS", "EUR", "AMR", "SAS"),
                          n.cases = c(50251, 88109, 242283, 29375, 16832),
                          n.controls = c(103909, 339395, 1569734, 59368, 33767))

################################################################################
#----------------------------------- MAIN -------------------------------------#
################################################################################
pop="EUR"
api.token <- "d1fe6afcb763"

# Run MR
inverse.ivs.MR(datasets.dt=info.datasets,
               path.code=path_code,
               ref.bfile=path_1kG)

# Run meta-analysis
for(t in unique(info.datasets[, if(.N>1) .SD, by=TRAIT]$TRAIT)){
  meta.MR(datasets.dt=info.datasets[TRAIT==t],
          trait=t)
}

# Process output
library(dplyr)
library(data.table)

dt <- data.table()
all.files <- list.files("ResultsMR/new/")[grepl(paste0("inverse_ivs_mr_"), list.files("ResultsMR/new/"))]
dt <- data.table::rbindlist(lapply(paste0("ResultsMR/new/", all.files), data.table::fread), fill=TRUE)
dt[, `:=`(plt.egger_intercept=NULL, plt.pval=NULL, plt.se=NULL)]
dt <- dt[outcome!="alzheimer"]
fwrite(dt, "ResultsMR/new/inverse_ivs_all_analysis_results.txt")

#-------------- T2D --> disease direction
processed.dt <- function.processMR(dt[outcome!="t2d"])
processed.dt[, .(outcome, cluster, beta, standard_error, p_value, p_value_fdr)]
processed.dt[, `:=`(egger_intercept_pval=NULL, egger_intercept=NULL, egger_intercept_se=NULL)]
unique(processed.dt)[p_value_plot==0] %>% .[, .(outcome, cluster, beta, standard_error, p_value)]
processed.dt[p_value_fdr<0.05, .(outcome, ancestry, cluster, beta, p_value)]

#------------- Add meta-analysis results
meta.dt <- data.table()
for (t in c("asthma", "OP", "COPD", "RA", "depression", "scz")){
  meta.dt <- rbind(meta.dt, data.table::fread(paste0("ResultsMR/new/inverse_ivs_meta_analysis_", t, "_results.txt")), fill=TRUE)
}

meta.dt <- process.meta.results(meta.dt)

final.dt <- rbind(processed.dt, meta.dt, fill=TRUE)
final.dt$p_value_fdr <- p.adjust(final.dt$p_value, method = "fdr")
final.dt[, .(outcome, ancestry, cluster, beta, p_value, p_value_fdr)]
final.dt[p_value_fdr<0.05, .(outcome, ancestry, cluster, beta, p_value)]

final.dt[ancestry!="Meta_Analysis", p_value_plot:=ifelse(p_value_fdr<0.05 & DiffDirection==FALSE & FlagHeterogeneity==FALSE & FlagPleiotropy==FALSE, 0, 2)]
final.dt[ancestry=="Meta_Analysis", p_value_plot:=ifelse(p_value_fdr<0.05 & FlagHeterogeneity==FALSE, 0, 2)]
final.dt[ancestry=="Meta_Analysis" & p_value_fdr<0.05]
final.dt[, y.axis.name:=paste(outcome, ancestry, sep="_")]
final.dt[p_value_fdr<0.05 & DiffDirection==TRUE]  # no consistent direction of effect across methods
final.dt[p_value_fdr<0.05 & FlagPleiotropy==TRUE]  # evidence of pleiotropy
final.dt[p_value_fdr<0.05 & FlagHeterogeneity==TRUE] # evidence of heterogeneity
final.dt[p_value_fdr<0.05 & FlagSteiger==TRUE] # evidence of reverse causation
final.dt[p_value_plot==0] # most robust results

data.table::fwrite(final.dt, "ResultsMR/new/inverse_ivs_all_analysis_results_w_meta_analysis.tsv")