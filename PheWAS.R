library(dplyr)
library(data.table)
pathMR="C:/Users/AnaLuizaArruda/OneDrive - Helmholtz Zentrum München/Projects/T2DGGI/T2DGGI_MR/project_2/"
setwd(pathMR)

mr.results <- fread("ResultsMR/new/all_analysis_results_w_meta_analysis.tsv")
sig.mr.results <- mr.results[p_value_fdr<0.05]
clusters <- c("Obesity", "Beta cell -PI", "Liver/lipid metabolism", "Lipodystrophy", "Beta cell +PI",
              "Metabolic syndrome", "Residual glycaemic", "Body fat")  

res.dt <- data.table::data.table()
all.dt <- data.table::data.table()
for(pop in c("eur", "eas", "amr", "afr", "sas")){
  full <- fread(paste0("PheWAS/stg015/", pop, "_cohort/grs__overall.csv"))
  full <- full[!(phecode_category %in% c("Cardiovascular", "Endocrine/Metab", "Gastrointestinal", "Infections", "Genitourinary", "Pregnancy", "Blood/Immune", "Dermatological", "Neoplasms", "Genetic", "Symptoms"))]
  all.dt <- rbind(all.dt, full[,`:=`(ancestry=pop, cluster="Full")])
  sig.full <- full[p_value<0.05/nrow(full)]
  out.dt <- sig.full[, .(phecode, cases, controls, phecode_string, phecode_category, beta, p_value, standard_error)][order(phecode_category)]
  if(nrow(out.dt)==0) next
  data.table::fwrite(out.dt, paste0("PheWAS/significant_overall_results_", pop, ".csv"))
  res.dt <- rbind(res.dt, out.dt[,`:=`(ancestry=pop, cluster="Full")])
  
  for(c in seq(1:8)){
    if(c==3) next
    cur.cluster <- clusters[c]
    
    phewas <- fread(paste0("PheWAS/stg015/", pop, "_cohort/grs__cluster", c, ".csv"))
    phewas <- phewas[!(phecode_category %in% c("Cardiovascular", "Endocrine/Metab", "Gastrointestinal", "Infections", "Genitourinary", "Pregnancy", "Blood/Immune", "Dermatological", "Neoplasms", "Genetic", "Symptoms"))]
    sig.phewas <- phewas
    sig.phewas[!(phecode_category %in% c("Cardiovascular", "Endocrine/Metab", "Gastrointestinal", "Infections", "Genitourinary", "Pregnancy", "Blood/Immune", "Dermatological", "Neoplasms", "Genetic", "Symptoms"))]
    out.dt <- sig.phewas[, .(phecode, cases, controls, phecode_string, phecode_category, beta, p_value, standard_error)][order(phecode_category)]
    if(nrow(out.dt)==0) next
    res.dt <- rbind(res.dt, out.dt[,`:=`(ancestry=pop, cluster=cur.cluster)])
    sig.phewas[phecode_category=="Muscloskeletal", .(phecode_string, beta, p_value)]
    sig.phewas[phecode_category=="Mental", .(phecode_string, beta, p_value)]
  }
}

cluster.dt <- data.table::data.table(cluster.name=clusters, cluster.number=seq(1:8))

### Meta-analysis
meta <- fread("PheWAS/meta_analysis/meta_analysis.tsv")
meta[, cluster.new:=cluster.dt[cluster.number==cluster, cluster.name], by=seq_len(nrow(meta))]
meta[, cluster:=cluster.new]
meta[, cluster.new:=NULL]
meta <- rbind(meta, fread("PheWAS/meta_analysis/meta_analysis_overall.tsv"))
meta <- meta[!(phecode_category %in% c("Cardiovascular", "Endocrine/Metab", "Gastrointestinal", "Infections", "Genitourinary", "Pregnancy", "Blood/Immune", "Dermatological", "Neoplasms", "Genetic", "Symptoms"))]
sig.meta <- meta  
out.dt <- sig.meta[, .(cluster, phecode, phecode_string, phecode_category, meta_beta, meta_pvalue, meta_se)][order(phecode_category)]
res.dt <- rbind(res.dt, out.dt[,.(cluster, phecode, phecode_string, phecode_category, ancestry="meta_analysis", beta=meta_beta, p_value=meta_pvalue, standard_error=meta_se)], fill=TRUE)
res.dt[, p_value_fdr:=p.adjust(p_value, method="fdr"), by="cluster"]
data.table::fwrite(res.dt, "PheWAS/all_results.csv")
data.table::fwrite(res.dt[p_value_fdr<0.05], "PheWAS/all_significant_results_meta_analysis.csv")


my.diseases <- c("Osteoarthritis", "Osteoporosis", "Carpal tunnel syndrome", "Chronic obstructive pulmonary disease", "Asthma", "Bipolar disorder",
                 "Major depressive disorder", "Back pain", "Cataract", "Mood")
similar.phenos <- c("joint", "Mononeuropathies", "Mononeuropathy", "Dyspnea", "Respiratory failure", "Chronic respiratory failure", "Pain in knee", "Pain in joint", "Fractures", "Polyneuropathy", "Polyneuropathies")
res.dt[p_value_fdr<0.05, .N, by=c("cluster", "ancestry")]
res.dt[grep(paste(c(my.diseases, similar.phenos), collapse="|"), phecode_string)][order(phecode_string)]
res.dt[grep(paste(c(my.diseases, similar.phenos), collapse="|"), phecode_string)][order(phecode_string) & cluster=="Full"]
res.dt[grep(paste(c(my.diseases, similar.phenos), collapse="|"), phecode_string)][order(phecode_string) & cluster=="Metabolic syndrome"]
res.dt[grep(paste(c(my.diseases, similar.phenos), collapse="|"), phecode_string)][order(phecode_string) & cluster=="Obesity"]