################################################################################
#--------------------------------- FUNCTIONS ----------------------------------#
################################################################################
get_value <- function(x, data){
  if (is.na(x)) return(NA)
  else return(data[, get(x)])
}

get_trait_data <- function(trait, tmp_data){
  dt <- data.table::fread(trait$FILE, tmpdir=tmp_data)
  output.dt <- dt[, .(p_value=get(trait$pval.col), 
                      beta=get(trait$beta.col),
                      standard_error=get(trait$se.col), 
                      rs_id=get(trait$rsid.col), 
                      effect_allele=toupper(get(trait$ea.col)), 
                      other_allele=toupper(get(trait$nea.col)))]
  output.dt$chr <- sapply(trait$chr.col, get_value, data=dt)
  output.dt$pos <- sapply(trait$pos.col, get_value, data=dt)
  output.dt$effect_allele_frequency <- sapply(trait$eaf.col, get_value, data=dt)
  output.dt$n <- sapply(trait$n.col, get_value, data=dt)
  output.dt$ncases <- sapply(trait$ncases.col, get_value, data=dt)

  if(!is.na(output.dt$chr)) output.dt[, chr:=as.integer(sub("chr", "", chr))]
  
  #Keep only SNVs with rsID
  output.dt <- output.dt[effect_allele %in% c("A", "T", "C", "G")]
  output.dt <- output.dt[other_allele %in% c("A", "T", "C", "G")]
  output.dt<- subset(output.dt, !is.na(rs_id))
  output.dt<- subset(output.dt, !grepl(",", rs_id))
  output.dt<- output.dt[rs_id!="."]
  output.dt<- output.dt[rs_id!=""]
  
  return(output.dt)
}

get_t2d_data <- function(info.t2d, path_T2DGGI, ancestry_T2DGGI){
  # Read ancestry-specific T2D data
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

clumping <- function(out.path, exposure.name, exposure.file, ref.bfile.path, pval.thr, r2_thr, kb_thr, plink, rsid.col="rs_id", pval.col="p_value"){
  #Clump the variants
  if(!file.exists(paste0(out.path, exposure.name, ".clumped"))) system(paste0(plink, " --bfile ", ref.bfile.path, " --clump-p1 ", pval.thr, " --clump-r2 ", r2_thr, " --clump-kb ", kb_thr, " --clump ", exposure.file, " --clump-snp-field ", rsid.col, " --clump-field ", pval.col, " --out ", out.path, exposure.name))
  
  #Import the clumping results
  clumping.res <- data.table::fread(paste0(out.path, exposure.name, ".clumped"))
  
  # clumped.ivs <- clumping.res$SNP
  return(clumping.res)
}

proxy_search <- function(need.proxy, outcome, clumping.res, ref.bfile.path){
  proxies <- NULL
  for (iv in need.proxy){
    cluster.snps <- unlist(strsplit(gsub("\\(1\\)", "", clumping.res[SNP==iv, SP2]), split=","))
    
    # Check LD between iv and potential.proxy
    ld.matrix <- ieugwasr::ld_matrix(c(iv, cluster.snps),
                                     plink_bin = genetics.binaRies::get_plink_binary(),
                                     bfile = ref.bfile.path)
    
    #Take the corresponding iv and remove the correlation with itself
    ld.vector <- ld.matrix[grep(rownames(ld.matrix), pattern = iv),-grep(rownames(ld.matrix), pattern = iv)]
    rm(ld.matrix)
    
    #Check if any of the proxies is in t2d
    all.potential.proxies = sub("_.*", "", names(which((ld.vector**2)>0.8)))
    if(nrow(subset(outcome[rs_id %in% all.potential.proxies]))>0){
      for (potential.proxy in sub("_.*", "", names(which((ld.vector**2)>0.8)))){
        # Check if potential.proxy is in T2D data
        if(nrow(outcome[rs_id==potential.proxy])!=0) break
      }
      proxies <- c(proxies, potential.proxy)
    }
  }
  return(proxies)
}

get_ivs <- function(exposure, exposure.name, exposure.file, outcome, out.path, ref.bfile.path, api.token, rsid.col, pval.col, ancestry){
  #Clumping
  clumping.res <- tryCatch(
    clumping(out.path, exposure.name, exposure.file, ref.bfile.path, pval.thr=5e-8, r2_thr = 0.001, kb_thr = 10000, plink = "plink", rsid.col, pval.col),
    error=function(e) e
  )
  if(inherits(clumping.res, "error")){
    print("No significant clumping results")
    file.create(paste0(out.path, exposure.name, "_", ancestry, "_ivs.txt"))
    return(NULL)
  }
  clumped.ivs <- clumping.res$SNP
  
  #F-stat
  exposure[, snp_fstat:=(exposure$beta)^2/(exposure$standard_error)^2]
  
  #Only keep the variants that are clumped
  exposure[, iv:=ifelse(rs_id %in% clumped.ivs & snp_fstat>=10, TRUE, FALSE)]
  
  #Check which clumped IVs are not available in T2D of matching ancestry
  need.proxy <- subset(exposure, iv)$rs_id[which(!(subset(exposure, iv)$rs_id %in% outcome$rs_id))]
  
  #Find proxies
  if(length(need.proxy)>0){
    proxies <- proxy_search(need.proxy, outcome, clumping.res, ref.bfile.path)
    ## Update IVs in sig.qtl
    exposure[rs_id %in% need.proxy, iv:=FALSE]
    exposure[rs_id %in% proxies, iv:=TRUE]
  }
  #Save list of IVs
  
  data.table::fwrite(data.table::data.table(clumped.ivs=exposure[iv==TRUE, rs_id]), paste0(out.path, exposure.name, "_", ancestry, "_ivs.txt"))
  
  return(exposure[iv==TRUE, rs_id])
}

proxy_lookup <- function(snp, api.token, population, outcome){
  proxies <- NULL
  all.potential.proxies <- data.table::as.data.table(LDlinkR::LDproxy(snp = snp, 
                                                                      pop = population, 
                                                                      r2d = "r2", 
                                                                      token = api.token,
                                                                      genome_build = "grch37"))
  all.potential.proxies <- all.potential.proxies[R2>0.8, RS_Number]
  chosen.proxy <- NA
  
  #Check if any of the proxies is included in the oa data
  if(nrow(outcome[rs_id %in% all.potential.proxies])>0){
    for (potential.proxy in all.potential.proxies){
      # Check if potential.proxy is in T2D data
      if(nrow(outcome[rs_id==potential.proxy])!=0) break
    }
    chosen.proxy <- potential.proxy
  }
  
  return(chosen.proxy)
}

get_t2d_ivs <- function(exposure, exposure.name, outcome, out.path, t2d.ivs.file, api.token, ancestry){
  t2d.ivs.dt <- data.table::fread(t2d.ivs.file)
  
  #F-stat
  exposure[, snp_fstat:=(exposure$beta)^2/(exposure$standard_error)^2]

  #Only keep the variants that are in the t2d.ivs list
  exposure[, iv:=ifelse(rs_id %in% t2d.ivs.dt$clumped.ivs & snp_fstat>=10, TRUE, FALSE)]
  
  #Check which clumped IVs are not available in T2D of matching ancestry
  need.proxy <- subset(exposure, iv)$rs_id[which(!(subset(exposure, iv)$rs_id %in% outcome$rs_id))]
  
  #Find proxies
  if(length(need.proxy)>0){
    for(i in need.proxy){
      chosen.proxy <- tryCatch(
        proxy_lookup(snp=i, api.token, population=ancestry, outcome),
        error=function(e) e
      )
      
      if(inherits(chosen.proxy, "error")) return(next)
      
      # create data table with columns: iv, OA trait, proxy
      if(!is.na(chosen.proxy)){
        exposure[rs_id==i, iv:=FALSE]
        exposure[rs_id==chosen.proxy, iv:=TRUE]
        t2d.ivs.dt[clumped.ivs==i, clumped.iv:=chosen.proxy]
      }
    }
  }
  
  #Save list of IVs
  data.table::fwrite(t2d.ivs.dt, paste0(out.path, "t2d_", ancestry, "_ivs.txt"))
  
  return(exposure[iv==TRUE, rs_id])
}

preprocess_t2d_ivs <- function(index.snps.file, one.per.locus.file){
  index.snps <- data.table::as.data.table(readxl::read_excel(index.snps.file))
  one_per_locus <- data.table::as.data.table(index.snps[, .SD[which.min(MR_MEGA_P)], by=Locus]$`Index SNV`)
  data.table::fwrite(one_per_locus, one.per.locus.file)
  return(one_per_locus)
}
  
preprocess.data <- function(datasets.dt, project.path, tmp_data, info.t2d, t2d.path, path.ref, api.token="d1fe6afcb763"){
  for(pop in unique(datasets.dt$ANCESTRY)){
    t2d <- get_t2d_data(info.t2d, t2d.path, pop)
    
    ref.bfile.path=paste0(path.ref, pop)
    
    for(t in datasets.dt[ANCESTRY==pop, TRAIT]){
      trait <- datasets.dt[ANCESTRY==pop & TRAIT==t]
      
      system(paste0("mkdir -p ", trait$PATH))
      #Read trait data and rename columns
      disease <- get_trait_data(trait, tmp_data)
      
      #Get IVs of trait data
      ivs.trait <- get_ivs(exposure=disease, exposure.name=trait$TRAIT, exposure.file=trait$FILE, outcome=t2d, out.path=trait$PATH, 
                           ref.bfile.path=ref.bfile.path, ancestry=pop, api.token=api.token, rsid.col=trait$rsid.col, pval.col=trait$pval.col)
      t2d.ivs.file <- "cluster_ivs.txt"
      
      #Check if t2d ivs need proxy
      t2d.ivs <- get_t2d_ivs(exposure=t2d, exposure.name="t2d", outcome=disease, out.path=trait$PATH, t2d.ivs.file, api.token, ancestry=pop)

      all.ivs <- unique(c(ivs.trait, t2d.ivs))
      
      #Output
      data.table::fwrite(disease[rs_id %in% all.ivs], paste0(trait$PATH,  trait$TRAIT, "_", pop, "_forMR.txt"))
      data.table::fwrite(t2d[rs_id %in% all.ivs], paste0(trait$PATH, "t2d_", pop, "_forMR.txt"))
      rm(disease)
      rm(ivs.trait)
    }
  }
}

comorbid.MR <- function(analysis, datasets.dt, path.code, ref.bfile){
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
      
      forward <- wrap.MR(exposure.file=paste0(trait$PATH, trait$TRAIT, "_", pop, "_forMR.txt"),
                         outcome.file=paste0(trait$PATH, "t2d_", pop, "_forMR.txt"),
                         exposure.name=trait$TRAIT,
                         outcome.name="t2d",
                         ivs=data.table::fread(paste0(trait$PATH, trait$TRAIT, "_", pop, "_ivs.txt"))$clumped.ivs,
                         ref.bfile.path)
      
      for(a in analysis){
        if(a=="all") t2d.ivs <- t2d.all.ivs$final.ivs
        else if (a=="one_per_locus") t2d.ivs <- t2d.all.ivs[clumped.ivs %in% one.per.locus$clumped.ivs, final.ivs]
        else if (a=="clumped") t2d.ivs <- t2d.all.ivs[clumped.ivs %in% clumped$clumped.ivs, final.ivs]
        else if (a=="cluster_specific") t2d.ivs <- t2d.all.ivs[clumped.ivs %in% clumped$clumped.ivs, final.ivs]   
        
        backward <-  wrap.MR(exposure.file=paste0(trait$PATH, "t2d_", pop, "_forMR.txt"),  
                             outcome.file=paste0(trait$PATH, trait$TRAIT, "_", pop, "_forMR.txt"),
                             exposure.name="t2d",
                             outcome.name=trait$TRAIT,
                             ivs=t2d.ivs,
                             ref.bfile.path,
                             correlated=ifelse(a=="all", TRUE, FALSE))
        
        cluster.mr <- data.table::data.table()
        for (c in unique(t2d.all.ivs$`Cluster assignment`)){
          if(c=="Liver/lipid metabolism") next  # do not run for the cluster that has only 2 SNPs

          if(a=="all") ivs.cluster <- t2d.all.ivs[`Cluster assignment`==c, final.ivs]
          else if (a=="one_per_locus") ivs.cluster <- t2d.all.ivs[`Cluster assignment`==c & clumped.ivs %in% one.per.locus$clumped.ivs, final.ivs]
          else if (a=="clumped") ivs.cluster <- t2d.all.ivs[`Cluster assignment`==c & clumped.ivs %in% clumped$clumped.ivs, final.ivs]
          else if (a=="cluster_specific") {
            if(!file.exists(paste0("t2d_ivs/", sub(" ", "_", sub("\\/", "_", c)), "_index_snps_clumped_t2d_", pop, ".txt"))) next
            cluster.ivs <- data.table::fread(paste0("t2d_ivs/", sub(" ", "_", sub("\\/", "_", c)), "_index_snps_clumped_t2d_", pop, ".txt"))

            ivs.cluster <- tryCatch(t2d.all.ivs[`Cluster assignment`==c & clumped.ivs %in% cluster.ivs$clumped.ivs, final.ivs],
                            error=function(e) e)
            if(inherits(ivs.cluster, "error")) {
              print("No clumped ivs for this cluster")
              next
            }
          }

          cluster.mr <- rbind(cluster.mr, wrap.MR(exposure.file=paste0(trait$PATH, "t2d_", pop, "_forMR.txt"),  # or "_all_forMR.txt"
                                                  outcome.file=paste0(trait$PATH, trait$TRAIT, "_", pop, "_forMR.txt"),  # or "_all_forMR.txt"
                                                  exposure.name="t2d",
                                                  outcome.name=trait$TRAIT,
                                                  ivs=ivs.cluster,
                                                  cluster.name=c,
                                                  ref.bfile.path,
                                                  correlated=ifelse(a=="all", TRUE, FALSE)), fill=TRUE)
        }
        res.dt <- rbind(forward, backward, cluster.mr, fill=TRUE)
        res.dt$ancestry <- pop
        if(nrow(res.dt)==0) return(res.dt)
        
        data.table::setnames(res.dt, old=c("b", "se", "pval", "or", "or_lci95", "or_uci95", "egger_intercept_se", "egger_intercept_pval"),
                                     new=c("beta", "standard_error", "p_value", "odds_ratio", "ci_lower", "ci_upper", "egger_intercept_standard_error", "egger_intercept_pvalue"))
        
        res.dt[is.na(cluster), cluster:="Full"]
        data.table::fwrite(res.dt[, `:=` (analysis=ifelse(is.na(a), "cluster_specific", a), ancestry=pop)], paste0("ResultsMR/new/all_mr_results_", trait$TRAIT, "_", pop, ifelse(is.na(a), "_cluster_specific", paste0("_", a)), ".txt"), na = "NA", sep = "\t", quote = F)
      }
    }
  }
}

meta.MR <- function(analysis, datasets.dt, trait){
  meta.dt <- data.table::data.table()
  for (pop in unique(datasets.dt$ANCESTRY)){
    for (a in analysis){
      res.dt <- data.table::fread(paste0("ResultsMR/new/all_mr_results_", trait, "_", pop, "_", a, ".txt"))
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
  
  data.table::fwrite(res.meta, paste0("ResultsMR/new/meta_analysis_", trait, "_results.txt"))
  return(res.meta)
}

meta.reverse.MR <- function(analysis, datasets.dt, trait){
  meta.dt <- data.table::data.table()
  for (pop in unique(datasets.dt$ANCESTRY)){
    for (a in analysis){
      res.dt <- data.table::fread(paste0("ResultsMR/new/all_mr_results_", trait, "_", pop, "_", a, ".txt"))
      meta.dt <- rbind(meta.dt, 
                       res.dt[method %in% c("Wald ratio", "Inverse variance weighted"), .(outcome, exposure, cluster, ancestry, cluster, analysis, method, nsnp, beta, standard_error, p_value)], 
                       fill=TRUE)
    }
  }
  
  reverse.res.meta <- data.table::data.table()
  for (a in analysis){
    dt <- meta.dt[exposure==trait & analysis==a]
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

    reverse.res.meta <- rbind(reverse.res.meta, dt, fill=TRUE)
  }
  data.table::fwrite(reverse.res.meta, paste0("ResultsMR/new/meta_analysis_reverse_", trait, "_results.txt"))
  return(reverse.res.meta)
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

#Number of cases and controls in each of the ancestries (numbers obtained from Kostas)
info.T2DGGI <- data.frame(ancestry = c("AFR", "EAS", "EUR", "AMR", "SAS"),
                          n.cases = c(50251, 88109, 242283, 29375, 16832),
                          n.controls = c(103909, 339395, 1569734, 59368, 33767))

################################################################################
#----------------------------------- MAIN -------------------------------------#
################################################################################
pop="EUR"
api.token <- "d1fe6afcb763"

# Preprocess data
preprocess.data(datasets.dt=info.datasets,
                project.path=pathMR,
                tmp_data=tmp_dir,
                info.t2d=info.T2DGGI,
                t2d.path=path_T2DGGI,
                path.ref=path_1kG)


# Run MR
comorbid.MR(analysis=c("clumped", "one_per_locus", "all", "cluster_specific"),
            datasets.dt=info.datasets[TRAIT=="COPD"],
            path.code=path_code,
            ref.bfile=path_1kG)

# Run meta-analyses
for(t in unique(info.datasets[, if(.N>1) .SD, by=TRAIT]$TRAIT)){
  meta.MR(analysis=c("all", "one_per_locus", "cluster_specific", "clumped"),
          datasets.dt=info.datasets[TRAIT==t],
          trait=t)
  
  meta.rerverse.MR(analysis=c("all"),
          datasets.dt=info.datasets[TRAIT==t],
          trait=t)
}