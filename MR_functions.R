library(dplyr)

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

my.harmonise_ld_dat <- function(x, ld){
  snpnames <- do.call(rbind, strsplit(rownames(ld), split="_"))
  i1 <- snpnames[,1] %in% x$SNP
  ld <- ld[i1,i1]
  snpnames <- snpnames[i1,]
  i2 <- x$SNP %in% snpnames[,1]
  x <- x[i2,]
  snpnames <- snpnames[order(match(snpnames[,1],x$SNP)),]
  stopifnot(all(snpnames[,1]== x$SNP))
  x$effect_allele.exposure <- as.character(x$effect_allele.exposure)
  x$other_allele.exposure <- as.character(x$other_allele.exposure)
  # Set1 x and ld alleles match
  snpnames <- data.frame(snpnames, stringsAsFactors=FALSE)
  snpnames <- merge(subset(x, select=c(SNP, effect_allele.exposure, other_allele.exposure)), snpnames, by.x="SNP", by.y="X1")
  snpnames <- snpnames[match(x$SNP, snpnames$SNP),]
  snpnames$keep <- (snpnames$X2 == snpnames$effect_allele.exposure & snpnames$X3 == snpnames$other_allele.exposure) |
    (snpnames$X3 == snpnames$effect_allele.exposure & snpnames$X2 == snpnames$other_allele.exposure)
  
  # What happens if everything is gone?
  if(nrow(x) == 0)
  {
    message(" - none of the SNPs could be aligned to the LD reference panel")
    return(NULL)
  }
  
  if(any(!snpnames$keep))
  {
    message(" - the following SNPs could not be aligned to the LD reference panel: \n- ", paste(subset(snpnames, !keep)$SNP, collapse="\n - "))
  }
  
  snpnames$flip1 <- snpnames$X2 != snpnames$effect_allele.exposure
  temp1 <- x$effect_allele.exposure[snpnames$flip1]
  temp2 <- x$other_allele.exposure[snpnames$flip1]
  x$beta.exposure[snpnames$flip1] <- x$beta.exposure[snpnames$flip1] * -1
  x$beta.outcome[snpnames$flip1] <- x$beta.outcome[snpnames$flip1] * -1
  x$effect_allele.exposure[snpnames$flip1] <- temp2
  x$other_allele.exposure[snpnames$flip1] <- temp1
  
  rownames(ld) <- snpnames$SNP
  colnames(ld) <- snpnames$SNP
  
  if(any(!snpnames$keep))
  {
    message("Removing ", sum(!snpnames$keep), " variants due to harmonisation issues")
    ld <- ld[snpnames$keep, snpnames$keep]
    x <- x[snpnames$keep, ]
    
  }
  return(list(x=x, ld=ld))
}

local.data.to.MRInput <- function(data, ref.bfile.path=NA, get_correlations=FALSE){
  res <- plyr::dlply(data, c("exposure", "outcome"), function(x){
    x <- plyr::mutate(x)
    message("Converting:")
    message(" - exposure: ", x$exposure[1])
    message(" - outcome: ", x$outcome[1])
    if(get_correlations){
      message(" - obtaining LD matrix")
      ld <- ieugwasr::ld_matrix(x$SNP,
                                plink_bin = genetics.binaRies::get_plink_binary(),
                                bfile = ref.bfile.path)
      out <- my.harmonise_ld_dat(x, ld)
      if(is.null(out)) return(NULL)
      x <- out$x
      ld <- out$ld
      
      MendelianRandomization::mr_input(
        bx = x$beta.exposure,
        bxse = x$se.exposure,
        by = x$beta.outcome,
        byse = x$se.outcome,
        exposure = x$exposure[1],
        outcome = x$outcome[1],
        snps = x$SNP,
        effect_allele=x$effect_allele.exposure,
        other_allele=x$other_allele.exposure,
        eaf = x$eaf.exposure,
        correlation = ld)
    }
    else {
      MendelianRandomization::mr_input(
        bx = x$beta.exposure,
        bxse = x$se.exposure,
        by = x$beta.outcome,
        byse = x$se.outcome,
        exposure = x$exposure[1],
        outcome = x$outcome[1],
        snps = x$SNP,
        effect_allele=x$effect_allele.exposure,
        other_allele=x$other_allele.exposure,
        eaf = x$eaf.exposure)
    }
  })
  return(res)
}

add_odds_ratios <- function(dt, beta="beta", se="standard_error"){
  dt$or <- exp(dt[, get(beta)])
  dt$or_lci95 <- exp(dt[, get(beta)] - 1.96 * dt[, get(se)])
  dt$or_uci95 <- exp(dt[, get(beta)] + 1.96 * dt[, get(se)])
  return(dt)
}

run.MRPresso <- function(my.data, nbdist=1000, beta.out="beta.outcome", beta.exp="beta.exposure", se.out="se.outcome", se.exp="se.exposure"){
  dt <- MRPRESSO::mr_presso(BetaOutcome=beta.out, 
                            BetaExposure=beta.exp, 
                            SdOutcome=se.out, 
                            SdExposure=se.exp, 
                            OUTLIERtest=TRUE, 
                            DISTORTIONtest=TRUE, 
                            data=my.data, 
                            NbDistribution=nbdist,  
                            SignifThreshold=0.05)
  
  #Results of the outlier tests
  #Also keep the disortion p-value and the number of snps
  res.presso <- data.table::as.data.table(dt$`Main MR results`)[`MR Analysis`=="Outlier-corrected"]
  res.presso <- add_odds_ratios(res.presso, "Causal Estimate", "Sd")
  
  res.presso <- res.presso[, .(b=`Causal Estimate`, se=`Sd`, or, or_lci95, or_uci95, pval=`P-value`,
                               method="MR_PRESSO_IVW_Outliers", exposure=my.data$exposure[1], outcome=my.data$outcome[1],
                               Distortion_pval = dt$`MR-PRESSO results`$`Distortion Test`$Pvalue,
                               nsnp = nrow(my.data)-length(dt$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`))]
  
  #Add the raw test
  res.presso.raw <- data.table::as.data.table(dt$`Main MR results`)[`MR Analysis`=="Raw"]
  res.presso.raw <- add_odds_ratios(res.presso.raw, "Causal Estimate", "Sd")
  
  res.presso.raw <- res.presso.raw[, .(b=`Causal Estimate`, se=`Sd`, or, or_lci95, or_uci95, pval=`P-value`,
                                       method="MR_PRESSO_IVW_Raw", exposure=my.data$exposure[1], outcome=my.data$outcome[1],
                                       nsnp = nrow(my.data))]
  
  # res.presso<- run.SensitivityAnalysis(my.data, res.presso)
  return(rbind(res.presso, res.presso.raw, fill = T))
}

run.SensitivityAnalysis <- function(data, res, Radial=FALSE, plt=TRUE, het=TRUE, presso=TRUE, steiger = TRUE){
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
    # #Run Steiger directionality test
    # steiger.test <- tryCatch(
    #   TwoSampleMR::directionality_test(data),
    #   error=function(e) e
    # )
    # 
    # if(inherits(steiger.test, "error")) print("Steiger directionality test yielded error")
    # else{
    #   res$directionality.test <- steiger.test
    # }
    
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
  if(presso){
    res.presso <- tryCatch(
      run.MRPresso(data), 
      error=function(e) e
    )
    
    if(inherits(res.presso, "error")) {
      print("MR-PRESSO failed first try. Probably not enough SNPs.")
      res.presso <- tryCatch(
        run.MRPresso(data, nbdist=1500), 
        error=function(e) e
      )
      if(inherits(res.presso, "error")) print("MR-PRESSO failed again.")
      else res <- rbind(res, res.presso, fill=TRUE)
    }else{
      res <- rbind(res, res.presso, fill=TRUE)
    }
  }
  if(Radial){
    #Run radial filtering, esp. for metabolites
    #Can have warning or error if model does not converge
    radial.outliers <- tryCatch(ivw_radial(data, alpha = 0.05, weights = 1), error = identity, warning = identity)
    #If there are SNPs left after filtering, proceed with the test
    if(is(radial.outliers, "error") | is(radial.outliers, "warning")){
      res.radial <- data.table(method = "Radial Inverse variance weighted", nsnp = NA)
    }else{
      if(nrow(data) == length(radial.outliers$outliers$SNP)){
        res.radial <- data.table(method = "Radial Inverse variance weighted", nsnp = 0)
      }else{
        if(is.null(ncol(radial.outliers$outliers))){
          #No outlier detected
          res.radial <- data.table(method = "Radial Inverse variance weighted", nsnp = nrow(data))
        }else{
          data.radial <- subset(data, !(SNP %in% radial.outliers$outliers$SNP))
          res.radial <- TwoSampleMR::mr(data.radial, method = "mr_ivw")
          res.radial <- data.table::as.data.table(TwoSampleMR::generate_odds_ratios(res.radial))
          res.radial[, `:=`(lo_ci=NULL, up_ci=NULL, id.exposure=NULL, id.outcome=NULL)]
          #Compute heterogeneity
          het <- data.table::as.data.table(TwoSampleMR::mr_heterogeneity(data.radial, method = "mr_ivw"))
          res.radial <- merge(res.radial, het[, .(exposure, outcome, method, Q, Q_df, Q_pval)],
                              by=c("exposure", "outcome", "method"), all.x=TRUE)
          res.radial$method <- "Radial Inverse variance weighted"
        }
      }
    }
    res <- rbind(res, res.radial, fill = TRUE)
  }
  
  return(res)
}

run.TwoSampleMR <- function(data, ref.bfile.path, Radial, correlated, cluster.name, presso){
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
  
    if(correlated & !is.na(cluster.name)){
      data2 <- local.data.to.MRInput(data, ref.bfile.path, get_correlations=TRUE)
      res.cor <- tryCatch(MendelianRandomization::mr_ivw(data2[[1]], correl = TRUE),
                          error=function(e) e)
      
      if(inherits(res.cor, "error")) return(res)
      
      res.cor <- data.table::data.table(outcome=res.cor$Outcome,
                                        exposure=res.cor$Exposure,
                                        method="Correlated IVW",
                                        nsnp=res.cor$SNPs,
                                        b=res.cor$Estimate,
                                        se=res.cor$StdError,
                                        pval=res.cor$Pvalue,
                                        or=exp(res.cor$Estimate),
                                        or_lci95=exp(res.cor$CILower),
                                        or_uci95=exp(res.cor$CIUpper))
      res <- rbind(res, res.cor)
    }
    
    # Run sensitivity analysis
    res <- run.SensitivityAnalysis(data, res, Radial, presso)
    
    return(res)
  }
}

wrap.MR <- function(exposure.file, outcome.file, exposure.name, outcome.name, ivs, ref.bfile.path, correlated=F, Radial=F, cluster.name=NA, presso=T){
  res.dt <- data.table::data.table()
  if(length(ivs)==0) return(res.dt)
  
  #----------- Get exposure -----------
  exposure <- get.data(data.table::fread(exposure.file), data.dir="exposure")
  exposure <- exposure[exposure$SNP %in% ivs,]
  if(nrow(exposure)==0) return(res.dt)
  exposure$exposure <- exposure.name
  
  #----------- Get outcome -----------
  outcome <- get.data(data.table::fread(outcome.file), data.dir="outcome")
  outcome$outcome <- outcome.name

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
    print("No SNP left after data harmonization, not able to run MR :(")
    return(res.dt)
  }
  data$id.exposure <- data$exposure
  data$id.outcome <- data$outcome      
  if(!is.na(cluster.name)) data$cluster <- cluster.name 
  else data$cluster <- "Full"
  data <- subset(data, mr_keep==TRUE)
  
  #----------- Run TwoSampleMR with sensitivity analyses -----------
  res.mr <- run.TwoSampleMR(data, ref.bfile.path, Radial, correlated, cluster.name, presso=presso)
  res.dt <- rbind(res.dt, res.mr, fill=TRUE)
  
  #Compute the I2 statistics
  res.dt$I2 <- (res.dt$Q-res.dt$Q_df)/res.dt$Q
  
  if(!is.na(cluster.name)) res.dt$cluster=cluster.name
  else res.dt$cluster="Full"
  
  return(res.dt)
}

