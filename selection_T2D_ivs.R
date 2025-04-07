################################################################################
#--------------------------------- FUNCTIONS ----------------------------------#
################################################################################
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

#Number of cases and controls in each of the ancestries
info.T2DGGI <- data.frame(ancestry = c("AFR", "EAS", "EUR", "AMR", "SAS"),
                          n.cases = c(50251, 88109, 242283, 29375, 16832),
                          n.controls = c(103909, 339395, 1569734, 59368, 33767))

################################################################################
#----------------------------------- MAIN -------------------------------------#
################################################################################
index.snps <- data.table::as.data.table(readxl::read_excel("index_snps_per_region.xlsx"))
merge(index.snps[, .(`Index SNV`)], cluster.snps)

########################## One IV per locus (611 IVs) #########################
# Method:
# - Select one IV per independent locus: 611 IVs
# - For the remaining index variants (1289-611):
#   - Check LD with the 611 IVs in all 5 major ancestries 
#   - If not in strict LD (r2>0.001) in any ancestry with any of the 611 IVs --> new IV

cluster.snps <- data.table::as.data.table(readxl::read_excel("cluster_ivs.xlsx"))
one.per.locus <- data.table::as.data.table(index.snps[, .SD[which.min(MR_MEGA_P)], by=Locus]$`Index SNV`)
one.per.locus <- merge(one.per.locus, cluster.snps, all.x=TRUE, by.x="V1", by.y="Index SNV")
colnames(one.per.locus) <- c("clumped.ivs", "cluster")
data.table::fwrite(one.per.locus, "one_ivs_per_locus.txt")

# Look at 678 remaining index SNPs
remaning.snps <- setdiff(index.snps$`Index SNV`, one.per.locus$clumped.ivs)
remaning.snps <- data.table::data.table(clumped.ivs=remaning.snps[!is.na(remaning.snps)])
remaning.snps <- merge(remaning.snps, cluster.snps, all.x=TRUE, by.x="clumped.ivs", by.y="Index SNV")
colnames(remaning.snps) <- c("clumped.ivs", "cluster")

not.in.ld <- NULL

for(pop in c("EUR", "AFR", "EAS", "AMR", "SAS")){
  # Produce LD matrix between remaning.snps and one.per.locus
  ld.matrix <- ieugwasr::ld_matrix_local(index.snps$`Index SNV`, 
                                         bfile = paste0(path_1kG, pop), 
                                         plink_bin = genetics.binaRies::get_plink_binary(), 
                                         with_alleles = TRUE)
  # Subset to have remaning.snps as rows and one.per.locus as columns
  selected_rows <- rownames(ld.matrix)[grep(paste(remaning.snps$clumped.ivs, collapse = "|"), rownames(ld.matrix))]
  selected_cols <- colnames(ld.matrix)[grep(paste(one.per.locus$clumped.ivs, collapse = "|"), colnames(ld.matrix))]
  ld.matrix2 <- ld.matrix[selected_cols, selected_cols]
  
  # Check if any row where all entries are <0.001
  indep.rows <- apply(ld.matrix2, 1, function(row) {
    diag_index <- match(rownames(ld.matrix2), colnames(ld.matrix2))
    any(row[-diag_index] > 0.01)
  })
  
  # Subset matrix based on selected rows
  result_matrix <- ld.matrix2[indep.rows, ]
  nrow(result_matrix)
  if(nrow(result_matrix)>0) not.in.ld <- c(not.in.ld, rownames(result_matrix))
}

########################## Clumping of IVs #########################
# Method:
# - Clump all 1,289 T2D IVs

for(pop in c("EUR", "AFR", "EAS", "AMR", "SAS")){ 
  t2d <- get_t2d_data(info.t2d=info.T2DGGI, path_T2DGGI, ancestry_T2DGGI=pop)
  all.ld.clumped <- ieugwasr::ld_clump(dplyr::tibble(rsid=t2d$rs_id, pval=t2d$p_value),
                                                 plink_bin = genetics.binaRies::get_plink_binary(),
                                                 bfile = paste0(path_1kG, pop),
                                                 clump_p=5e-8)
  data.table::fwrite(data.table::data.table(clumped.ivs=all.ld.clumped$rsid), paste0("t2d_ivs/all_snps_clumped_t2d_", pop, ".txt"))
 
  index.ld.clumped <- ieugwasr::ld_clump(dplyr::tibble(rsid=t2d[rs_id %in% index.snps$`Index SNV`, rs_id], pval=t2d[rs_id %in% index.snps$`Index SNV`, p_value]),
                                                 plink_bin = genetics.binaRies::get_plink_binary(),
                                                 bfile = paste0(path_1kG, pop),
                                                 clump_p=5e-8)
  clumped.t2d.ivs <- data.table::data.table(clumped.ivs=index.ld.clumped$rsid)
  clumped.t2d.ivs <- merge(clumped.t2d.ivs, cluster.snps, by.x="clumped.ivs", by.y="Index SNV", all.x=TRUE)
  data.table::fwrite(clumped.t2d.ivs, paste0("t2d_ivs/index_snps_clumped_t2d_", pop, ".txt"))
}


dt <- data.table::fread("data/oa/t2d_EUR_clumped_ivs.txt")
dt <- merge(dt, cluster.snps, by.x="clumped.ivs", by.y="Index SNV", all.x=TRUE)
data.table::fwrite(dt, "data/oa/t2d_EUR_clumped_ivs.txt")

########################## Cluster-stratified clumping of IVs #########################
# Method:
# - Group IVs per cluster
# - Clump IVs per cluster

clusters <- data.table::as.data.table(readxl::read_excel("cluster_ivs.xlsx"))

for(pop in c("EUR", "AFR", "EAS", "AMR", "SAS")){
  t2d <- get_t2d_data(info.t2d=info.T2DGGI, path_T2DGGI, ancestry_T2DGGI=pop)
  for(c in unique(clusters$`Cluster assignment`)){
    cluster.snps <- clusters[`Cluster assignment`==c, `Index SNV`]

    index.ld.clumped <- tryCatch(ieugwasr::ld_clump(dplyr::tibble(rsid=t2d[rs_id %in% cluster.snps, rs_id], pval=t2d[rs_id %in% cluster.snps, p_value]),
                                         plink_bin = genetics.binaRies::get_plink_binary(),
                                         bfile = paste0(path_1kG, pop),
                                         clump_p=5e-8), error=function(e) e)
    if(inherits(index.ld.clumped, "error")){
      print("Clumping removes all SNPs")
      next
    }
    clumped.t2d.ivs <- data.table::data.table(clumped.ivs=index.ld.clumped$rsid)
    clumped.t2d.ivs <- merge(clumped.t2d.ivs, clusters[`Cluster assignment`==c], by.x="clumped.ivs", by.y="Index SNV", all.x=TRUE)
    data.table::fwrite(clumped.t2d.ivs, paste0("t2d_ivs/", sub(" ", "_", sub("\\/", "_", c)), "_index_snps_clumped_t2d_", pop, ".txt"))
  }
}