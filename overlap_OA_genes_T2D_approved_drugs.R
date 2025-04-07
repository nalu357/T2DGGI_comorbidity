library(data.table)

indications <- as.data.table(readRDS(file="OT_move_to_lustre/data_parsed/OT_indications.rds", refhook=NULL))
indications[grepl("type 2 diabetes mellitus", approved_indications)]
t2d.approved.drugs <- unique(unlist(indications[grepl("type 2 diabetes mellitus", approved_indications)]$ID))

OT <- as.data.table(readRDS(file="OT_move_to_lustre/analysis/OT_by_gene.rds", refhook=NULL))
OT[ID %in% t2d.approved.drugs, .N]
t2d.drugs.genes <- unique(OT[ID %in% t2d.approved.drugs, Gene])  # 69 unique genes targeted by T2D approved drugs
t2d.drugs.ensg <- unique(OT[ID %in% t2d.approved.drugs, ENSG_id])  # 69 unique genes targeted by T2D approved drugs
unique(OT[ID %in% t2d.approved.drugs, drug_name])  # 45 unique approved drugs for T2D

setwd("C:/Users/AnaLuizaArruda/OneDrive - Helmholtz Zentrum München/Projects/T2DGGI/T2DGGI_MR/project_2/OpenTargets/")
result.dt <- data.table()
disease.dt <- as.data.table(readxl::read_excel("prio_genes_table.xlsx"))

for(i in 1:nrow(disease.dt)){
  d <- disease.dt[i]
  if(d$Disease %in% c("Depression EUR", "Osteoporosis", "Polycystic ovary syndrome", "Rheumatoid arthritis", "Back pain", "Carpal tunnel syndrome", "Obsessive compulsory disorder")) next
  cur.disease <- as.data.table(readxl::read_excel(paste0(d$File, ".xlsx")))
  if(!is.na(d$gene.id.col)) cur.disease.gene.id <- cur.disease[get(d$gene.id.col) %in% t2d.drugs.ensg, get(d$gene.id.col)]
  else cur.disease.gene.id <- NULL
  
  if(!is.na(d$gene.name.col)) cur.disease.gene.name <- cur.disease[get(d$gene.name.col) %in% t2d.drugs.genes, get(d$gene.name.col)]
  else cur.disease.gene.name <- NULL
  
  result.dt <- rbind(result.dt, data.table(disease=d$Disease, gene.id=cur.disease.gene.id, gene.name=cur.disease.gene.name), fill=TRUE)
}

result.dt[, drugs:=OT[ID %in% t2d.approved.drugs & Gene %in% gene.name, drug_name], by=gene.name]
result.dt <- merge(result.dt, OT[ID %in% t2d.approved.drugs & Gene %in% unique(result.dt$gene.name), .(Gene, ENSG_id, drug_name)], by.x="gene.name", by.y="Gene", all.x=TRUE)

OT[drug_name %in% unique(result.dt$drug_name)]