wd <- getwd()
library(Seurat)
library(dplyr)
library(ggplot2)
library(sessioninfo)
session_info(to_file = "session.log")

s.1.obj <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Blood/5619-RG-14/filtered_feature_bc_matrix") %>% CreateSeuratObject(project="PBMC_Naive_First")
s.2.obj <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Blood/6421_RG_3/filtered_feature_bc_matrix") %>% CreateSeuratObject(project="PBMC_Naive_Second")

obj.list <- list(s.1.obj, s.2.obj)

## SCTransform
sum(unlist(lapply(obj.list, function(X){ncol(X)})))
obj.list <- lapply(obj.list, function(X){
    X <- PercentageFeatureSet(X, pattern = "mt-", col.name="percent.mt")
    X <- subset(X, subset=nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)
    X <- SCTransform(X,vars.to.regress='percent.mt',verbose=T)
})
sum(unlist(lapply(obj.list, function(X){ncol(X)})))

## CCA Integration
obj.features <- SelectIntegrationFeatures(obj.list, nfeatures=3000)
obj.list <- PrepSCTIntegration(object.list=obj.list, anchor.features=obj.features)

anchors <- FindIntegrationAnchors(object.list=obj.list, normalization.method="SCT",
                                  anchor.features=obj.features)
Final_Seurat <- IntegrateData(anchorset=anchors, normalization.method="SCT")

## Stats for Table
median(Final_Seurat$nCount_RNA)
median(Final_Seurat$nFeature_RNA)

Final_Seurat <- Final_Seurat %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:25, return.model=TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:25) %>%
    FindClusters(resolution = 0.5)

load('mouse.cc.genes.R')
s.genes <- mouse.cc.genes$s.genes
g2m.genes <- mouse.cc.genes$g2m.genes
DefaultAssay(Final_Seurat) <- 'RNA'
Final_Seurat <- NormalizeData(Final_Seurat)
Final_Seurat <- CellCycleScoring(Final_Seurat, s.features = s.genes,
                                 g2m.features = g2m.genes, set.ident = TRUE)

## Cluster Numbers in order for colors to roughly match in Manuscript,
## i.e. change seurat cluster number order to Cell Type number order
c.type <- Final_Seurat$seurat_clusters %>% as.character() %>%
    recode(
        '0'='0', '1'='1', '2'='2', '3'='11', '4'='10', '5'='5', '6'='6', '7'='13', '8'='7',
        '9'='8', '10'='3', '11'='14', '12'='9', '13'='20', '14'='16', '15'='12', '16'='15',
        '17'='4', '18'='17', '19'='18', '20'='19',
        )
Final_Seurat$Clus_Num_Vis <- factor(c.type, levels=c('0','1','2','3','4','5','6','7','8','9','10',
                                                     '11','12','13','14','15','16','17','18','19','20'))

c.type <- Final_Seurat$Clus_Num_Vis %>% as.character() %>%
    recode(
        '0' = 'BC1',
        '1' = 'BC2',
        '2' = 'BC3',
        '3' = 'BC4',
        '4' = 'Platelet',
        '5' = 'TC1',
        '6' = 'TC2',
        '7' = 'TC3',
        '8' = 'TC4',
        '9' = 'TC5',
        '10' = 'NK',
        '11' = 'iGC',
        '12' = 'mGC',
        '13' = 'Mo_1',
        '14' = 'Mo_2',
        '15' = 'Mac_0',
        '16' = 'DC',
        '17' = 'PB',
        '18' = 'prol.Lk',
        '19' = 'Mast',
        '20' = 'UI')
Final_Seurat$Manuscript <- factor(c.type,
                                  levels= c('BC1', 'BC2', 'BC3', 'BC4', 'Platelet',
                                            'TC1', 'TC2', 'TC3', 'TC4', 'TC5',
                                            'NK', 'iGC', 'mGC', 'Mo_1', 'Mo_2', 'Mac_0', 'DC', 'PB',
                                            'prol.Lk', 'Mast', 'UI'))

## ## Save full object - Which will be downloadable on GEO
## setwd('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit')
## saveRDS(Final_Seurat, file='Peripheral_Blood.RDS')
## setwd(wd)


