wd <- getwd()
library(Seurat)
library(dplyr)
library(ggplot2)

WT_1 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_3/Sample_WT_sciatic_3_days_1/filtered_feature_bc_matrix")
WT_2 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_3/Sample_WT_sciatic_3_days_2/filtered_feature_bc_matrix")
WT_3 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_3/Sample_WT2/filtered_feature_bc_matrix")
s1_obj <- CreateSeuratObject(counts = WT_1, project = "Giger_sciatic_Batch1")
s2_obj <- CreateSeuratObject(counts = WT_2, project = "Giger_sciatic_Batch1.2")
s3_obj <- CreateSeuratObject(counts = WT_3, project = "Giger_sciatic_Batch2")
## names and treatment to match  Manuscript 1
WT_B1 <- merge(s1_obj, y=s2_obj, add.cell.ids=c("batch1_WT_1","batch1_WT_2"))

S1 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_3/Sample_1081-RG-5/filtered_feature_bc_matrix")
s4_obj <- CreateSeuratObject(counts = S1, project = "MPZ_myelin_3d")


S1 <-Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_3/Sample_1181-RG-1/filtered_feature_bc_matrix")
S2 <-Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_3/Sample_1181-RG-2/filtered_feature_bc_matrix")
s6_obj <- CreateSeuratObject(counts = S1, project = "WT_NoDrug_1")
s7_obj <- CreateSeuratObject(counts = S2, project = "WT_NoDrug_2")

obj.list <- list(WT_B1,s3_obj,s4_obj,s6_obj,s7_obj)

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

rm(list=ls(pattern="^a|^o|^s|^S|^W"))

Final_Seurat <- Final_Seurat %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, return.model=TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 0.5)

load('mouse.cc.genes.R')
s.genes <- mouse.cc.genes$s.genes
g2m.genes <- mouse.cc.genes$g2m.genes
DefaultAssay(Final_Seurat) <- 'RNA'
Final_Seurat <- NormalizeData(Final_Seurat)
Final_Seurat <- CellCycleScoring(Final_Seurat, s.features = s.genes,
                                 g2m.features = g2m.genes, set.ident = TRUE)

## Save the UMAP coordinates to aid in splitting some clusters
temp <- as.data.frame(Embeddings(Final_Seurat[["umap"]]))
Final_Seurat[["UMAP_1"]] <- temp$UMAP_1
Final_Seurat[["UMAP_2"]] <- temp$UMAP_2

### Macro cluster IDs - that match across all time points
c.type <- Final_Seurat$seurat_clusters %>% as.character() %>%
    recode(
        '0'='Mac',
        '1'='eMES',
        '2'='Mac',
        '3'='SC',
        '4'='Mac',
        '5'='prol.MES',
        '6'='MoDC',
        '7'='EC',
        '8'='SC',
        '9'='dMES',
        '10'='Mac',
        '11'='Fb',
        '12'='vSMC_PC',
        '13'='pMES',
        '14'='Mo',
        '15'='NK', ## split out TC
        '16'='Mac',
        '17'='DCx', ## split out pDC and Mast
        '18'='cDC',
        '19'='GC',
        '20'='vSMC_PC',
        '21'='EC')
Final_Seurat$Cell.Type <- c.type

## Bottom pennisula of DCx is pDC clusters
## split out as those that express Siglech or CCr9
#temp <- subset(Final_Seurat, subset = Cell.Type == "DCx" & (Siglech > 0 | Ccr9 > 0))
temp <- subset(Final_Seurat, subset = Cell.Type == "DCx" & UMAP_1 > 3 & UMAP_2 < -9)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'pDC'

## Some cells from DCx cluster map with T_NK cluster in UMAP
## These cells look like Mast cells
temp <- subset(Final_Seurat, subset = Cell.Type == "DCx" & UMAP_1 < 3)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'Mast'

## TC and NK are clustering together (more NK this time)
## split out TC as those that express either Cd3g or Lat
temp <- subset(Final_Seurat, subset = Cell.Type == "NK" & (Cd3g > 0 | Lat > 0))
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'TC'

Final_Seurat$Cell.Type <-factor(c.type,
                                levels= c('GC','Mo','Mac','MoDC','cDC','pDC','DCx','Mast','NK','TC',
                                          'eMES','pMES','dMES','prol.MES','Fb','SC','EC','vSMC_PC'))

p1 <- DimPlot(Final_Seurat, reduction = "umap", label = TRUE, pt.size=.25, group.by='Cell.Type') +
    theme(legend.position = "none")

## Cluster Numbers in order for colors to roughly match in Manuscript,
## i.e. change seurat cluster number order to Cell Type number order
### Make cluster numbers for pDC, MAST, and prol.SC
c.type <- as.character(Final_Seurat$seurat_clusters)

## split out pDC
temp <- subset(Final_Seurat, subset = seurat_clusters == "17" & UMAP_1 > 3 & UMAP_2 < -9)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '90'

## split out Mast
temp <- subset(Final_Seurat, subset = seurat_clusters == "17" & UMAP_1 < 3)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '100'

## split out prol.SC
temp <- subset(Final_Seurat, subset = seurat_clusters == "3" & Final_Seurat$Phase %in% c('G2M','S'))
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '230'

c.type <- c.type  %>%
    recode('0'='1', '1'='16', '2'='2', '3'='22', '4'='3', '5'='15', '6'='6', '7'='18', '8'='24',
           '9'='14', '10'='4', '11'='13', '12'='20', '13'='17', '14'='0', '15'='12', '16'='5',
           '17'='8', '18'='7', '19'='11', '20'='21', '21'='19',
           '90'='9', '100'='10', '230'='23')

Final_Seurat$Clus_Num_Vis <-factor(c.type, levels=c('0','1','2','3','4','5','6','7','8','9',
                                                    '10','11','12','13','14','15','16','17',
                                                    '18','19','20','21','22','23','24'))

### Unique ID for each cluster
c.type <- Final_Seurat$Clus_Num_Vis %>% as.character() %>%
    recode(
        '0'='Mo',
        '1'='Mac1',
        '2'='Mac2',
        '3'='Mac3',
        '4'='Mac4',
        '5'='prol.Mac',
        '6'='MoDC',
        '7'='cDC',
        '8'='DCx',
        '9'='pDC',
        '10'='MAST',
        '11'='GC',
        '12'='T_NK',
        '13'='Fb',
        '14'='dMES',
        '15'='prol.MES',
        '16'='eMES',
        '17'='pMES',
        '18'='EC1',
        '19'='EC2',
        '20'='vSMC_PC1',
        '21'='vSMC_PC2',
        '22'='nmSC',
        '23'='prol.SC',
        '24'='rSC')
Final_Seurat$Manuscript <-factor(c.type,
                                 levels= c('Mo','Mac1','Mac2','Mac3','Mac4','prol.Mac',
                                           'MoDC','cDC','DCx','pDC','MAST','GC','T_NK',
                                           'Fb','dMES','prol.MES','eMES','pMES',
                                           'EC1','EC2','vSMC_PC1','vSMC_PC2','nmSC','prol.SC','rSC'))

Idents(Final_Seurat) <- 'Manuscript'
p2 <- DimPlot(Final_Seurat, reduction = "umap", label = TRUE, pt.size=.25,
              group.by='seurat_clusters') +
    theme(legend.position = "none")
p3 <- DimPlot(Final_Seurat, reduction = "umap", label = TRUE, pt.size=.25,
              group.by='Manuscript') +
    theme(legend.position = "none")
p4 <- FeaturePlot(Final_Seurat, features='Mki67')
## pdf('UMAPs.pdf', width=18, height=6)
print((p1 | p2) / (p3 | p4))
## dev.off()

## ## Save full object - Which will be downloadable on GEO
## saveRDS(Final_Seurat, file='Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/post_injury_day_3_test.RDS')
