wd <- getwd()
library(Seurat)
library(dplyr)
library(ggplot2)

RG_3088_01 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_1/3088-RG-1/filtered_feature_bc_matrix")
RG_3088_02 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_1/3088-RG-2/filtered_feature_bc_matrix")
RG_3088_03 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_1/3088-RG-3/filtered_feature_bc_matrix")
RG_3088_04 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_1/3088-RG-4/filtered_feature_bc_matrix")
RG_4884_05 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_1/4884-RG-5/filtered_feature_bc_matrix")
s_3088_01_obj <- CreateSeuratObject(counts = RG_3088_01, project = "Giger_UCLA-3088_01")
s_3088_02_obj <- CreateSeuratObject(counts = RG_3088_02, project = "Giger_UCLA-3088_02")
s_3088_03_obj <- CreateSeuratObject(counts = RG_3088_03, project = "Giger_UCLA-3088_03")
s_3088_04_obj <- CreateSeuratObject(counts = RG_3088_04, project = "Giger_UCLA-3088_04")
s_4884_05_obj <- CreateSeuratObject(counts = RG_4884_05, project = "Giger_UCLA-4884_05")
obj.list <- lapply(ls(pattern="_obj$"), function(X){get(X)})

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

rm(list=ls(pattern="^a|^o|^s|^R"))

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
        '2'='GC',
        '3'='Mo',
        '4'='Mo',
        '5'='MoDC', ## Need pDC
        '6'='UI',
        '7'='Mac',
        '8'='Mac',
        '9'='Mac',
        '10'='dMES',
        '11'='Mac',
        '12'='Mo',
        '13'='prol.MES',
        '14'='SC',
        '15'='pMES',
        '16'='TC',
        '17'='EC',
        '18'='GC',
        '19'='DCx',
        '20'='vSMC_PC',
        '21'='cDC',
        '22'='Mast')
Final_Seurat$Cell.Type <- c.type

## TC and NK are clustering together
## split out NK as those that express either Ncr1 or Klra7
temp <- subset(Final_Seurat, subset = Cell.Type == "TC" & (Ncr1 > 0 | Klra7 > 0))
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'NK'

## Clear Fibroblast cells in the dMES
## Cutoff using the UMAP, couldn't determine a better way
temp <- subset(Final_Seurat, subset = Cell.Type == 'dMES' & UMAP_2 > 6)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'Fb'

## Small population of pDC cells clustering off of the MoDC cluster
temp <- subset(Final_Seurat, subset = Cell.Type == "MoDC" & UMAP_2 < -1.5)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'pDC'

Final_Seurat$Cell.Type <-factor(c.type,
                                levels= c('GC','Mo','Mac','MoDC','cDC','pDC','DCx','Mast','NK','TC',
                                          'eMES','pMES','dMES','prol.MES','Fb',
                                          'SC','EC','vSMC_PC'))
p1 <- DimPlot(Final_Seurat, reduction = "umap", label = TRUE, pt.size=.25, group.by='Cell.Type') +
    theme(legend.position = "none")

## Cluster Numbers in order for colors to roughly match in Manuscript,
## i.e. change seurat cluster number order to Cell Type number order
### Make cluster numbers for the dMES and pDC as above
c.type <- as.character(Final_Seurat$seurat_clusters)

## split out pDC
temp <- subset(Final_Seurat, subset = seurat_clusters == "5" & UMAP_2 < -1.5)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '50'

## split out Fb
temp <- subset(Final_Seurat, subset = seurat_clusters == "10" & UMAP_2 > 6)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '100'

c.type <- c.type  %>%
    recode('0'='1', '1'='18', '2'='11', '3'='0', '4'='0', '5'='6', '6'='14', '7'='5', '8'='2',
           '9'='3', '10'='16', '11'='4', '12'='0', '13'='17', '14'='22', '15'='19', '16'='13',
           '17'='20', '18'='12', '19'='8', '20'='21', '21'='7', '22'='10',
           '50'='9', '100'='15')

Final_Seurat$Clus_Num_Vis <-factor(c.type,
                                   levels=c('0','1','2','3','4','5','6','7','8','9','10',
                                            '11','12','13','14','15','16','17','18',
                                            '19','20','21','22'))

### Unique ID for each cluster
c.type <- Final_Seurat$Clus_Num_Vis %>% as.character() %>%
    recode(
        '0'='Mo',
        '1'='Mac-I',
        '2'='Mac-II',
        '3'='Mac-III',
        '4'='Mac-IV',
        '5'='Mac-V',
        '6'='MoDC',
        '7'='cDC',
        '8'='DCx',
        '9'='pDC',
        '10'='MAST',
        '11'='GC',
        '12'='GC-L',
        '13'='T_NK',
        '14'='UI-M1',
        '15'='Fb',
        '16'='dMES',
        '17'='prol.MES',
        '18'='eMES',
        '19'='pMES',
        '20'='EC',
        '21'='vSMC_PC',
        '22'='SC')
Final_Seurat$Manuscript <- c.type
Final_Seurat$Manuscript <-factor(c.type,
                                 levels= c('Mo','Mac-I','Mac-II','Mac-III','Mac-IV','Mac-V',
                                           'MoDC','DCx','cDC','pDC','MAST','GC','GC-L',
                                           'T_NK','UI-M1',
                                           'eMES','dMES','prol.MES','pMES','Fb',
                                           'EC','vSMC_PC','SC'))

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
## saveRDS(Final_Seurat, file='Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/post_injury_day_1.RDS')
