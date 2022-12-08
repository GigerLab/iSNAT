wd <- getwd()
library(Seurat)
library(dplyr)
library(ggplot2)

RG_5 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_7/2405-RG-5/filtered_feature_bc_matrix")
RG_6 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_7/2405-RG-6/filtered_feature_bc_matrix")
RG_7 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_7/2405-RG-7/filtered_feature_bc_matrix")
RG_8 <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_7/2405-RG-8/filtered_feature_bc_matrix")
ss5_obj <- CreateSeuratObject(counts = RG_5, project = "Giger_Sciatic_1")
ss6_obj <- CreateSeuratObject(counts = RG_6, project = "Giger_Sciatic_2")
ss7_obj <- CreateSeuratObject(counts = RG_7, project = "Giger_Sciatic_3")
ss8_obj <- CreateSeuratObject(counts = RG_8, project = "Giger_Sciatic_4")

S2 <-Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_7/Sample_1081-RG-7/filtered_feature_bc_matrix")
S3 <-Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_7/Sample_1081-RG-8/filtered_feature_bc_matrix")
s253_2_obj <- CreateSeuratObject(counts = S2, project = "MPZ_myelin_7d")
s253_3_obj <- CreateSeuratObject(counts = S3, project = "MPZ_myelin_depleted_7d")

RG_13 <-Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_7/2856-RG-13/filtered_feature_bc_matrix")
RG_14 <-Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_7/2856-RG-14/filtered_feature_bc_matrix")
s13_obj <- CreateSeuratObject(counts = RG_13, project = "Giger_UCLA-24_13")
s14_obj <- CreateSeuratObject(counts = RG_14, project = "Giger_UCLA-24_14")

obj.list <- lapply(ls(pattern="_obj$"), function(X){get(X)})

## SCTransform
sum(unlist(lapply(obj.list,function(X){length(X$orig.ident)})))
obj.list <- lapply(obj.list, function(X){
    X <- PercentageFeatureSet(X, pattern = "mt-", col.name="percent.mt")
    X <- subset(X, subset=nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)
    X <- SCTransform(X,vars.to.regress='percent.mt',verbose=T)
})
sum(unlist(lapply(obj.list,function(X){length(X$orig.ident)})))

## CCA Integration
obj.features <- SelectIntegrationFeatures(obj.list, nfeatures=3000)
obj.list <- PrepSCTIntegration(object.list=obj.list, anchor.features=obj.features)

anchors <- FindIntegrationAnchors(object.list=obj.list, normalization.method="SCT",
                                  anchor.features=obj.features)
Final_Seurat <- IntegrateData(anchorset=anchors, normalization.method="SCT")

## Stats for Table
median(Final_Seurat$nCount_RNA)
median(Final_Seurat$nFeature_RNA)

rm(list=ls(pattern="^a|^o|^s|^S|^R"))

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
        '0'='SC',
        '1'='dMES',
        '2'='dMES',
        '3'='MoDC',
        '4'='Mac',
        '5'='Mac',
        '6'='TC',
        '7'='Fb',
        '8'='Mac',
        '9'='Mac',
        '10'='SC',
        '11'='SC',
        '12'='NK',
        '13'='EC',
        '14'='vSMC_PC',
        '15'='SC',
        '16'='pMES',
        '17'='vSMC_PC',
        '18'='eMES',
        '19'='MoDC',
        '20'='Mo', ## GC needs labeled
        '21'='vSMC_PC',
        '22'='cDC',
        '23'='Mac',
        '24'='dMES',
        '25'='NoLabel',
        '26'='pDC',
        '27'='NK',
        '28'='MoDC',
        '29'='Mast',
        '30'='vSMC_PC',
        '31'='cDC')
Final_Seurat$Cell.Type <- c.type

## Some GC cells in the Mo cluster
## Section of Mo cluster expresses Cxcr2, S100a9, or Lcn2. Label as GC
temp <- subset(Final_Seurat, subset = Cell.Type == "Mo" & UMAP_1 > -6.1 & UMAP_2 < 1)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'GC'

## SC cluster region with expression of Mki67 will be labeled as proliferating
temp <- subset(Final_Seurat, subset = Cell.Type == "SC" & UMAP_2 < -11)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'prol.MES'

## A section of the unlabeled area appears to be DCx cells
temp <- subset(Final_Seurat, subset = Cell.Type == "NoLabel" & UMAP_2 > -3)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'DCx'

## NoLabel will be grey
Final_Seurat$Cell.Type <-factor(c.type,
                                levels= c('GC','Mo','Mac','MoDC','cDC','pDC','DCx','Mast','NK','TC',
                                          'eMES','pMES','dMES','prol.MES','Fb',
                                          'SC','EC','vSMC_PC'))
p1 <- DimPlot(Final_Seurat, reduction = "umap", label = TRUE, pt.size=.25, group.by='Cell.Type') +
    theme(legend.position = "none")

### Make cluster labels in the order of the manuscript tables
### Same order between time points but not same cluster numbers because clusters come and go
### split DCx, GC, prol.MES as above
c.type <- as.character(Final_Seurat$seurat_clusters)

## A section of cluster 25 (NoLabel) looks like DCx cells
temp <- subset(Final_Seurat, subset = seurat_clusters == "25" & UMAP_2 > -3)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '100'

## Section of cluster 20 (Monocytes) expresses Cxcr2, S100a9, or Lcn2. Label as GC
temp <- subset(Final_Seurat, subset = seurat_clusters == "20" & UMAP_1 > -6.1 & UMAP_2 < 1)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '140'

## prol.MES
temp <- subset(Final_Seurat, subset = seurat_clusters == "11" & UMAP_2 < -11)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '230'

c.type <- c.type  %>%
    recode('0'='30', '1'='19', '2'='20', '3'='7', '4'='1', '5'='3', '6'='16', '7'='18', '8'='4',
           '9'='2', '10'='33', '11'='31', '12'='15', '13'='25', '14'='26', '15'='32', '16'='24',
           '17'='27', '18'='23', '19'='9', '20'='0', '21'='28', '22'='8', '23'='6', '24'='21',
           '25'='34', '26'='11', '27'='17', '28'='5', '29'='13', '30'='29', '31'='12',
           '100'='10', '140'='14', '230'='22')


Final_Seurat$Clus_Num_Vis <-factor(c.type, levels=c('0','1','2','3','4','5','6','7','8','9',
                                                    '10','11','12','13','14','15','16','17',
                                                    '18','19','20','21','22','23','24','25',
                                                    '26','27','28','29','30','31','32','33',
                                                    '34'))

### Unique ID for each cluster
c.type <- Final_Seurat$Clus_Num_Vis %>% as.character() %>%
    recode(
        '0'='Mo',
        '1'='Mac-a',
        '2'='Mac-b',
        '3'='Mac-c',
        '4'='Mac-d',
        '5'='Mac-e',
        '6'='prol.Mac',
        '7'='MoDC',
        '8'='cDC1',
        '9'='cDC2',
        '10'='DCx',
        '11'='pDC',
        '12'='prol.DC',
        '13'='MAST',
        '14'='GC',
        '15'='NK',
        '16'='TC',
        '17'='UI-M',
        '18'='Fb',
        '19'='dMES1',
        '20'='dMES2',
        '21'='dMES3',
        '22'='prol.MES',
        '23'='eMES',
        '24'='pMES',
        '25'='EC',
        '26'='vPC1',
        '27'='vPC2',
        '28'='aPC',
        '29'='prol.PC',
        '30'='nmSC',
        '31'='prol.SC',
        '32'='rSC',
        '33'='mSC',
        '34'='UI-SC')
Final_Seurat$Manuscript <-factor(c.type,
                                 levels= c('Mo','Mac-a','Mac-b','Mac-c','Mac-d','Mac-e','prol.Mac',
                                           'MoDC','cDC1','cDC2','DCx','pDC','prol.DC','MAST','GC','NK','TC',
                                           'UI-M','Fb','dMES1','dMES2','dMES3','prol.MES','eMES','pMES',
                                           'EC','vPC1','vPC2','aPC','prol.PC',
                                           'nmSC','prol.SC','rSC','mSC','UI-SC'))
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
## saveRDS(Final_Seurat, file='Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/post_injury_day_7.RDS')
