wd <- getwd()
library(Seurat)
library(dplyr)
library(ggplot2)

RG_4884_03 <- Read10X('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_0/4884-RG-3/filtered_feature_bc_matrix')
RG_4884_04 <- Read10X('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Day_0/4884-RG-4/filtered_feature_bc_matrix')
s_4884_01_obj <- CreateSeuratObject(counts = RG_4884_03, project = "Giger_UCLA-4884_03")
s_4884_02_obj <- CreateSeuratObject(counts = RG_4884_04, project = "Giger_UCLA-4884_04")

combined.obj <- merge(s_4884_01_obj, y = s_4884_02_obj,
                      add.cell.ids = c("IM", "O"), project = "UCLA-4884")

## SCTransform
ncol(combined.obj)
X <- PercentageFeatureSet(combined.obj, pattern = "mt-", col.name="percent.mt")
X <- subset(X, subset=nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)
X <- SCTransform(X,vars.to.regress='percent.mt',verbose=T)
Final_Seurat <- X
rm(list=ls(pattern="^c|^X|^s|^R"))

## Stats for Table
length(which(Final_Seurat$orig.ident=="Giger_UCLA-4884_03"))
length(which(Final_Seurat$orig.ident=="Giger_UCLA-4884_04"))
median(Final_Seurat$nCount_RNA[which(Final_Seurat$orig.ident=="Giger_UCLA-4884_03")])
median(Final_Seurat$nCount_RNA[which(Final_Seurat$orig.ident=="Giger_UCLA-4884_04")])
median(Final_Seurat$nFeature_RNA[which(Final_Seurat$orig.ident=="Giger_UCLA-4884_03")])
median(Final_Seurat$nFeature_RNA[which(Final_Seurat$orig.ident=="Giger_UCLA-4884_04")])

Final_Seurat <- Final_Seurat %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, return.model=TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20) %>%
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
        '0'='eMES',
        '1'='Fb',
        '2'='Mac',
        '3'='Fb',
        '4'='SC',
        '5'='eMES',
        '6'='pMES',
        '7'='EC',
        '8'='MoDC',
        '9'='vSMC_PC',
        '10'='pMES',
        '11'='vSMC_PC',
        '12'='EC',
        '13'='Fb',
        '14'='EC',
        '15'='Mac',
        '16'='eMES',
        '17'='eMES',
        '18'='Mast',
        '19'='TC',
        '20'='prol.MES',
        '21'='TC',
        '22'='SC')
Final_Seurat$Cell.Type <- c.type

## TC and NK are clustering together
## split out NK as those that express either Ncr1 or Klra7
temp <- subset(Final_Seurat, subset = Cell.Type == "TC" & (Ncr1 > 0 | Klra7 > 0))
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'NK'

Final_Seurat$Cell.Type <-factor(c.type,
                                levels= c('GC','Mo','Mac','MoDC','cDC','pDC','DCx','Mast','NK','TC',
                                          'eMES','pMES','dMES','prol.MES','Fb',
                                          'SC','EC','vSMC_PC'))

p1 <- DimPlot(Final_Seurat, reduction = "umap", label = TRUE, pt.size=.25, group.by='Cell.Type') +
    theme(legend.position = "none")

## Cluster Numbers in order for colors to roughly match in Manuscript,
## i.e. change seurat cluster number order to Cell Type number order
## Proliferation cells in cluster 20 are split in the UMAP with half being Mes and half being Mac ...
## Make a new cluster number for them using the UMAP coordinates
c.type <- as.character(Final_Seurat$seurat_clusters)
temp <- subset(Final_Seurat, subset = seurat_clusters == "20" & UMAP_1 < -9)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '200' #prol.mac

## FeaturePlot(Final_Seurat, features='Fn1') + geom_hline(yintercept =2) ## Mac0
## FeaturePlot(Final_Seurat, features='Cx3cr1') + geom_hline(yintercept =2) ## endo-mac
## FeaturePlot(Final_Seurat, features='Napsa') + geom_hline(yintercept =2) ## DC
temp <- subset(Final_Seurat, subset = seurat_clusters == "8" & UMAP_2 < 2)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '201' #endo

temp <- subset(Final_Seurat, subset = seurat_clusters == "8" & UMAP_2 > 2 & Napsa > 0)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '202' #DC

temp <- subset(Final_Seurat, subset = seurat_clusters == "8" & UMAP_2 > 2 & Fn1 > 0)
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- '203' #Mac0

c.type <- c.type  %>%
    recode('0'='11', '1'='8', '2'='2', '3'='9', '4'='23', '5'='12', '6'='14', '7'='17', '8'='7',
           '9'='21', '10'='15', '11'='22', '12'='18', '13'='10', '14'='19', '15'='1', '16'='20',
           '17'='13', '18'='4', '19'='5', '20'='16', '21'='7', '22'='24', '200'='6', '201'='2',
           '202'='3', '203'='0')

Final_Seurat$Clus_Num_Vis <-factor(c.type,
                                   levels=c('0','1','2','3','4','5','6','7','8','9','10',
                                            '11','12','13','14','15','16','17','18',
                                            '19','20','21','22','23','24'))
### Unique ID for each cluster
c.type <- Final_Seurat$Clus_Num_Vis %>% as.character() %>%
    recode(
        '0'='Mac0',
        '1'='epi-Mac',
        '2'='endo-Mac',
        '3'='DC',
        '4'='MAST',
        '5'='TC',
        '6'='prol.Mac',
        '7'='UI1',
        '8'='Fb1',
        '9'='Fb2',
        '10'='Fb3',
        '11'='eMES1',
        '12'='eMES2',
        '13'='eMES3',
        '14'='pMES1',
        '15'='pMES2',
        '16'='prol.MES',
        '17'='EC1',
        '18'='EC2',
        '19'='EC3',
        '20'='LyEC',
        '21'='vSMC',
        '22'='PC',
        '23'='nmSC',
        '24'='mSC')
temp <- subset(Final_Seurat, subset = Clus_Num_Vis == "5" & (Ncr1 > 0 | Klra7 > 0))
ix <- which(colnames(Final_Seurat) %in% colnames(temp))
c.type[ix] <- 'NK'

Final_Seurat$Manuscript <-factor(c.type,
                                 levels= c('Mac0','epi-Mac','endo-Mac','DC','MAST','NK','TC','prol.Mac','UI1',
                                           'Fb1','Fb2','Fb3','eMES1','eMES2','eMES3','pMES1','pMES2',
                                           'prol.MES','EC1','EC2','EC3','LyEC','vSMC','PC','nmSC','mSC'))
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
## saveRDS(Final_Seurat, file='Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/post_injury_day_0.RDS')
