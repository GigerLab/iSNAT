wd <- getwd()
library(Seurat)
library(dplyr)
library(ggplot2)

samples <- paste0("4361-RG-",1:4)
s.1.obj <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Distal_Injury/4361-RG-1/filtered_feature_bc_matrix") %>% CreateSeuratObject(project=samples[1])
s.2.obj <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Distal_Injury/4361-RG-2/filtered_feature_bc_matrix") %>% CreateSeuratObject(project=samples[2])
s.3.obj <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Distal_Injury/4361-RG-3/filtered_feature_bc_matrix") %>% CreateSeuratObject(project=samples[3])
s.4.obj <- Read10X("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Giger_2022_CellRangerFiles/Distal_Injury/4361-RG-4/filtered_feature_bc_matrix") %>% CreateSeuratObject(project=samples[4])

obj.list <- lapply(ls(pattern=".obj$"), function(X){get(X)})

## SCTransform
unlist(lapply(obj.list, function(X){ncol(X)}))
obj.list <- lapply(obj.list, function(X){
    X <- PercentageFeatureSet(X, pattern = "mt-", col.name="percent.mt")
    X <- subset(X, subset=nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)
    X <- SCTransform(X,vars.to.regress='percent.mt',verbose=T)
})
unlist(lapply(obj.list, function(X){ncol(X)}))

## CCA Integration
obj.features <- SelectIntegrationFeatures(obj.list, nfeatures=3000)
obj.list <- PrepSCTIntegration(object.list=obj.list, anchor.features=obj.features)
anchors <- FindIntegrationAnchors(object.list=obj.list, normalization.method="SCT",
                                  anchor.features=obj.features)
Combined <- IntegrateData(anchorset=anchors, normalization.method="SCT")

## Stats for Table
median(Combined$nCount_RNA[which(Combined$orig.ident %in% c('4361-RG-1', '4361-RG-2'))])
median(Combined$nCount_RNA[which(Combined$orig.ident %in% c('4361-RG-3', '4361-RG-4'))])
median(Combined$nFeature_RNA[which(Combined$orig.ident %in% c('4361-RG-1', '4361-RG-2'))])
median(Combined$nFeature_RNA[which(Combined$orig.ident %in% c('4361-RG-3', '4361-RG-4'))])

Combined <- Combined %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, return.model=TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20) %>%
    FindClusters(resolution = 0.5)

## Limiting to only the Myeloid cells
cells <- colnames(Combined)[which(!Combined$seurat_clusters %in% c(8,9,10,16))]

myeloid <- subset(Combined, cells=cells)
DefaultAssay(myeloid) <- 'RNA'
myeloid <- SCTransform(myeloid, vars.to.regress='percent.mt', verbose=T)
myeloid <- myeloid %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, return.model=TRUE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20) %>%
    FindClusters(resolution = 0.5)

DefaultAssay(myeloid) <- 'RNA'
myeloid <- NormalizeData(myeloid)

c.type <- myeloid$orig.ident %>% as.character() %>%
    recode('4361-RG-1'='Injury', '4361-RG-2'='Injury', '4361-RG-3'='Distal', '4361-RG-4'='Distal')
myeloid$Location <- factor(c.type, levels=c('Injury','Distal'))

### For Shiny
c.type <- myeloid$seurat_clusters %>% as.character() %>%
    recode(
        '0'='Mac3',
        '1'='Mac4',
        '2'='MoDC',
        '3'='Mac1',
        '4'='Mo',
        '5'='UI-M1',
        '6'='GC',
        '7'='Mac1',
        '8'='prol.Mac',
        '9'='cDC',
        '10'='UI-M2',
        '11'='NK',
        '12'='DCx',
        '13'='epiMac',
        '14'='TC')
myeloid$CellType <- c.type
myeloid$CellType <-factor(c.type,
                          levels= c('Mo', 'Mac1', 'Mac3', 'Mac4', 'epiMac', 'prol.Mac', 'GC',
                                    'MoDC', 'cDC', 'DCx', 'TC', 'NK', 'UI-M1', 'UI-M2'))

## Cluster Numbers in order for colors to roughly match in Manuscript,
## i.e. change seurat cluster number order to Cell Type number order
c.type <- as.character(myeloid$seurat_clusters)
c.type <- c.type  %>%
    recode('0'='2', '1'='3', '2'='7', '3'='1', '4'='0', '5'='12', '6'='6', '7'='1', '8'='5',
           '9'='8', '10'='13', '11'='11', '12'='9', '13'='4', '14'='10')
myeloid$Clus_Num_Vis <-factor(c.type,
                              levels=c('0','1','2','3','4','5','6','7','8','9','10',
                                       '11','12','13'))

## ## Save full object - Which will be availble in GEO
##  saveRDS(myeloid, file='Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Myeloid_Subset_3Day_Distal_Injury_Site.RDS')

