library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(glmGamPoi)

### This file will do the Immune integration and output files needed for the integration section of the Shiny
### These 5 RDS files can be downloaded at GEO Accession = GSE198582
wd <- getwd()
setwd("Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit")
blood <- readRDS('Peripheral_Blood.RDS')
day0 <- readRDS('post_injury_day_0.RDS')
day1 <- readRDS('post_injury_day_1.RDS')
day3 <- readRDS('post_injury_day_3.RDS')
day7 <- readRDS('post_injury_day_7.RDS')
setwd(wd)

## for DE store genes that are in every dataset
genes.in.all.times <- Reduce(intersect, list(rownames(blood), rownames(day0),
                                             rownames(day1), rownames(day3), rownames(day7)))
##

myeloid.cells.blood <- colnames(blood)[which(blood$Manuscript %in%
                                             c("TC1", "TC2", "TC3", "TC4", "TC5", "NK",
                                               "iGC", "mGC", "Mo_1", "Mo_2", "Mac_0",
                                               "DC", "Mast"))]

macro.cell.types.to.keep <- c("GC", "Mo", "Mac", "MoDC", "cDC", "pDC", "DCx", "Mast", "NK", "TC")
myeloid.cells.0d <- colnames(day0)[which(day0$Cell.Type %in% macro.cell.types.to.keep)]
myeloid.cells.1d <- colnames(day1)[which(day1$Cell.Type %in% macro.cell.types.to.keep)]
myeloid.cells.3d <- colnames(day3)[which(day3$Cell.Type %in% macro.cell.types.to.keep)]
myeloid.cells.7d <- colnames(day7)[which(day7$Cell.Type %in% macro.cell.types.to.keep)]

DefaultAssay(blood) <- 'RNA'
DefaultAssay(day0) <- 'RNA'
DefaultAssay(day1) <- 'RNA'
DefaultAssay(day3) <- 'RNA'
DefaultAssay(day7) <- 'RNA'
blood <- subset(blood, cells=myeloid.cells.blood)
blood.v2 <- SCTransform(blood, vst.flavor = "v2", vars.to.regress='percent.mt', verbose=T)
day0 <- subset(day0, cells=myeloid.cells.0d)
day0.v2 <- SCTransform(day0, vst.flavor = "v2", vars.to.regress='percent.mt', verbose=T)
day1 <- subset(day1, cells=myeloid.cells.1d)
day1.v2 <- SCTransform(day1, vst.flavor = "v2", vars.to.regress='percent.mt', verbose=T)
day3 <- subset(day3, cells=myeloid.cells.3d)
day3.v2 <- SCTransform(day3, vst.flavor = "v2", vars.to.regress='percent.mt', verbose=T)
day7 <- subset(day7, cells=myeloid.cells.7d)
day7.v2 <- SCTransform(day7, vst.flavor = "v2", vars.to.regress='percent.mt', verbose=T)
rm(blood, day0, day1, day3, day7)

## Integrate together
all.data.list <- list(Blood = blood.v2, Niave = day0.v2, Day1 = day1.v2, Day3 = day3.v2, Day7 = day7.v2 )
rm(blood.v2, day1.v2, day3.v2, day7.v2)
features <- SelectIntegrationFeatures(all.data.list, nfeatures = 3000)
all.data.list <- PrepSCTIntegration(all.data.list, anchor.features = features)

anchors <- FindIntegrationAnchors(all.data.list, normalization.method = "SCT",
                                  anchor.features = features)
combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(anchors, all.data.list)

combined <- RunPCA(combined, verbose = FALSE)
ElbowPlot(combined, ndims=50)
combined <- combined %>%
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:20) %>%
    FindClusters(resolution = 0.5)

combined$Time <- combined$orig.ident %>%
    recode( "PBMC_Naive_First"='Blood', "PBMC_Naive_Second"='Blood', "Giger_UCLA-4884_03"='Naive',
           "Giger_UCLA-3088_01" = 'Day1',
           "Giger_UCLA-3088_02" = 'Day1',  "Giger_UCLA-3088_03" = 'Day1', "Giger_UCLA-3088_04" = 'Day1',
           "Giger_UCLA-4884_05" = 'Day1', "Giger_sciatic_Batch1" = 'Day3',
           "Giger_sciatic_Batch1.2" = 'Day3',
           "Giger_sciatic_Batch2" = 'Day3', "MPZ_myelin_3d" = 'Day3',  "WT_NoDrug_1" = 'Day3',
           "WT_NoDrug_2" = 'Day3', "Giger_UCLA-24_13" = 'Day7',  "Giger_UCLA-24_14" = 'Day7',
           "MPZ_myelin_7d" = 'Day7', "MPZ_myelin_depleted_7d" = 'Day7', "Giger_Sciatic_1" = 'Day7',
           "Giger_Sciatic_2" = 'Day7', "Giger_Sciatic_3" = 'Day7', "Giger_Sciatic_4" = 'Day7')
combined$Time <- factor(combined$Time, levels=c('Blood', 'Naive', 'Day1', 'Day3', 'Day7'))

DefaultAssay(combined) <- 'RNA'
combined <- NormalizeData(combined)

## Find markers by seurat cluster for iSNAT
TLM <- FindAllMarkers(combined, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
saveRDS(TLM, 'Myeloid_Integration_Seurat_Cluster_Markers.RDS')


## meta.table <- as_tibble(combined@meta.data)
## tally <- meta.table %>%
##     group_by(seurat_clusters, Manuscript) %>%
##     tally() %>%
##     group_by(seurat_clusters) %>%
##     mutate(per=round(n/sum(n),4)) %>%
##     group_by(Manuscript) %>%
##     mutate(per2=round(n/sum(n),4)) %>%
##     arrange(seurat_clusters, desc(n))
## write.table(data.frame(tally), "Mapping_SeuratClusters_to_ManuscriptLabels.csv",
##             quote=F, sep=',', row.names=F)

## Macro labels for Trajectory
temp1 <- DimPlot(combined, group.by='seurat_clusters', label=T) + NoLegend()
temp2 <- DimPlot(combined, group.by='Manuscript', label=T) + NoLegend()
c.type <- combined$seurat_clusters %>% as.character() %>%
    recode(
        '0'='Mac',
        '1'='Mac',
        '2'='Mo',
        '3'='GC',
        '4'='MoDC',
        '5'='TC',
        '6'='Mac',
        '7'='TC',
        '8'='Mac',
        '9'='NK',
        '10'='TC',
        '11'='Mac',
        '12'='GC',
        '13'='pDC',
        '14'='TC',
        '15'='DCx',
        '16'='cDC',
        '17'='Mac',
        '18'='prol.Mac',
        '19'='Mast')
combined$Trajectory_CellTypes <-factor(c.type,
                                       levels= c('GC','Mo','Mac','prol.Mac','MoDC','cDC',
                                                 'pDC','DCx','Mast','NK','TC'))

### Slingshot
library(slingshot)
library(viridis)

sds <- slingshot(Embeddings(combined, "pca")[,1:4],
                 clusterLabels = combined$seurat_clusters,
                 stretch = 0, start.clus = '2', allow.breaks=F)
saveRDS(sds, file='sds.RDS')

## random forest
top_hvg <- combined@assays$integrated@var.features

## Subset highly variable
dat_use <- t(as.matrix(GetAssayData(combined, slot = "data"))[top_hvg,])
marker.list <- list()

## Genes important to pseudotime
library(ranger)
library(tidymodels)
for(i in c(6,9)){
    curve.to.use <- i  ## change to curve you want [9 or 6]
    dat_use_df <- cbind(slingPseudotime(true.sds)[,curve.to.use], dat_use)
    ## First column is pseudotime
    colnames(dat_use_df)[1] <- "pseudotime"
    ## Remove the cells not in this curve
    dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

    ## Split into training and validataion step using rsample
    set.seed(42)
    dat_split <- initial_split(dat_use_df)
    dat_train <- training(dat_split)
    dat_val <- testing(dat_split)

    ## ranger doesn't like colnames(dat_train)
    ## must be leading numbers or something wrong with gene name
    ## Fix them with R base
    colnames(dat_train) <- make.names(colnames(dat_train),unique=T)
    colnames(dat_val) <- make.names(colnames(dat_val),unique=T)

    ## Use tidymodels for the random forests
    model <- rand_forest(mtry=200, trees=1400, min_n=15, mode="regression") %>%
        set_engine("ranger", importance="impurity", num.threads=6, seed=42) %>%
        fit(pseudotime ~ ., data=dat_train)

    ## Look at some model stats
    val_results <- dat_val %>%
        mutate(estimate = predict(model, .[,-1]) %>% pull()) %>%
        dplyr::select(truth = pseudotime, estimate)
    metrics(data = val_results, truth, estimate)

    ## Get top genes with their "variable importance"
    var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
    var_imp_rank <- 1:length(var_imp)
    names(var_imp_rank) <- names(var_imp)
    if(i==6){ saveRDS(var_imp_rank, file='Variable_Importance_Mo_Mac.RDS')}
    if(i==9){ saveRDS(var_imp_rank, file='Variable_Importance_Mo_DC.RDS')}
}

## Shrink File for RShiny
DefaultAssay(combined) <- 'RNA'
format(object.size(combined), units='Gb')
combined_Diet <- DietSeurat(combined, counts=FALSE, data=TRUE, scale.data=FALSE,
                            assays=c('RNA'), dimreducs=c('pca','umap'), graphs=NULL )
combined_Diet@tools <- list()
## deleted pca even though I told it not to
combined_Diet@reductions$pca <- combined@reductions$pca
format(object.size(combined_Diet), units='Gb')
saveRDS(combined_Diet, file='Myeloid_Integration_Diet.RDS')

library(sessioninfo)
session_info(to_file = "Immune_Integration_session.log")

