library(Seurat)
library(dplyr)
library(future)
plan("multisession", workers = 8)

### This file will load the GEO seurat objects,
## shrink them to make the shiny reasonable,
## calculate cluster markers,
## run cell chat
## output files needed for the shiny

## Load the full PBMC Seurat Dataset. Available for download at GEO Accession = GSE198582
## Shrink it for RShiny
s.obj.blood <- readRDS('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Peripheral_Blood.RDS')
format(object.size(s.obj.blood), units='Gb')
s.obj.blood <- DietSeurat(s.obj.blood, counts=FALSE, data=TRUE, scale.data=FALSE,
                          assays='RNA', dimreducs=c('pca','umap'), graphs=NULL )
format(object.size(s.obj.blood), units='Gb')
saveRDS(s.obj.blood, file='Peripheral_Blood_Diet.RDS')

## Load the Day 0 full Seurat Dataset. Available on GEO
## Shrink it for RShiny
s.obj.0d <- readRDS('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/post_injury_day_0.RDS')
format(object.size(s.obj.0d), units='Gb')
s.obj.0d <- DietSeurat(s.obj.0d, counts=FALSE, data=TRUE, scale.data=FALSE,
                       assays='RNA', dimreducs=c('pca','umap'), graphs=NULL )
format(object.size(s.obj.0d), units='Gb')
saveRDS(s.obj.0d, file='post_injury_day_0_Diet.RDS')

## Load the Day 1 full Seurat Dataset. Available on GEO
## Shrink it for RShiny
s.obj.1d <- readRDS('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/post_injury_day_1.RDS')
format(object.size(s.obj.1d), units='Gb')
s.obj.1d <- DietSeurat(s.obj.1d, counts=FALSE, data=TRUE, scale.data=FALSE,
                       assays='RNA', dimreducs=c('pca','umap'), graphs=NULL )
format(object.size(s.obj.1d), units='Gb')
saveRDS(s.obj.1d, file='post_injury_day_1_Diet.RDS')

## Load the Day 3 full Seurat Dataset. Available on GEO
## Shrink it for RShiny
s.obj.3d <- readRDS('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/post_injury_day_3.RDS')
format(object.size(s.obj.3d), units='Gb')
s.obj.3d <- DietSeurat(s.obj.3d, counts=FALSE, data=TRUE, scale.data=FALSE,
                       assays='RNA', dimreducs=c('pca','umap'), graphs=NULL )
format(object.size(s.obj.3d), units='Gb')
saveRDS(s.obj.3d, file='post_injury_day_3_Diet.RDS')

## Load the Day 7 full Seurat Dataset. Available on GEO
## Shrink it for RShiny
s.obj.7d <- readRDS('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/post_injury_day_7.RDS')
format(object.size(s.obj.7d), units='Gb')
s.obj.7d <- DietSeurat(s.obj.7d, counts=FALSE, data=TRUE, scale.data=FALSE,
                       assays='RNA', dimreducs=c('pca','umap'), graphs=NULL )
format(object.size(s.obj.7d), units='Gb')
saveRDS(s.obj.7d, file='post_injury_day_7_Diet.RDS')

## Load the Location Distal / Injury Seurat Dataset. Available on GEO
## Shrink it for RShiny
s.obj.loc <- readRDS('Y:/CDBLabs-Bioinformatics/Some_Place_Where_Big_Files_Can_Sit/Myeloid_Subset_3Day_Distal_Injury_Site.RDS')
format(object.size(s.obj.loc), units='Gb')
s.obj.loc <- DietSeurat(s.obj.loc, counts=FALSE, data=TRUE, scale.data=FALSE,
                        assays=c('RNA'), dimreducs=c('pca','umap'), graphs=NULL )
format(object.size(s.obj.loc), units='Gb')
saveRDS(s.obj.loc, 'Location_MyeloidOnly_Diet.RDS')

## s.obj.0d <- readRDS("post_injury_day_0_Diet.RDS")
## s.obj.1d <- readRDS("post_injury_day_1_Diet.RDS")
## s.obj.3d <- readRDS("post_injury_day_3_Diet.RDS")
## s.obj.7d <- readRDS("post_injury_day_7_Diet.RDS")

## made this file for now, should go back and do the "observe" technique I used for the PBMC page
all.genes <- unique(c(rownames(s.obj.0d), rownames(s.obj.1d), rownames(s.obj.3d), rownames(s.obj.7d)))
saveRDS(all.genes, file='All_Gene_Names.RDS')

## made this file for now, should go back and do the "observe" technique I used for the PBMC page
saveRDS(rownames(s.obj.loc), file='Location_Gene_Names.RDS')

## Set Default Assay and annotation
DefaultAssay(s.obj.blood) <- 'RNA'
Idents(s.obj.blood) <- 'Manuscript'
DefaultAssay(s.obj.0d) <- 'RNA'
Idents(s.obj.0d) <- 'Cell.Type'
DefaultAssay(s.obj.1d) <- 'RNA'
Idents(s.obj.1d) <- 'Cell.Type'
DefaultAssay(s.obj.3d) <- 'RNA'
Idents(s.obj.3d) <- 'Cell.Type'
DefaultAssay(s.obj.7d) <- 'RNA'
Idents(s.obj.7d) <- 'Cell.Type'
DefaultAssay(s.obj.loc) <- 'RNA'
Idents(s.obj.loc) <- 'CellType'
## Blood data is a stand-alone in iSNAT ..
## There is no overall celltype because it isn't sciatic nerve

#############################
##### Cell Type Markers #####
#############################
## Day 0
TLM <- FindAllMarkers(s.obj.0d, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_top_lvl_0d.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)

## Day 1
TLM <- FindAllMarkers(s.obj.1d, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_top_lvl_1d.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)

## Day 3
TLM <- FindAllMarkers(s.obj.3d, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_top_lvl_3d.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)

## Day 7
TLM <- FindAllMarkers(s.obj.7d, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_top_lvl_7d.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)

## Blood (No Cell Type, not sciatic nerve)
TLM <- FindAllMarkers(s.obj.blood, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Blood.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)

## Location
Markers <- FindAllMarkers(s.obj.loc, assay='RNA', only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25,
                          random.seed=42, verbose=T)
Markers <- Markers %>% select(-p_val) %>% filter(!grepl('^Rpl|^Rps|^Mt-|*Rik', gene)) %>% relocate(gene)
saveRDS(Markers, 'Location_Markers.RDS')

## Distal vs Injury
Distal.pos <- NULL
for(i in 1:length(levels(s.obj.loc))){
    Distal.pos[[i]] <- FindMarkers(s.obj.loc, assay='RNA', only.pos=TRUE, min.pct=0.25,
                                   logfc.threshold=0.25, ransom.seed=42, verbose=T,
                                   group.by='Location', subset.ident=levels(s.obj.loc)[i],
                                   ident.1='Distal', ident.2='Injury')
}
names(Distal.pos) <- levels(s.obj.loc)

## Injury vs Distal
Injury.pos <- NULL
for(i in 1:length(levels(s.obj.loc))){
    Injury.pos[[i]] <- FindMarkers(s.obj.loc, assay='RNA', only.pos=TRUE, min.pct=0.25,
                                   logfc.threshold=0.25, ransom.seed=42, verbose=T,
                                   group.by='Location', subset.ident=levels(s.obj.loc)[i],
                                   ident.1='Injury', ident.2='Distal')
}
names(Injury.pos) <- levels(s.obj.loc)

saveRDS(Distal.pos, file='Distal.DE.RDS')
saveRDS(Injury.pos, file='Injury.DE.RDS')

#############################
##### Immune Cells only #####
#############################
#### Subset Naive to include Immune Cells
#### Find all DE genes relative to this subset gene universe
sub <- subset(s.obj.0d, subset = Manuscript %in% c("epi-Mac", "endo-Mac", "Mac0", "DC", "MAST", "NK",
                                                   "TC", "prol.Mac", "UI1"))
Idents(sub) <- 'Manuscript'

TLM <- FindAllMarkers(sub, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Subset_Immune_0d.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)


#### Subset 1 Day to include Immune Cells
#### Find all DE genes relative to this subset gene universe
sub <- subset(s.obj.1d, subset = Manuscript %in% c("Mo", "Mac-I", "Mac-II", "Mac-III", "Mac-IV", "Mac-V",
                                                   "MoDC", "cDC", "DCx", "pDC", "MAST", "GC", "GC-L",
                                                   "T_NK", "UI-M1"))
Idents(sub) <- 'Manuscript'

TLM <- FindAllMarkers(sub, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Subset_Immune_1d.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)


#### Subset 3 Day to include Immune Cells
#### Find all DE genes relative to this subset gene universe
sub <- subset(s.obj.3d, subset = Manuscript %in% c("Mo", "Mac1", "Mac2", "Mac3", "Mac4", "prol.Mac",
                                                   "MoDC", "cDC", "DCx", "pDC", "MAST", "GC", "T_NK"))
Idents(sub) <- 'Manuscript'

TLM <- FindAllMarkers(sub, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Subset_Immune_3d.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)

#### Subset 7 Day to include Immune Cells
#### Find all DE genes relative to this subset gene universe
sub <- subset(s.obj.7d, subset = Manuscript %in% c("Mo", "Mac-a", "Mac-b", "Mac-c", "Mac-d", "Mac-e",
                                                   "prol.Mac", "MoDC", "cDC1", "cDC2", "DCx", "pDC",
                                                   "prol.DC", "MAST", "GC", "NK", "TC", "TC_NK"))
Idents(sub) <- 'Manuscript'

TLM <- FindAllMarkers(sub, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>%
                        group_by(cluster) %>%
                        top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Subset_Immune_7d.txt'),
            quote=FALSE, sep="\t", row.names=FALSE)

#################################
##### Structural Cells Only #####
#################################
#### Subset Naive to have only "Structure Cells""
#### Find all DE genes relative to this subset gene universe
sub <- subset(s.obj.0d, subset = Cell.Type %in% c("eMES","pMES","prol.MES","Fb"))
Idents(sub) <- 'Manuscript'

TLM <- FindAllMarkers(sub, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Subset_Structure_0d.txt'), quote=FALSE, sep="\t", row.names=FALSE)


#### Subset 1 Day to have only "Structure Cells"
#### Find all DE genes relative to this subset gene universe
sub <- subset(s.obj.1d, subset = Cell.Type %in% c("eMES", "dMES", "prol.MES", "pMES", "Fb"))
Idents(sub) <- 'Manuscript'

TLM <- FindAllMarkers(sub, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Subset_Structure_1d.txt'), quote=FALSE, sep="\t", row.names=FALSE)


#### Subset 3 Day to have only "Structure Cells"
#### Find all DE genes relative to this subset gene universe
sub <- subset(s.obj.3d, subset = Cell.Type %in% c("Fb", "dMES", "prol.MES", "eMES", "pMES"))
Idents(sub) <- 'Manuscript'

TLM <- FindAllMarkers(sub, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Subset_Structure_3d.txt'), quote=FALSE, sep="\t", row.names=FALSE)

#### Subset 7 Day to have only "Structure Cells"
#### Find all DE genes relative to this subset gene universe
sub <- subset(s.obj.7d, subset = Cell.Type %in% c("Fb", "dMES", "prol.MES", "eMES", "pMES"))
Idents(sub) <- 'Manuscript'

TLM <- FindAllMarkers(sub, assay="RNA", only.pos = TRUE, min.pct = 0.25,
                      logfc.threshold = 0.25, random.seed=42, verbose=T)
TLM.df <- as.data.frame(TLM %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC))
write.table(TLM.df, file=paste0('Top_100_Markers_Subset_Structure_7d.txt'), quote=FALSE, sep="\t", row.names=FALSE)


####################
##### CellChat #####
####################
library(CellChat)
options(stringsAsFactors = FALSE)
plan("multisession", workers = 6)

## 0 Day
Idents(s.obj.0d) <- 'Manuscript'
## levels(s.obj.0d)
## Remove proliferating and Unidentified
s.obj.0d.sub <- subset(s.obj.0d, subset=Manuscript %in% c("epi-Mac", "endo-Mac", "Mac0", "DC", "MAST",
                                                          "NK", "TC", "Fb1", "Fb2", "Fb3",
                                                          "eMES1", "eMES2", "eMES3", "pMES1", "pMES2",
                                                          "EC1", "EC2", "EC3", "lyEC", "vSMC", "PC",
                                                          "nmSC", "mSC"))
s.obj.0d.sub$Manuscript <- droplevels(s.obj.0d.sub$Manuscript)
## 1 group per cell type
c.type <- s.obj.0d.sub$Manuscript %>% as.character() %>%
    recode(
        'TC'='T_NK', 'NK'='T_NK',
        'Fb1'='Fb', 'Fb2'='Fb', 'Fb3'='Fb',
        'eMES1'='eMES', 'eMES2'='eMES', 'eMES3'='eMES',
        'pMES1'='pMES', 'pMES2'='pMES',
        'EC1'='EC', 'EC2'='EC', 'EC3'='EC',
        'mSC'='SC', 'nmSC'='SC')
s.obj.0d.sub$CT.Groups <-factor(c.type,
                                levels= c('epi-Mac','endo-Mac','Mac0','DC','MAST','T_NK',
                                          'Fb','eMES','pMES','EC','vSMC','PC','SC'))
Idents(s.obj.0d.sub) <- 'CT.Groups'
cellchat.0d <- createCellChat(s.obj.0d.sub, group.by='CT.Groups')
cellchat.0d@DB <- CellChatDB.mouse
cellchat.0d <- cellchat.0d %>% subsetData() %>% identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>% projectData(adjMatrix=PPI.mouse) %>%
    computeCommunProb(type =  "truncatedMean", trim = 0.25) %>%
    filterCommunication(min.cells = 10) %>% computeCommunProbPathway() %>% aggregateNet()
df.net.0d <- subsetCommunication(cellchat.0d)
pathways.0d <- unique(df.net.0d$pathway_name)

### 1 Day
Idents(s.obj.1d) <- 'Manuscript'
levels(s.obj.1d)
## Remove proliferating and Unidentified
s.obj.1d.sub <- subset(s.obj.1d, subset=Manuscript %in% c("Mo", "Mac-I", "Mac-II", "Mac-III",
                                                          "Mac-IV", "Mac-V", "MoDC", "DCx",
                                                          "cDC", "pDC", "MAST", "GC", "GC-L", "T_NK",
                                                          "eMES", "dMES", "pMES",
                                                          "Fb", "EC", "vSMC_PC",  "SC"))
s.obj.1d.sub$Manuscript <- droplevels(s.obj.1d.sub$Manuscript)
s.obj.1d.sub$CT.Groups <- s.obj.1d.sub$Manuscript
Idents(s.obj.1d.sub) <- 'CT.Groups'
cellchat.1d <- createCellChat(s.obj.1d.sub, group.by='CT.Groups')
cellchat.1d@DB <- CellChatDB.mouse
cellchat.1d <- cellchat.1d %>% subsetData() %>% identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>% projectData(adjMatrix=PPI.mouse) %>%
    computeCommunProb(type =  "truncatedMean", trim = 0.15) %>%
    filterCommunication(min.cells = 10) %>% computeCommunProbPathway() %>% aggregateNet()


### 3 Day
Idents(s.obj.3d) <- 'Manuscript'
levels(s.obj.3d)
## Remove proliferating and Unidentified
s.obj.3d.sub <- subset(s.obj.3d, subset=Manuscript %in% c("Mo", "Mac1", "Mac2", "Mac3", "Mac4",
                                                          "MoDC", "cDC", "DCx", "pDC", "MAST",
                                                          "GC", "T_NK", "Fb", "dMES", "eMES",
                                                          "pMES", "EC1", "EC2", "vSMC_PC1",
                                                          "vSMC_PC2", "nmSC", "rSC"))
s.obj.3d.sub$Manuscript <- droplevels(s.obj.3d.sub$Manuscript)
## 1 group per cell type
c.type <- s.obj.3d.sub$Manuscript %>% as.character() %>%
    recode(
        'vSMC_PC1'='vSMC_PC', 'vSMC_PC2'='vSMC_PC',
        'EC1'='EC', 'EC2'='EC',
        'nmSC'='SC', 'rSC'='SC')
s.obj.3d.sub$CT.Groups <-factor(c.type,
                                levels= c("Mo", "Mac1", "Mac2", "Mac3", "Mac4", "MoDC", "cDC",
                                          "DCx", "pDC", "MAST", "GC", "T_NK", "Fb", "dMES",
                                          "eMES", "pMES", "EC", "vSMC_PC", "SC"))
Idents(s.obj.3d.sub) <- 'CT.Groups'
cellchat.3d <- createCellChat(s.obj.3d.sub, group.by='CT.Groups')
cellchat.3d@DB <- CellChatDB.mouse
cellchat.3d <- cellchat.3d %>% subsetData() %>% identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>% projectData(adjMatrix=PPI.mouse) %>%
    computeCommunProb(type =  "truncatedMean", trim = 0.15) %>%
    filterCommunication(min.cells = 10) %>% computeCommunProbPathway() %>% aggregateNet()


## 7 Day
Idents(s.obj.7d) <- 'Manuscript'
levels(s.obj.7d)
## Remove proliferating and Unidentified
s.obj.7d.sub <- subset(s.obj.7d, subset=Manuscript %in% c("Mo", "Mac-a", "Mac-b", "Mac-c", "Mac-d",
                                                          "Mac-e", "MoDC", "cDC1", "cDC2", "DCx",
                                                          "pDC", "MAST", "GC", "NK", "TC", "Fb",
                                                          "dMES1", "dMES2", "dMES3", "eMES", "pMES",
                                                          "EC", "vPC1", "vPC2", "aPC", "nmSC", "rSC",
                                                          "mSC"))
s.obj.7d.sub$Manuscript <- droplevels(s.obj.7d.sub$Manuscript)
## 1 group per cell type
c.type <- s.obj.7d.sub$Manuscript %>% as.character() %>%
    recode(
        'cDC1'='cDC', 'cDC2'='cDC',
        'dMES1'='dMES', 'dMES2'='dMES', 'dMES3'='dMES',
        'vPC1'='vSMC_PC', 'vPC2'='vSMC_PC',
        'nmSC'='SC', 'rSC'='SC', 'mSC'='SC')
s.obj.7d.sub$CT.Groups <-factor(c.type,
                                levels= c("Mo", "Mac-a", "Mac-b", "Mac-c", "Mac-d",
                                          "Mac-e", "MoDC", "cDC", "DCx",
                                          "pDC", "MAST", "GC", "NK", "TC", "Fb",
                                          "dMES", "eMES", "pMES",
                                          "EC", "vSMC_PC", "aPC", "SC"))
Idents(s.obj.7d.sub) <- 'CT.Groups'
cellchat.7d <- createCellChat(s.obj.7d.sub, group.by='CT.Groups')
cellchat.7d@DB <- CellChatDB.mouse
cellchat.7d <- cellchat.7d %>% subsetData() %>% identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>% projectData(adjMatrix=PPI.mouse) %>%
    computeCommunProb(type =  "truncatedMean", trim = 0.15) %>%
    filterCommunication(min.cells = 10) %>% computeCommunProbPathway() %>% aggregateNet()

saveRDS(cellchat.0d, file='CellChat_day_0.RDS')
saveRDS(cellchat.1d, file='CellChat_day_1.RDS')
saveRDS(cellchat.3d, file='CellChat_day_3.RDS')
saveRDS(cellchat.7d, file='CellChat_day_7.RDS')

### Create some variables I need in the Shiny ui
df.net.0d <- subsetCommunication(cellchat.0d)
pathways.0d <- unique(df.net.0d$pathway_name)
df.net.1d <- subsetCommunication(cellchat.1d)
pathways.1d <- unique(df.net.1d$pathway_name)
df.net.3d <- subsetCommunication(cellchat.3d)
pathways.3d <- unique(df.net.3d$pathway_name)
df.net.7d <- subsetCommunication(cellchat.7d)
pathways.7d <- unique(df.net.7d$pathway_name)

possible.source.clusters.0d <- levels(df.net.0d$source)
possible.source.clusters.1d <- levels(df.net.1d$source)
possible.source.clusters.3d <- levels(df.net.3d$source)
possible.source.clusters.7d <- levels(df.net.7d$source)
possible.pathways <- unique(c(pathways.0d, pathways.1d, pathways.3d, pathways.7d))
possible.pathways <- possible.pathways[order(possible.pathways)]

temp <- list("possible.source.clusters.0d" = possible.source.clusters.0d,
             "possible.source.clusters.1d" = possible.source.clusters.1d,
             "possible.source.clusters.3d" = possible.source.clusters.3d,
             "possible.source.clusters.7d" = possible.source.clusters.7d,
             "possible.pathways" = possible.pathways)
saveRDS(temp, 'CellChat_Choices.RDS')
