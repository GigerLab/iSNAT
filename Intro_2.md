# Specifics about the Shiny Application
1. [Expression Analysis](#Exp_Analysis)  
2. [Two Gene](#TwoGene")  
3. [Immune Cells](#immune_cells)  
4. [Structural Cells](#structure)  
5. [Integrated Immune Cells](#integrated)  
6. [Trajectory](#trajectory)  
7. [Cell Chat](#CellChat)  
8. [Ligand-Receptor](#LR)  
9. [Spatial Distribution](#spatial)  
10. [PBMC](#PBMC)
11. [General Tips for Using the Application](#general)

## Expression Analysis <a name="Exp_Analysis"></a>

This function allows gene expression analysis in the naïve nerve and at 1 day, 3 days, and 7 days post-injury. Enter the MGI symbol of your gene of interest in the autofill text box. If the gene is in the dataset, it will appear, click to view feature plots.

As an example, Arg1 (for arginase-1) is shown in the text box by default. Feature plots show that Arg1 is expressed by few cells in the naïve nerve, strongly upregulated in many immune cells at 1 Day and 3 Days Post-Injury. However, at 7 Days Post-Injury, expression has greatly declined.

Under the feature plots there is a pulldown menu where you can select any of the cell clusters to see the top genes (Cluster Positive Markers). Markers for the selected cell cluster are shown in naïve nerve (Naïve markers), 1 Day, 3 Days, and 7 Days post-injury (time points are tabs in the box).

Inside the box you will find the top markers determined using the Wilcoxon Rank Sum test.

Displayed genes are positively expressed in the cluster relative to all other cells.

Stats in the box include:
- avg_log2FC = log 2 fold change of cells in cluster compared to cells out of cluster
- pct.1 = percent of cells in the cluster expressing the gene
- pct.2 = percent of cells outside the cluster expressing the gene
- p_val_adj = The FDR adjusted P-value from the Wilcoxon Rank Sum Test

There is a second box (Cluster Positive Markers Box 2) identical to the first for side-by-side comparison of gene expression at different time points or in different cell clusters.

## Two Gene <a name="TwoGene"></a>
This tab allows you to compare expression of two genes of interest. 

Enter names of Gene 1 and Gene 2 and click “GO” to view feature plots of Gene 1 in red and Gene 2 in green (be patient as plots render). Cells co-expressing both genes are labeled yellow.  Note that each individual gene is scaled so its expression in the dataset is between 0 and 10 with 10 being the cells that most strongly express the gene. A yellow cell therefor expresses both genes at high levels respective to expression throughout the dataset. Thus, this feature allows for comparison of genes with different absolute expression levels.  

The scatter plot shows the percentile of cell types (color coded) that co-express the two selected genes.

## Immune Cells <a name="immune_cells"></a>
This tab allows you to assess gene expression in immune cell subpopulations. Enter the gene name and obtain feature plots at each of the four different time points.

Under the feature plots there is a pulldown menu where you can select any immune cell clusters to see the top cluster enriched genes at different time points.

Inside the box you will find the top markers determined using the Wilcoxon Rank Sum test.

Resulting genes are positively expressed in the cluster relative to all other immune cells

Stats in the box include:
- avg_log2FC = log 2 fold change of cells in cluster compared to cells out of cluster
- pct.1 = percent of cells in the cluster expressing the gene
- pct.2 = percent of cells outside the cluster expressing the gene
- p_val_adj = The FDR adjusted P-value from the Wilcoxon Rank Sum Test

## Structural Cells <a name="structure"></a>
This tab allows you to assess gene expression in structural cells, including Fb (fibroblasts), perineurial mesenchymal cells (pMES), endoneurial mesenchymal cells (eMES), and differentiating mesenchymal cells (dMES).

Enter the gene name to obtain feature plots at each of the four different time points.

Stats in the box include:
- avg_log2FC = log 2 fold change of cells in cluster compared to cells out of cluster
- pct.1 = percent of cells in the cluster expressing the gene
- pct.2 = percent of cells outside the cluster expressing the gene
- p_val_adj = The FDR adjusted P-value from the Wilcoxon Rank Sum Test

## Integrated Immune Cells <a name="integrated"></a>
For longitudinal studies of gene expression, we used computational methods to isolate and integrate all immune cell scRNAseq datasets. Click tab at top of page to load integrated immune dataset  

The top raw shows color coded UMAP plots of naïve sciatic nerve, PBMC in cardiac blood, Day 1, Day 3, and Day 7 injured sciatic nerve. For longitudinal analysis of gene expression, enter a gene name to view feature plots.

The integrated UMAP plots show numbered and color-coded Seurat clusters for analysis of cluster enriched genes (markers).  Select the Seurat cluster of interest (0 – 22) to display cluster enriched genes. 

Stats in the box include:
- avg_log2FC = log 2 fold change of cells in cluster compared to cells out of cluster
- pct.1 = percent of cells in the cluster expressing the gene
- pct.2 = percent of cells outside the cluster expressing the gene
- p_val_adj = The FDR adjusted P-value from the Wilcoxon Rank Sum Test

## Trajectory <a name="trajectory"></a>
Click tab at top of **Integrated Immune Cells** tab to load integration data  
Click tab at top of trajectory page to load Trajectory data  

Trajectory analysis was carried out using slingshot from the first 5 Principal Components of the integrated Immune dataset. 

From the "Trajectory" pull-down menu select the trajectory you are interested in. 
From the "Label Clusters By" select whether to label plots by the Macro Trajectory_CellTypes or by the integrated Suerat Clusters.  

From left to right you will see:  
a PCA plot with cells labeled with the selected labels,  
a UMAP plot with cells labeled with the selected labels,  
the PCA plot with cells of the trajectory colored by pseudotime (from slingshot),  
the UMAP plot with cells of the trajectory colored by psuedotime.  

Under these plots is a table showing genes and their rank of importance in predicted pseudotime based on the random forests machine learning technique.  
The table can be sorted by either trajectories column; The search field will allow you to look for your favorite gene.  

Entering a gene in the "Gene:" field will generate feature plots of that gene using the PCA and UMAP. 

## Cell Chat <a name="CellChat"></a>
Click tab at top of page to load CellChat data  

Select a "source" cell type of interest at a specific time point to view with which cells in the nerve an interaction may occur. The ligand and receptor pairs through which interactions may occur are listed. The CellChat derived probability of the interaction is shown, and the pathway name is listed.

## Ligand-Receptor <a name="LR"></a>
Click tab at top of **CellChat page** load CellChat data  

CellChat software was used for identification of ligand-receptor pairs in the naïve and injured nerve. CellChat calculates the probability for the listed cell-cell communications to occur through the selected signaling pathway. The circular plot shows sender cells and predicts receiver cells. The bar graph shows the contribution of specific ligand-receptor pairs to the signaling pathway under investigation.

## Spatial distribution <a name="spatial"></a>

At 3 days following sciatic nerve crush, innate immune cells (anti-CD11b) were purified from a 3 mm nerve segment that contains the injury site, and from the same mice, from a 3 mm nerve segment distal to injury site. Cells were isolated and subjected to scRNA-seq. Feature plots of "injury" versus "distal" reveal if a gene is preferentially expressed by Mac subpopulations at the nerve crush site versus the distal nerve.

## PBMC <a name="PBMC"></a>

Click tab at top of page to load PBMC data.  

Enter gene name to view feature plots. 

## General Tips for Using the Application <a name="general"></a>

Note that the number of cells is high for our humble resources. For this reason much of the data is loaded only when necessary. 
You will often see a "Load Data" button in the top-left of the page, please click the link to load the data for that tab. 

If you have a smaller screen you can collapse the left panel with the three lines in the blue header next to the title iSNAT.  

MGI nomenclature is used for gene symbols. If your gene does not show in the autofill box then it is not in the dataset (at least with the name you tried). 


