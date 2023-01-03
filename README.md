# This Repository contains the code to locally run the **Injured Sciatic Nerve Atlas (iSNAT) described in this [eLife publication](https://doi.org/10.7554/eLife.80881 "The Injured Sciatic Nerve Atlas (iSNAT), insights into the cellular and molecular basis of neural tissue degeneration and regeneration")**  

![Workflow for peripheral nervous tissue analysis](cartoon.png)

# Order of events needed to run application

## Pre-Process
The Pre-Process folder has the scripts used to take CellRanger output and make the Seruat Objects that are uploaded to GEO: GSE198582. You can download these objects and skip this step

## Markers_and_RShiny_Files.r
This file will take the Seurat Objects from GEO and shrink them for use in the shiny application (if you know how to make them smaller please share). It will also produce the files needed for the application including top cluster markers (sorry, all flat files) and CellChat.

## Immune_Integration_SCTransform_V2.r
This file will also start at the full Seurat objects from GEO but will subset the immune cells and integrate all nerve time points, and PBMC, into a single object. It will also conduct the trajectory analysis and output all files needed for the shiny integration and trajectory pages (Shiny can run without these). Note that this step was done with the development version of Seurat and a more current version of R, see Session_Logs for exact versions.  

## app.r
After the previous two files are run, app.r should start iSNAT on your computer.
