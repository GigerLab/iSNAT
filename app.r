library(CellChat)
library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(shinydashboard)
library(stringr)
library(tableHTML)
library(waiter)
library(shinyjs)
library(markdown)
library(DT)
library(slingshot)
library(viridis)
library(shinycssloaders)

waiting_screen <- tagList(
    spin_loaders(12),
    h4("Loading iSNAT Datasets...")
)

### Clusters to use in pulldown
clusters <- c("GC", "Mo", "Mac", "MoDC", "cDC", "pDC", "DCx", "Mast", "NK", "TC",
              "eMES", "pMES", "dMES", "prol.MES", "Fb", "SC", "EC", "vSMC_PC")

top.0d <- read.table("Top_100_Markers_top_lvl_0d.txt", header=T)
top.1d <- read.table("Top_100_Markers_top_lvl_1d.txt", header=T)
top.3d <- read.table("Top_100_Markers_top_lvl_3d.txt", header=T)
top.7d <- read.table("Top_100_Markers_top_lvl_7d.txt", header=T)
top.loc <- readRDS('Location_Markers.RDS')
top.int <- readRDS('Myeloid_Integration_Seurat_Cluster_Markers.RDS')
top.blood <- read.table("Top_100_Markers_Blood.txt", header=T)

distal.DE <- readRDS('Distal.DE.RDS')
injury.DE <- readRDS('Injury.DE.RDS')

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

### Going to swap TC and SC to show differences between TC and NK
gg_cols <- c(GC = "#F8766D", Mo = "#E88526", Mac = "#D39200", MoDC = "#B79F00", cDC = "#93AA00",
             pDC = "#5EB300", DCx = "#00BA38", Mast = "#00BF74", NK = "#00C19F", TC = "#F564E3",
             eMES = "#00B9E3", pMES = "#00ADFA", dMES = "#619CFF", prol.MES = "#AE87FF",
             Fb = "#DB72FB", SC = "#00BFC4", EC = "#FF61C3", vSMC_PC = "#FF699C")

v.plot <- function(m){
    ggplot(m, aes(group,counts)) +
        geom_violin(scale='width', aes(fill=factor(group))) +
        scale_fill_discrete(drop=FALSE) +
        scale_x_discrete(drop=FALSE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.title.x = element_blank()) +
        geom_jitter(width=0.2,size=.5) +
        geom_boxplot(width=.2,outlier.shape=NA) +
        NoLegend()
}

my.sub.top.0d <- read.table("Top_100_Markers_Subset_Immune_0d.txt", header=T)
my.sub.top.1d <- read.table("Top_100_Markers_Subset_Immune_1d.txt", header=T)
my.sub.top.3d <- read.table("Top_100_Markers_Subset_Immune_3d.txt", header=T)
my.sub.top.7d <- read.table("Top_100_Markers_Subset_Immune_7d.txt", header=T)

my.cluster.0d <- unique(my.sub.top.0d$cluster)
my.cluster.1d <- unique(my.sub.top.1d$cluster)
my.cluster.3d <- unique(my.sub.top.3d$cluster)
my.cluster.7d <- unique(my.sub.top.7d$cluster)

st.sub.top.0d <- read.table("Top_100_Markers_Subset_Structure_0d.txt", header=T)
st.sub.top.1d <- read.table("Top_100_Markers_Subset_Structure_1d.txt", header=T)
st.sub.top.3d <- read.table("Top_100_Markers_Subset_Structure_3d.txt", header=T)
st.sub.top.7d <- read.table("Top_100_Markers_Subset_Structure_7d.txt", header=T)
st.cluster.0d <- unique(st.sub.top.0d$cluster)
st.cluster.1d <- unique(st.sub.top.1d$cluster)
st.cluster.3d <- unique(st.sub.top.3d$cluster)
st.cluster.7d <- unique(st.sub.top.7d$cluster)

temp <- readRDS('CellChat_Choices.RDS')
possible.source.clusters.0d <- temp[[1]]
possible.source.clusters.1d <- temp[[2]]
possible.source.clusters.3d <- temp[[3]]
possible.source.clusters.7d <- temp[[4]]
possible.pathways <- temp[[5]]

summary.df <- data.frame("Total Number of Cells Sequenced" =
                             c('54,413', '5,121', '18,001', '32,907', '17,384',
                               '10,986', '8,349', '9,055', '36,508'),
                         "Number of Cells Analyzed after applying excluding criteria" =
                             c('34,386', '4,539', '17,434', '29,070', '14,581',
                               '10,091', '6,971', '7,370', '32,967'),
                         "Number of Biological Replicates (number of total replicates)" =
                             c('n = 2(2)', 'n = 1(1)', 'n = 1(1)', 'n = 2(5)', 'n = 2(3)',
                               'n = 2(3)', 'n = 1(2)', 'n = 1(2)', 'n = 3(8)'),
                         "Median UMI Counts per Cell" =
                             c('5,174', '6,217', '5,585', '4,677', '8,322',
                               '4,154','7,703','8,784', '7,868'),
                         "Median Genes per Cell" =
                             c('1,520','2,159', '2,099', '1,887', '2,754',
                               '1,687', '2,549', '2,734', '2,549'),
                         "Number of Reads" =
                             c('1.98B', '0.38B', '1.45B', '2.1B', '1.05B',
                               '1.64B', '0.77B', '1.1B',  '2.8B'),
                         "Sequencing Saturation" =
                             c('63.31%', '82.7%', '83%', '78.31%', '61.49%',
                               '89.23%', '78.5%', '83.5%', '78.5%'),
                         row.names=c('PBMC (naive mouse blood)', 'Naive SN (CD45(+) immune)',
                                     'Naive SN (CD45(-) non-immune)', '1d Injured SN',
                                     '3d Injured SN (Kalinski et al., 2020)', '3d Injured SN',
                                     '3d SN CD11b(+), Injury Site', '3d SN CD11b(+), Distal Nerve',
                                     '7d Injured SN'),
                         check.names=F)

ui <- dashboardPage(
    dashboardHeader(title = "iSNAT Mouse"),
    dashboardSidebar(
          sidebarMenu(
              menuItem("Project Summary", tabName = 'Home', icon = icon('cube', lib='font-awesome')),
              menuItem("Expression Analysis", tabName = 'DataClus',
                       icon = icon('cubes', lib='font-awesome')),
              menuItem("Two Genes", tabName = 'TwoGenes', icon = icon('dice-d20', lib='font-awesome')),
              menuItem("Immune Cells", tabName = 'Immune', icon = icon('cubes', lib='font-awesome')),
              menuItem("Structural Cells", tabName='Str', icon = icon('cubes', lib='font-awesome')),
              menuItem("Integrated Immune Cells", tabName = 'int_immune', icon = icon('cubes', lib='font-awesome')),
              menuItem("Trajectory", tabName = 'Trajectory', icon = icon('cubes', lib='font-awesome')),
              menuItem("Cell Chat", tabName='Signal', icon = icon('cubes', lib='font-awesome')),
              menuItem("Ligand-Receptor", tabName='SigPath', icon = icon('cubes', lib='font-awesome')),
              menuItem("Spatial Distribution", tabName='Location', icon = icon('cubes', lib='font-awesome')),
              menuItem("PBMC", tabName='PBMC', icon = icon('cubes', lib='font-awesome'))
          )
      ),
      dashboardBody(
          tabItems(
              tabItem(tabName = 'Home',
                      fluidPage(
                          useWaiter(),
                          waiter_preloader(html=waiting_screen),
                          includeMarkdown('Intro_1.md'),
                          br(),
                          br(),
                          tableHTML_output("summaryTable"),
                          br(),
                          h2("Data Composition"),
                          p('UMAP plots for each time point were generated. Clusters were assigned to specific cell types.'),
                          p('The tables underneath the UMAPs show the number of cells given each label.'),
                          fluidRow(column(3, withSpinner(plotOutput("DimPlot.0d"), type=6),
                                          style='padding:4px;'),
                                   column(3, withSpinner(plotOutput("DimPlot.1d"), type=6),
                                          style='padding:4px;'),
                                   column(3, withSpinner(plotOutput("DimPlot.3d"), type=6),
                                          style='padding:4px;'),
                                   column(3, withSpinner(plotOutput("DimPlot.7d"), type=6),
                                          style='padding:4px;')
                                   ),
                          fluidRow(column(3, box(width=NULL, tableOutput("CellComp.0d"))),
                                   column(3, box(width=NULL, tableOutput("CellComp.1d"))),
                                   column(3, box(width=NULL, tableOutput("CellComp.3d"))),
                                   column(3, box(width=NULL, tableOutput("CellComp.7d")))
                                   ),
                          p(strong('Abbreviations:')),
                          p('GC (granulocytes), Mo (monocytes), Mac (macrophages), MoDC (monocyte-derived dendritic cells), cDC (conventional dendritic cells), pDC (plasmocytoid dendritic cells), DCx (dendritic cells destined for homing),', style='font-size: small; margin-bottom: 3px;'),
                          p('Mast (Mast cells), NK (natural killer cells), TC (T cells), eMES (endoneurial mesenchymal cells), pMES (perineurial mesenchymal cells), dMES (differentiating mesenchymal cells), prol.MES (proliferating mesenchymal cells)', style='font-size: small; margin-bottom: 3px;'),
                          p('Fb (fibroblasts), SC (Schwann cells), EC (endothelial cells), vSMC_PC (vascular smooth muscle cells/pericytes).', style='font-size: small;'),
                          br(),
                          includeMarkdown('Intro_2.md')
                      )
                      ),
              tabItem(tabName = 'DataClus',
                    fluidPage(
                        fluidRow(selectizeInput("Gene1", label = "Gene 1", choices=NULL)),
                        fluidRow(column(3, withSpinner(plotOutput('AllDimPlot.0d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('AllDimPlot.1d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('AllDimPlot.3d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('AllDimPlot.7d'), type=6), style='padding:4px;'))
                       ,
                        br(),
                        fluidRow(column(3, withSpinner(plotOutput('AllFeatPlot.0d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('AllFeatPlot.1d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('AllFeatPlot.3d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('AllFeatPlot.7d'), type=6), style='padding:4px;')
                                 ),
                        fluidRow(column(4, selectInput("Cluster", label="Cluster Box One",
                                                       choices=clusters)),
                                 column(4, offset=2, selectInput("Cluster2", label="Cluster Box Two",
                                                                 choices=clusters))),
                        fluidRow(column(4,
                                        tabBox(width=NULL,
                                            title="Cluster Positive Markers",
                                            tabPanel("Naive Markers", tableOutput("Topmark0d")),
                                            tabPanel("1 Day Markers", tableOutput("Topmark1d")),
                                            tabPanel("3 Day Markers", tableOutput("Topmark3d")),
                                            tabPanel("7 Day Markers", tableOutput("Topmark7d")))),
                                 column(4, offset=2,
                                        tabBox(width=NULL,
                                            title="Cluster Positive Markers Box 2",
                                            tabPanel("Naive Markers", tableOutput("Topmark0d2")),
                                            tabPanel("1 Day Markers", tableOutput("Topmark1d2")),
                                            tabPanel("3 Day Markers", tableOutput("Topmark3d2")),
                                            tabPanel("7 Day Markers", tableOutput("Topmark7d2"))))
                                 )
                    )
                    ),
           tabItem(tabName = 'TwoGenes',
                   fluidPage(
                       h2("Enter two genes and push the Go Button"),
                       fluidRow(column(2, selectizeInput("CompGene1", label = "Gene 1", choices=NULL)),
                                column(2, selectizeInput("CompGene2", label = "Gene 2", choices=NULL)),
                                column(2, actionButton("LoadTwoGenes", label='GO'))),
                       h2("Naive"),
                       fluidRow(column(9, withSpinner(plotOutput('Blend.0d'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Fscatter.0d'), type=6), style='padding:4px;')),
                       h2("1 Day Post Crush"),
                       fluidRow(column(9, withSpinner(plotOutput('Blend.1d'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Fscatter.1d'), type=6), style='padding:4px;')),
                       h2("3 Day Post Crush"),
                       fluidRow(column(9, withSpinner(plotOutput('Blend.3d'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Fscatter.3d'), type=6), style='padding:4px;')),
                       h2("7 Day Post Crush"),
                       fluidRow(column(9, withSpinner(plotOutput('Blend.7d'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Fscatter.7d'), type=6), style='padding:4px;'))
                   )
                   ),
           tabItem(tabName = 'Immune',
                   fluidPage(
                       h2("Zoomed in view of Immune cells From the original UMAP"),
                       fluidRow(column(3, withSpinner(plotOutput("My.DimPlot.0d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("My.DimPlot.1d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("My.DimPlot.3d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("My.DimPlot.7d"), type=6), style='padding:4px;')
                                ),
                       fluidRow(selectizeInput("my.gene", label = "Gene For Feature Plot", choices=NULL)),
                       fluidRow(column(3, withSpinner(plotOutput("My.FeatPlot.0d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("My.FeatPlot.1d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("My.FeatPlot.3d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("My.FeatPlot.7d"), type=6), style='padding:4px;')
                                ),
                       h2("Tables show positively expressed genes in chosen cluster relative to other immune clusters"),
                       fluidRow(column(3, selectInput("my.clus.0d", label="Cluster Naive",
                                                      choices=c("", my.cluster.0d))),
                                column(3, selectInput("my.clus.1d", label="Cluster 1 Day",
                                                      choices=c("", my.cluster.1d))),
                                column(3, selectInput("my.clus.3d", label="Cluster 3 Day",
                                                      choices=c("", my.cluster.3d))),
                                column(3, selectInput("my.clus.7d", label="Cluster 7 Day",
                                                      choices=c("", my.cluster.7d)))
                                ),
                         fluidRow(column(3,
                                         box(width=NULL, title="Naive Clusters",
                                             tableOutput("TopMarks.my.0d"))),
                                  column(3,
                                         box(width=NULL, title="1 Day Clusters",
                                             tableOutput("TopMarks.my.1d"))),
                                  column(3,
                                         box(width=NULL, title="3 Day Clusters",
                                             tableOutput("TopMarks.my.3d"))),
                                  column(3,
                                         box(width=NULL, title="7 Day Clusters",
                                             tableOutput("TopMarks.my.7d"))))
                     )),
           tabItem(tabName = 'Str',
                   fluidPage(
                       h2("Zoomed in view of Structural cells From the original UMAP"),
                       fluidRow(column(3, withSpinner(plotOutput("st.DimPlot.0d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("st.DimPlot.1d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("st.DimPlot.3d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("st.DimPlot.7d"), type=6), style='padding:4px;')
                                ),
                       fluidRow(selectizeInput("st.gene", label = "Gene For Feature Plot", choices=NULL)),
                       fluidRow(column(3, withSpinner(plotOutput("st.FeatPlot.0d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("st.FeatPlot.1d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("st.FeatPlot.3d"), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput("st.FeatPlot.7d"), type=6), style='padding:4px;')
                                ),
                       h2("Tables show positively expressed genes in chosen cluster relative to other structural clusters"),
                       fluidRow(column(3, selectInput("st.clus.0d", label="Cluster Naive",
                                                      choices=c("", st.cluster.0d))),
                                column(3, selectInput("st.clus.1d", label="Cluster 1 Day",
                                                      choices=c("", st.cluster.1d))),
                                column(3, selectInput("st.clus.3d", label="Cluster 3 Day",
                                                      choices=c("", st.cluster.3d))),
                                column(3, selectInput("st.clus.7d", label="Cluster 7 Day",
                                                      choices=c("", st.cluster.7d)))
                                ),
                       fluidRow(column(3,
                                       box(width=NULL, title="Naive Clusters",
                                           tableOutput("TopMarks.st.0d"))),
                                column(3,
                                       box(width=NULL, title="1 Day Clusters",
                                           tableOutput("TopMarks.st.1d"))),
                                column(3,
                                       box(width=NULL, title="3 Day Clusters",
                                           tableOutput("TopMarks.st.3d"))),
                                column(3,
                                       box(width=NULL, title="7 Day Clusters",
                                           tableOutput("TopMarks.st.7d"))))
                   )),
           tabItem(tabName = 'int_immune',
                   fluidPage(
                       fluidRow(actionButton("LoadIntegrated",'Load Integrated Immune Data')),
                       br(),
                       h2("Integrated UMAPs with cell colors and labels from Individual analysis"),
                       fluidRow(column(4, style='padding:4px;',
                                       column(6, withSpinner(plotOutput('Int.DimPlot.0d'), type=6), style='padding:4px;'),
                                       column(6, withSpinner(plotOutput('Int.DimPlot.b'), type=6), style='padding:4px;')),
                                column(8, style='padding:4px;',
                                       column(4, withSpinner(plotOutput('Int.DimPlot.1d'), type=6), style='padding:4px;'),
                                       column(4, withSpinner(plotOutput('Int.DimPlot.3d'), type=6), style='padding:4px;'),
                                       column(4, withSpinner(plotOutput('Int.DimPlot.7d'), type=6), style='padding:4px;'))
                                ),
                       fluidRow(selectizeInput("IntGene", label = "Gene for FeaturePlot", choices=NULL)),
                       fluidRow(column(4, style='padding:4px;',
                                       column(6, withSpinner(plotOutput('Int.0d.feature'), type=6), style='padding:4px;'),
                                       column(6, withSpinner(plotOutput('Int.b.feature'), type=6), style='padding:4px;')
                                       ),
                                column(8, style='padding:4px;',
                                       column(4, withSpinner(plotOutput('Int.1d.feature'), type=6), style='padding:4px;'),
                                       column(4, withSpinner(plotOutput('Int.3d.feature'), type=6), style='padding:4px;'),
                                       column(4, withSpinner(plotOutput('Int.7d.feature'), type=6), style='padding:4px;')
                                       )
                                ),
                       br(),
                       h2("Integrated UMAPs showing Seurat Cluster labels on integrated anlaysis"),
                       fluidRow(column(4, style='padding:4px;',
                                       column(6, withSpinner(plotOutput('Int.0d.seurat'), type=6), style='padding:4px;'),
                                       column(6, withSpinner(plotOutput('Int.b.seurat'), type=6), style='padding:4px;')
                                       ),
                                column(8, style='padding:4px;',
                                       column(4, withSpinner(plotOutput('Int.1d.seurat'), type=6), style='padding:4px;'),
                                       column(4, withSpinner(plotOutput('Int.3d.seurat'), type=6), style='padding:4px;'),
                                       column(4, withSpinner(plotOutput('Int.7d.seurat'), type=6), style='padding:4px;')
                                       ),
                                fluidRow(column(2, selectInput('Int.Cluster', label='Seurat Cluster',
                                                               choices=as.character(0:22),
                                                               selected='0'))),
                                ## fluidRow(column(4, box(width=NULL, title='Markers',
                                ##                        tableOutput("Int_Seurat_Markers"))))
                                fluidRow(column(6, dataTableOutput("Int_Seurat_Markers"))
                                         )
                                )
                   )),
           tabItem(tabName = 'Trajectory',
                   fluidPage(
                       h4("First: if you have not already, load the integrated data in the 'Integrated Immune Cells' tab"),
                       h4("Second: Load Trajectory data with this button"),
                       fluidRow(actionButton("LoadTrajectory",'Load Trajectory Data')),
                       br(),
                       fluidRow(column(2, selectInput('Traj', label='Trajectory', selected='Mo_to_Mac',
                                                      choices=c('Mo_to_Mac','Mo_to_DC'))),
                                column(2, selectInput('LabelBy', label='Label Clusters By',
                                                      selected='Trajectory_CellTypes',
                                                      choices=c('Trajectory_CellTypes', 'seurat_clusters')))
                                ),
                       fluidRow(column(3, withSpinner(plotOutput('Traj_PCA'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Traj_UMAP'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Traj_PCA_pseudo'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Traj_UMAP_pseudo'), type=6), style='padding:4px;')
                                ),
                       br(),
                       fluidRow(column(4, dataTableOutput("Var_Imp")),
                                column(2, selectizeInput('TrajGene', label = 'Gene:', selected='Arg1',
                                                         choices=NULL)),
                                column(3, withSpinner(plotOutput('Traj_PCA_Feature'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Traj_UMAP_Feature'), type=6), style='padding:4px;')
                                )
                   )),
           
            tabItem(tabName = 'Signal',
                    fluidPage(
                        fluidRow(actionButton("LoadCellChat",'Load Cell Chat Data')),
                        br(),
                        fluidRow(column(3, selectInput("source.0d", label="Day 0 Source Cell Type",
                                                       choices=possible.source.clusters.0d)),
                                 column(3, selectInput("source.1d", label="Day 1 Source Cell Type",
                                                       choices=possible.source.clusters.1d)),
                                 column(3, selectInput("source.3d", label="Day 3 Source Cell Type",
                                                       choices=possible.source.clusters.3d)),
                                 column(3, selectInput("source.7d", label="Day 7 Source Cell Type",
                                                       choices=possible.source.clusters.7d))
                                 ),
                        fluidRow(column(3, box(width=NULL, title="L-R Day 0", tableOutput("L_R_Day_0"))),
                                 column(3, box(width=NULL, title="L-R Day 1", tableOutput("L_R_Day_1"))),
                                 column(3, box(width=NULL, title="L-R Day 3", tableOutput("L_R_Day_3"))),
                                 column(3, box(width=NULL, title="L-R Day 7", tableOutput("L_R_Day_7")))
                                 ),
                        fluidRow(column(3, withSpinner(plotOutput('NetVisCircle.0d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('NetVisCircle.1d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('NetVisCircle.3d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('NetVisCircle.7d'), type=6), style='padding:4px;')
                        )
                    )),
            tabItem(tabName = 'SigPath',
                    fluidPage(
                        h4('Data for this page must be loaded in the Cell Chat tab'),
                        br(),
                        fluidRow(selectInput("signal.path", label=" Signal Pathway",
                                             choices=possible.pathways, selected='CXCL')),
                        fluidRow(column(3, h5("0 Day Naive")),
                                 column(3, h5("1 Day post-injury")),
                                 column(3, h5("3 Day post-injury")),
                                 column(3, h5("7 Day post-injury"))),
                        fluidRow(column(3, withSpinner(plotOutput('heatmap.path.0d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('heatmap.path.1d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('heatmap.path.3d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('heatmap.path.7d'), type=6), style='padding:4px;')
                                 ),
                        br(),
                        fluidRow(column(3, withSpinner(plotOutput('chord.path.0d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('chord.path.1d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('chord.path.3d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('chord.path.7d'), type=6), style='padding:4px;')
                                 ),
                        br(),
                        fluidRow(column(3, withSpinner(plotOutput('contrib.path.0d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('contrib.path.1d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('contrib.path.3d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('contrib.path.7d'), type=6), style='padding:4px;')
                                 ),
                        br(),
                        fluidRow(column(3, selectizeInput("LR_Pair_0d", label="L-R Pair 0d", choices=NULL)),
                                 column(3, selectizeInput("LR_Pair_1d", label="L-R Pair 1d", choices=NULL)),
                                 column(3, selectizeInput("LR_Pair_3d", label="L-R Pair 3d", choices=NULL)),
                                 column(3, selectizeInput("LR_Pair_7d", label="L-R Pair 7d", choices=NULL))
                                 ),
                        fluidRow(column(3, withSpinner(plotOutput('LR_Pair.0d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('LR_Pair.1d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('LR_Pair.3d'), type=6), style='padding:4px;'),
                                 column(3, withSpinner(plotOutput('LR_Pair.7d'), type=6), style='padding:4px;')
                                 )
                    )),
           tabItem(tabName = 'Location',
                   fluidPage(
                       fluidRow(column(3, withSpinner(plotOutput('Location.UMAP'), type=6), style='padding:4px;'),
                                column(6, withSpinner(plotOutput('Location.UMAP.Split'), type=6), style='padding:4px;'),
                                column(3, box(width=NULL, tableOutput("Loc.Cluster.Size")))
                                ),
                       fluidRow(column(3, selectizeInput("loc.gene", label = "Gene For Feature Plot",
                                                         choices=NULL))
                                ),
                       fluidRow(column(3, withSpinner(plotOutput('Location.Feature'), type=6), style='padding:4px;'),
                                column(6, withSpinner(plotOutput('Location.Feature.Split'), type=6), style='padding:4px;')
                                ),
                       fluidRow(column(3, selectInput("loc.clus", label="Cluster",
                                                      choices=c("Mo", "Mac1", "Mac3", "Mac4",
                                                                "epiMac", "prol.Mac", "GC", "MoDC",
                                                                "cDC", "DCx", "TC", "NK", "UI-M1",
                                                                "UI-M2")))
                                ),
                       fluidRow(column(4, box(width=NULL, title="Cluster Markers",
                                              tableOutput("loc.markers"))),
                                column(4, box(width=NULL, title="Distal", tableOutput("DE.distal"))),
                                column(4, box(width=NULL, title="Injury", tableOutput("DE.injury")))
                                )
                   )),
           tabItem(tabName = 'PBMC',
                   fluidPage(
                       fluidRow(actionButton("LoadPBMCData",'Load PBMC Data')),
                       fluidRow(selectizeInput("Gene.pbmc", label = "Gene 1", choices=NULL)),
                       fluidRow(column(3, withSpinner(plotOutput('AllDimPlot.blood'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('AllFeaturePlot.blood'), type=6), style='padding:4px;')
                                ),
                       br(),
                       fluidRow(selectizeInput("Gene2.pbmc", label = "Gene 2", choices=NULL)),
                       fluidRow(column(9, withSpinner(plotOutput('Blend.blood'), type=6), style='padding:4px;'),
                                column(3, withSpinner(plotOutput('Fscatter.blood'), type=6), style='padding:4px;')
                                ),
                       fluidRow(column(4, selectInput("Cluster.blood", label="Cluster for Markers",
                                                      choices=clusters))
                                ),
                       fluidRow(column(4, box(width=NULL, title="Cluster Markers",
                                              tableOutput("blood.markers")))
                                )
                   ))
          )
      )
)

server <- function(input, output, session){
### selectize
    genes.to.select <- readRDS('All_Gene_Names.RDS')
    loc.genes.to.select <- readRDS('Location_Gene_Names.RDS')
    updateSelectizeInput(session, 'Gene1', choices = genes.to.select, selected='Arg1', server = TRUE)
    updateSelectizeInput(session, 'CompGene1', choices = genes.to.select, selected=character(0), server = TRUE)
    updateSelectizeInput(session, 'CompGene2', choices = genes.to.select, selected=character(0), server = TRUE)
    updateSelectizeInput(session, 'my.gene', choices = genes.to.select, selected=character(0), server = TRUE)
    updateSelectizeInput(session, 'sc.gene', choices = genes.to.select, selected=character(0), server = TRUE)
    updateSelectizeInput(session, 'st.gene', choices = genes.to.select, selected=character(0), server = TRUE)
    updateSelectizeInput(session, 'loc.gene', choices = loc.genes.to.select, selected=character(0), server = TRUE)
    updateSelectizeInput(session, 'IntGene', choices = genes.to.select, selected='Arg1', server = TRUE)
    updateSelectizeInput(session, 'TrajGene', choices = genes.to.select, selected='Arg1', server = TRUE)

    ## Big files
    s.obj.0d <- readRDS("post_injury_day_0_Diet.RDS")
    s.obj.1d <- readRDS("post_injury_day_1_Diet.RDS")
    s.obj.3d <- readRDS("post_injury_day_3_Diet.RDS")
    s.obj.7d <- readRDS("post_injury_day_7_Diet.RDS")
    s.obj.loc <- readRDS('Location_MyeloidOnly_Diet.RDS')

    Idents(s.obj.0d) <- 'Cell.Type'
    Idents(s.obj.1d) <- 'Cell.Type'
    Idents(s.obj.3d) <- 'Cell.Type'
    Idents(s.obj.7d) <- 'Cell.Type'

    orig.UMAP.0d <- DimPlot(s.obj.0d, cols=gg_cols, group.by='Cell.Type', pt.size=.2) +
        NoLegend() + ggtitle("Naive Sciatic Nerve")
    orig.UMAP.0d <- LabelClusters(orig.UMAP.0d, id='Cell.Type', size=4, repel=1, box.padding=1)

    orig.UMAP.1d <- DimPlot(s.obj.1d, cols=gg_cols, group.by='Cell.Type', pt.size=.2) +
        NoLegend() + ggtitle("1 Day Post-Injury")
    orig.UMAP.1d <- LabelClusters(orig.UMAP.1d, id='Cell.Type', size=4, repel=1, box.padding=1)

    orig.UMAP.3d <- DimPlot(s.obj.3d, cols=gg_cols, group.by='Cell.Type', pt.size=.2) +
        NoLegend() + ggtitle("3 Day Post-Injury")
    orig.UMAP.3d <- LabelClusters(orig.UMAP.3d, id='Cell.Type', size=4, repel=1, box.padding=1)

    orig.UMAP.7d <- DimPlot(s.obj.7d, cols=gg_cols, group.by='Cell.Type', pt.size=.2) +
        NoLegend() + ggtitle("7 Day Post-Injury")
    orig.UMAP.7d <- LabelClusters(orig.UMAP.7d, id='Cell.Type', size=4, repel=1, box.padding=1)

    CC.Data <- eventReactive(input$LoadCellChat, {
        w.load.CC$show()
        CellChat.0d <- readRDS('CellChat_day_0.RDS')
        df.net.0d <- subsetCommunication(CellChat.0d)
        pathways.0d <- unique(df.net.0d$pathway_name)
        CellChat.1d <- readRDS('CellChat_day_1.RDS')
        df.net.1d <- subsetCommunication(CellChat.1d)
        pathways.1d <- unique(df.net.1d$pathway_name)
        CellChat.3d <- readRDS('CellChat_day_3.RDS')
        df.net.3d <- subsetCommunication(CellChat.3d)
        pathways.3d <- unique(df.net.3d$pathway_name)
        CellChat.7d <- readRDS('CellChat_day_7.RDS')
        df.net.7d <- subsetCommunication(CellChat.7d)
        pathways.7d <- unique(df.net.7d$pathway_name)
        w.load.CC$hide()
        list('CC.0d'=CellChat.0d, 'df.0d'=df.net.0d, 'path.0d'=pathways.0d,
             'CC.1d'=CellChat.1d, 'df.1d'=df.net.1d, 'path.1d'=pathways.1d,
             'CC.3d'=CellChat.3d, 'df.3d'=df.net.3d, 'path.3d'=pathways.3d,
             'CC.7d'=CellChat.7d, 'df.7d'=df.net.7d, 'path.7d'=pathways.7d)
    })
    Int.Data <- eventReactive(input$LoadIntegrated, {
        w.load.Int$show()
        Int <- readRDS('Myeloid_Integration_Diet.RDS')
        w.load.Int$hide()
        Int
    })

    PBMC.Data <- eventReactive(input$LoadPBMCData, {
        w.load.blood$show()
        temp.data <- readRDS('Peripheral_Blood_Diet.RDS')
        w.load.blood$hide()
        temp.data
    })

    ### Waiters
    w.load.CC <- Waiter$new(html=h1('Loading Cell Chat Data'))
    w.load.Int <- Waiter$new(html=h1('Loading Integrated Immune Data'))
    w.load.TG <- Waiter$new(html=h1(HTML(paste('Allow 10 seconds for plot rendering',
                                               'after this screen disappears',
                                               sep='<br/>'))))
    w.load.Traj <- Waiter$new(html=h1('Loading Trajectory Data'))
    w.load.blood <- Waiter$new(html=h1('Loading PBMC Data'))
##    w1 <- Waiter$new(html=h1('Hold It!'))
##    w2 <- Waiter$new(html=h1('It takes me a while to do this, Grab a cup of coffee!'))
    output$summaryTable <- render_tableHTML(
        tableHTML(summary.df, collapse = 'separate_shiny', caption="scRNA-seq DataSets",
                  widths=c(190, rep(150,3), rep(100,4))) %>%
        add_css_column(css=list('text-align', 'center'), columns=names(summary.df)) %>%
        add_css_header(css = list(c('height', 'background-color', 'text-align'), c('30px', '#e6e6e6', 'center')),
                       headers = 1:8) %>%
        add_css_column(css = list(c('font-weight', 'background-color'), c('bold', '#e6e6e6')),
                       columns = 'rownames')
    )

    Gene1 <- reactive(str_to_title(input$Gene1))
    CompGene1 <- reactive(str_to_title(input$CompGene1))
    CompGene2 <- reactive(str_to_title(input$CompGene2))
    my.gene <- reactive(str_to_title(input$my.gene))
    sc.gene <- reactive(str_to_title(input$sc.gene))
    st.gene <- reactive(str_to_title(input$st.gene))
    loc.gene <- reactive(str_to_title(input$loc.gene))

###################################################################################################
###### Summary Data Composition #####
###### UMAPS in first Row  #####
    output$DimPlot.0d <- renderPlot({ orig.UMAP.0d })
    output$DimPlot.1d <- renderPlot({ orig.UMAP.1d })
    output$DimPlot.3d <- renderPlot({ orig.UMAP.3d })
    output$DimPlot.7d <- renderPlot({ orig.UMAP.7d })

###### Cell Count Table  #####

    output$CellComp.0d <- renderTable({
        cell.comp.0d <- table(s.obj.0d$Cell.Type)
        cell.comp.0d <- data.frame(cbind(names(cell.comp.0d)[1:9], cell.comp.0d[1:9],
                                         names(cell.comp.0d)[10:18], cell.comp.0d[10:18]),
                                   row.names=NULL)
        colnames(cell.comp.0d) <- c('CellType','Count','CellType','Count')
        cell.comp.0d
    })
    output$CellComp.1d <- renderTable({
        cell.comp.1d <- table(s.obj.1d$Cell.Type)
        cell.comp.1d <- data.frame(cbind(names(cell.comp.1d)[1:9], cell.comp.1d[1:9],
                                         names(cell.comp.1d)[10:18], cell.comp.1d[10:18]),
                                    row.names=NULL)
        colnames(cell.comp.1d) <- c('CellType','Count','CellType','Count')
        cell.comp.1d
    })
    output$CellComp.3d <- renderTable({
        cell.comp.3d <- table(s.obj.3d$Cell.Type)
        cell.comp.3d <- data.frame(cbind(names(cell.comp.3d)[1:9], cell.comp.3d[1:9],
                                         names(cell.comp.3d)[10:18], cell.comp.3d[10:18]),
                                    row.names=NULL)
        colnames(cell.comp.3d) <- c('CellType','Count','CellType','Count')
        cell.comp.3d
    })
    output$CellComp.7d <- renderTable({
        cell.comp.7d <- table(s.obj.7d$Cell.Type)
        cell.comp.7d <- data.frame(cbind(names(cell.comp.7d)[1:9], cell.comp.7d[1:9],
                                         names(cell.comp.7d)[10:18], cell.comp.7d[10:18]),
                                    row.names=NULL)
        colnames(cell.comp.7d) <- c('CellType','Count','CellType','Count')
        cell.comp.7d
    })

###################################################################################################
###### Expression Analysis #####
###### UMAPS #####
    output$AllDimPlot.0d <- renderPlot({ orig.UMAP.0d })
    output$AllDimPlot.1d <- renderPlot({ orig.UMAP.1d })
    output$AllDimPlot.3d <- renderPlot({ orig.UMAP.3d })
    output$AllDimPlot.7d <- renderPlot({ orig.UMAP.7d })

###################################################################################################
###### Expression Analysis - FeaturePlots 2nd Row #####
    output$AllFeatPlot.0d <- renderPlot({
        req(Gene1() %in% rownames(s.obj.0d))
        print(FeaturePlot(s.obj.0d, slot='data', features=Gene1(), order=T, min.cutoff=1))
    })
    output$AllFeatPlot.1d <- renderPlot({
        req(Gene1() %in% rownames(s.obj.1d))
        print(FeaturePlot(s.obj.1d, slot='data', features=Gene1(), order=T, min.cutoff=1))
    })
    output$AllFeatPlot.3d <- renderPlot({
        req(Gene1() %in% rownames(s.obj.3d))
        print(FeaturePlot(s.obj.3d, slot='data', features=Gene1(), order=T, min.cutoff=1))
    })
    output$AllFeatPlot.7d <- renderPlot({
        req(Gene1() %in% rownames(s.obj.7d))
        print(FeaturePlot(s.obj.7d, slot='data', features=Gene1(), order=T, min.cutoff=1))
    })
###################################################################################################
###### Expression Analysis - Top Markers 3rd Row #####
    output$Topmark0d <- renderTable({
        req(input$Cluster %in% clusters)
        marks.0d <- top.0d %>% filter(cluster==input$Cluster) %>%
            select(-p_val) %>% as.data.frame()
        if(nrow(marks.0d)==0){
            paste0('There are either no cells or no markers for cluster ',
                   input$Cluster,' at this time point')
        } else {
            marks.0d %>% relocate(gene)
        }
    })
    output$Topmark1d <- renderTable({
        req(input$Cluster %in% clusters)
        marks.1d <- top.1d %>% filter(cluster==input$Cluster) %>%
            select(-p_val) %>% as.data.frame()
        if(nrow(marks.1d)==0){
            paste0('There are either no cells or no markers for cluster ',
                   input$Cluster,' at this time point')
        } else {
            marks.1d %>% relocate(gene)
        }
    })
    output$Topmark3d <- renderTable({
        req(input$Cluster %in% clusters)
        marks.3d <- top.3d %>% filter(cluster==input$Cluster) %>%
            select(-p_val) %>% as.data.frame()
        if(nrow(marks.3d)==0){
            paste0('There are either no cells or no markers for cluster ',
                   input$Cluster,' at this time point')
        } else {
            marks.3d %>% relocate(gene)
        }
    })
    output$Topmark7d <- renderTable({
        req(input$Cluster %in% clusters)
        marks.7d <- top.7d %>% filter(cluster==input$Cluster) %>%
            select(-p_val) %>% as.data.frame()
        if(nrow(marks.7d)==0){
            paste0('There are either no cells or no markers for cluster ',
                   input$Cluster,' at this time point')
        } else {
            marks.7d %>% relocate(gene)
        }
    })
###################################################################################################
###### Expression Analysis - Top Markers 3rd Row, 2nd Column #####
    output$Topmark0d2 <- renderTable({
        req(input$Cluster2 %in% clusters)
        marks.0d2 <- top.0d %>% filter(cluster==input$Cluster2) %>% as.data.frame()
        if(nrow(marks.0d2)==0){
            paste0('There are either no cells or no markers for cluster ',
                   input$Cluster2,' at this time point')
        } else {
            marks.0d2 %>% relocate(gene)
        }
    })
    output$Topmark1d2 <- renderTable({
        req(input$Cluster2 %in% clusters)
        marks.1d2 <- top.1d %>% filter(cluster==input$Cluster2) %>% as.data.frame()
        if(nrow(marks.1d2)==0){
            paste0('There are either no cells or no markers for cluster ',
                   input$Cluster2,' at this time point')
        } else {
            marks.1d2 %>% relocate(gene)
        }
    })
    output$Topmark3d2 <- renderTable({
        req(input$Cluster2 %in% clusters)
        marks.3d2 <- top.3d %>% filter(cluster==input$Cluster2) %>% as.data.frame()
        if(nrow(marks.3d2)==0){
            paste0('There are either no cells or no markers for cluster ',
                   input$Cluster2,' at this time point')
        } else {
            marks.3d2 %>% relocate(gene)
        }
    })
    output$Topmark7d2 <- renderTable({
        req(input$Cluster2 %in% clusters)
        marks.7d2 <- top.7d %>% filter(cluster==input$Cluster2) %>% as.data.frame()
        if(nrow(marks.7d2)==0){
            paste0('There are either no cells or no markers for cluster ',
                   input$Cluster2,' at this time point')
        } else {
            marks.7d2 %>% relocate(gene)
        }
    })
###################################################################################################
###### Two Genes #####
    TwoGene.Data <- eventReactive(input$LoadTwoGenes, {
        w.load.TG$show()
        blend.0d <- tryCatch(FeaturePlot(s.obj.0d, min.cutoff=1, order=T, slot='data',
                                         features=c(input$CompGene1, input$CompGene2), blend=T, pt.size=0.2),
                             error = function(e){e$message})
        scatter.0d <- FeatureScatter(s.obj.0d, feature1=input$CompGene1, feature2=input$CompGene2,
                                     group.by='Cell.Type', plot.cor=T, pt.size=0.2)
        blend.1d <- tryCatch(FeaturePlot(s.obj.1d, min.cutoff=1, order=T, slot='data',
                                         features=c(input$CompGene1, input$CompGene2), blend=T, pt.size=0.2),
                             error = function(e){e$message})
        scatter.1d <- FeatureScatter(s.obj.1d, feature1=input$CompGene1, feature2=input$CompGene2,
                                     group.by='Cell.Type', plot.cor=T, pt.size=0.2)
        blend.3d <- tryCatch(FeaturePlot(s.obj.3d, min.cutoff=1, order=T, slot='data',
                                         features=c(input$CompGene1, input$CompGene2), blend=T, pt.size=0.2),
                             error = function(e){e$message})
        scatter.3d <- FeatureScatter(s.obj.3d, feature1=input$CompGene1, feature2=input$CompGene2,
                                     group.by='Cell.Type', plot.cor=T, pt.size=0.2)
        blend.7d <- tryCatch(FeaturePlot(s.obj.7d, min.cutoff=1, order=T, slot='data',
                                         features=c(input$CompGene1, input$CompGene2), blend=T, pt.size=0.2),
                             error = function(e){e$message})
        scatter.7d <- FeatureScatter(s.obj.7d, feature1=input$CompGene1, feature2=input$CompGene2,
                                     group.by='Cell.Type', plot.cor=T, pt.size=0.2)
        test <- any(is.character(blend.0d), is.character(blend.1d), is.character(blend.3d), is.character(blend.7d))
        output <- list('Blend.0d'=blend.0d, 'Scatter.0d'=scatter.0d,
                       'Blend.1d'=blend.1d, 'Scatter.1d'=scatter.1d,
                       'Blend.3d'=blend.3d, 'Scatter.3d'=scatter.3d,
                       'Blend.7d'=blend.7d, 'Scatter.7d'=scatter.7d,
                       'Test'=test)
        w.load.TG$hide()
        output
    })

    ## This test and subsequent 'if' statement puts plot generation at 10 seconds when there are no errors.
    ## 25 seconds with errors.
    twogene.test <- reactive({ TwoGene.Data()[['Test']] })
    output$Blend.0d <- renderPlot({
        if(twogene.test() == FALSE){
            print(TwoGene.Data()[['Blend.0d']])
        } else {
            validate(
                need(!is.character(TwoGene.Data()[['Blend.0d']]), paste(TwoGene.Data()[['Blend.0d']]))
            )
            print(TwoGene.Data()[['Blend.0d']])
        }
    })
    output$Fscatter.0d <- renderPlot({
        if(twogene.test() == FALSE){
            print(TwoGene.Data()[['Scatter.0d']])
        } else {
            validate(
                need(!is.character(TwoGene.Data()[['Blend.0d']]), paste(TwoGene.Data()[['Blend.0d']]))
            )
            print(TwoGene.Data()[["Scatter.0d"]])
        }
    })
    output$Blend.1d <- renderPlot({
        if(twogene.test() == FALSE){
            print(TwoGene.Data()[['Blend.1d']])
        } else {
            validate(
                need(!is.character(TwoGene.Data()[['Blend.1d']]), paste(TwoGene.Data()[['Blend.1d']]))
            )
            print(TwoGene.Data()[["Blend.1d"]])
        }
    })
    output$Fscatter.1d <- renderPlot({
        if(twogene.test() == FALSE){
            print(TwoGene.Data()[['Scatter.1d']])
        } else {
            validate(
                need(!is.character(TwoGene.Data()[['Blend.1d']]), paste(TwoGene.Data()[['Blend.1d']]))
            )
            print(TwoGene.Data()[["Scatter.1d"]])
        }
    })
    output$Blend.3d <- renderPlot({
        if(twogene.test() == FALSE){
            print(TwoGene.Data()[['Blend.3d']])
        } else {
            validate(
                need(!is.character(TwoGene.Data()[['Blend.3d']]), paste(TwoGene.Data()[['Blend.3d']]))
            )
            print(TwoGene.Data()[["Blend.3d"]])
        }
    })
    output$Fscatter.3d <- renderPlot({
        if(twogene.test() == FALSE){
            print(TwoGene.Data()[['Scatter.3d']])
        } else {
            validate(
                need(!is.character(TwoGene.Data()[['Blend.3d']]), paste(TwoGene.Data()[['Blend.3d']]))
            )
            print(TwoGene.Data()[["Scatter.3d"]])
        }
    })
    output$Blend.7d <- renderPlot({
        if(twogene.test() == FALSE){
            print(TwoGene.Data()[['Blend.7d']])
        } else {
            validate(
                need(!is.character(TwoGene.Data()[['Blend.7d']]), paste(TwoGene.Data()[['Blend.7d']]))
            )
            print(TwoGene.Data()[["Blend.7d"]])
        }
    })
    output$Fscatter.7d <- renderPlot({
        if(twogene.test() == FALSE){
            print(TwoGene.Data()[['Scatter.7d']])
        } else {
            validate(
                need(!is.character(TwoGene.Data()[['Blend.7d']]), paste(TwoGene.Data()[['Blend.7d']]))
            )
            print(TwoGene.Data()[["Scatter.7d"]])
        }
    })
    
###################################################################################################
###### Immune Cell Clusters #####
    my.sub.0d <- subset(s.obj.0d, subset = Manuscript %in% c('Mac0','epi-Mac','endo-Mac','DC','MAST',
                                                             'NK','TC','prol.Mac','UI1'))
    my.sub.1d <- subset(s.obj.1d, subset = Manuscript %in% c("Mo", "Mac-I", "Mac-II", "Mac-III",
                                                             "Mac-IV", "Mac-V", "MoDC", "cDC", "DCx",
                                                             "pDC", "MAST", "GC", "GC-L", "T_NK",
                                                             "UI-M1"))
    my.sub.3d <- subset(s.obj.3d, subset = Manuscript %in% c("Mo", "Mac1", "Mac2", "Mac3", "Mac4",
                                                             "prol.Mac", "MoDC", "cDC", "DCx", "pDC",
                                                             "MAST", "GC", "T_NK"))
    my.sub.7d <- subset(s.obj.7d, subset = Manuscript %in% c("Mo", "Mac-a", "Mac-b", "Mac-c", "Mac-d",
                                                             "Mac-e", "prol.Mac", "MoDC", "cDC1",
                                                             "cDC2", "DCx", "pDC", "prol.DC", "MAST",
                                                             "GC", "NK", "TC", "TC_NK"))
    output$My.DimPlot.0d <- renderPlot({
        dim.0d.mps <- DimPlot(my.sub.0d, group.by='Manuscript') + xlim(-14, -5) + ylim(-10,13) +
            NoLegend() + ggtitle("Naive Sciatic Nerve")
        suppressWarnings(print(LabelClusters(dim.0d.mps, id='Manuscript', size=4, repel=1, box.padding=1)))
    })
    output$My.DimPlot.1d <- renderPlot({
        dim.1d.mps <- DimPlot(my.sub.1d, group.by='Manuscript') + xlim(-13, 5) + ylim(-10, 16) +
            NoLegend() + ggtitle("1 Day Post-Injury")
        suppressWarnings(print(LabelClusters(dim.1d.mps, id='Manuscript', size=4, repel=1, box.padding=1)))
    })
    output$My.DimPlot.3d <- renderPlot({
        dim.3d.mps <- DimPlot(my.sub.3d, group.by='Manuscript') + xlim(-1, 11) + ylim(-15, 10) +
            NoLegend() + ggtitle("3 Day Post-Injury")
        suppressWarnings(print(LabelClusters(dim.3d.mps, id='Manuscript', size=4, repel=1, box.padding=1)))
    })
    output$My.DimPlot.7d <- renderPlot({
        dim.7d.mps <- DimPlot(my.sub.7d, group.by='Manuscript') + xlim(-15, 1) + ylim(-9, 12) +
            NoLegend() + ggtitle("7 Day Post-Injury")
        suppressWarnings(print(LabelClusters(dim.7d.mps, id='Manuscript', size=4, repel=1, box.padding=1)))
    })

    output$My.FeatPlot.0d <- renderPlot({
        req(my.gene() %in% rownames(s.obj.0d))
        suppressWarnings(print(FeaturePlot(my.sub.0d, slot='data', features=my.gene(),
                                           order=T, min.cutoff=1) +
                               xlim(-14, -5) + ylim(-10,13)))
    })
    output$My.FeatPlot.1d <- renderPlot({
        req(my.gene() %in% rownames(s.obj.1d))
        suppressWarnings(print(FeaturePlot(my.sub.1d, slot='data', features=my.gene(),
                                           order=T, min.cutoff=1) +
                               xlim(-13, 5) + ylim(-10, 16)))
    })
    output$My.FeatPlot.3d <- renderPlot({
        req(my.gene() %in% rownames(s.obj.3d))
        suppressWarnings(print(FeaturePlot(my.sub.3d, slot='data', features=my.gene(),
                                           order=T, min.cutoff=1) +
                               xlim(-1, 11) + ylim(-15, 10)))
    })
    output$My.FeatPlot.7d <- renderPlot({
        req(my.gene() %in% rownames(s.obj.7d))
        suppressWarnings(print(FeaturePlot(my.sub.7d, slot='data', features=my.gene(),
                                           order=T, min.cutoff=1) +
                               xlim(-15, 1) + ylim(-9, 12)))
    })
    output$TopMarks.my.0d <- renderTable({
        my.sub.top.0d %>% ungroup() %>% filter(cluster==input$my.clus.0d) %>%
            select(-p_val, -p_val_adj, -cluster) %>% relocate(gene) %>% as.data.frame()
    })
    output$TopMarks.my.1d <- renderTable({
        my.sub.top.1d %>% ungroup() %>% filter(cluster==input$my.clus.1d) %>%
            select(-p_val, -p_val_adj, -cluster) %>% relocate(gene) %>% as.data.frame()
    })
    output$TopMarks.my.3d <- renderTable({
        my.sub.top.3d %>% ungroup() %>% filter(cluster==input$my.clus.3d) %>%
            select(-p_val, -p_val_adj, -cluster) %>% relocate(gene) %>% as.data.frame()
    })
    output$TopMarks.my.7d <- renderTable({
        my.sub.top.7d %>% ungroup() %>% filter(cluster==input$my.clus.7d) %>%
            select(-p_val, -p_val_adj, -cluster) %>% relocate(gene) %>% as.data.frame()
    })
    
###################################################################################################
###### Structure Cell Clusters #####
    st.sub.0d <- subset(s.obj.0d, subset = Cell.Type %in% c("Fb", "dMES", "prol.MES", "eMES", "pMES"))
    st.sub.1d <- subset(s.obj.1d, subset = Cell.Type %in% c("Fb", "dMES", "prol.MES", "eMES", "pMES"))
    st.sub.3d <- subset(s.obj.3d, subset = Cell.Type %in% c("Fb", "dMES", "prol.MES", "eMES", "pMES"))
    st.sub.7d <- subset(s.obj.7d, subset = Cell.Type %in% c("Fb", "dMES", "prol.MES", "eMES", "pMES"))
    output$st.DimPlot.0d <- renderPlot({
        dim.0d.st <- DimPlot(st.sub.0d, group.by='Manuscript') + xlim(-13, 12) + ylim(-9, 11) +
            NoLegend() + ggtitle("Naive Sciatic Nerve")
        suppressWarnings(print(LabelClusters(dim.0d.st, id='Manuscript', size=4, repel=1, box.padding=1)))
    })
    output$st.DimPlot.1d <- renderPlot({
        dim.1d.st <- DimPlot(st.sub.1d, group.by='Manuscript') + xlim(6,14) + ylim(-7,9) +
            NoLegend() + ggtitle("1 Day Post-Injury")
        suppressWarnings(print(LabelClusters(dim.1d.st, id='Manuscript', size=4, repel=1, box.padding=1)))
    })
    output$st.DimPlot.3d <- renderPlot({
        dim.3d.st <- DimPlot(st.sub.3d, group.by='Manuscript') + xlim(-15, -3) + ylim(-13, 7) +
            NoLegend() + ggtitle("3 Day Post-Injury")
        suppressWarnings(print(LabelClusters(dim.3d.st, id='Manuscript', size=4, repel=1, box.padding=1)))
    })
    output$st.DimPlot.7d <- renderPlot({
        dim.7d.st <- DimPlot(st.sub.7d, group.by='Manuscript') +
            NoLegend() + ggtitle("7 Day Post-Injury")
        suppressWarnings(print(LabelClusters(dim.7d.st, id='Manuscript', size=4, repel=1, box.padding=1)))
    })

    output$st.FeatPlot.0d <- renderPlot({
        req(st.gene() %in% rownames(s.obj.0d))
        suppressWarnings(print(FeaturePlot(st.sub.0d, slot='data', features=st.gene(), order=T, min.cutoff=1) +
                               xlim(-13, 12) + ylim(-9, 11)))
    })
    output$st.FeatPlot.1d <- renderPlot({
        req(st.gene() %in% rownames(s.obj.1d))
        suppressWarnings(print(FeaturePlot(st.sub.1d, slot='data', features=st.gene(), order=T, min.cutoff=1) +
                               xlim(6,14) + ylim(-7,9)))
    })
    output$st.FeatPlot.3d <- renderPlot({
        req(st.gene() %in% rownames(s.obj.3d))
        suppressWarnings(print(FeaturePlot(st.sub.3d, slot='data', features=st.gene(), order=T, min.cutoff=1) +
                               xlim(-15, -3) + ylim(-13, 7)))
    })
    output$st.FeatPlot.7d <- renderPlot({
        req(st.gene() %in% rownames(s.obj.7d))
        suppressWarnings(print(FeaturePlot(st.sub.7d, slot='data', features=st.gene(), order=T, min.cutoff=1)))
    })
    output$TopMarks.st.0d <- renderTable({
        req(input$st.clus.0d %in% st.cluster.0d)
        st.sub.top.0d %>% ungroup() %>% filter(cluster==input$st.clus.0d) %>%
            select(-p_val, -p_val_adj, -cluster) %>% relocate(gene) %>% as.data.frame()
    })
    output$TopMarks.st.1d <- renderTable({
        req(input$st.clus.1d %in% st.cluster.1d)
        st.sub.top.1d %>% ungroup() %>% filter(cluster==input$st.clus.1d) %>%
            select(-p_val, -p_val_adj, -cluster) %>% relocate(gene) %>% as.data.frame()
    })
    output$TopMarks.st.3d <- renderTable({
        req(input$st.clus.3d %in% st.cluster.3d)
        st.sub.top.3d %>% ungroup() %>% filter(cluster==input$st.clus.3d) %>%
            select(-p_val, -p_val_adj, -cluster) %>% relocate(gene) %>% as.data.frame()
    })
    output$TopMarks.st.7d <- renderTable({
        req(input$st.clus.7d %in% st.cluster.7d)
        st.sub.top.7d %>% ungroup() %>% filter(cluster==input$st.clus.7d) %>%
            select(-p_val, -p_val_adj, -cluster) %>% relocate(gene) %>% as.data.frame()
    })
    
###################################################################################################
###### Integrated Immune Cells #####
    output$Int.DimPlot.0d <- renderPlot({
        req(exists("Int.Data"))
        x <- DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Naive'),
                     group.by='Manuscript', pt.size=.2) + NoLegend() + ggtitle("Naive Sciatic Nerve")
        x <- LabelClusters(x, id='Manuscript', size=4, repel=1, box.padding=1)
        x
    })
    output$Int.DimPlot.b <- renderPlot({
        req(exists("Int.Data"))
        x <- DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Blood'),
                     group.by='Manuscript', pt.size=.2) + NoLegend() + ggtitle("Blood")
        x <- LabelClusters(x, id='Manuscript', size=4, repel=1, box.padding=1)
        x
    })
    output$Int.DimPlot.1d <- renderPlot({
        req(exists("Int.Data"))
        x <- DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Day1'),
                     group.by='Manuscript', pt.size=.2) + NoLegend() + ggtitle("Day1 Sciatic Nerve")
        x <- LabelClusters(x, id='Manuscript', size=4, repel=1, box.padding=1)
        x
    })
    output$Int.DimPlot.3d <- renderPlot({
        req(exists("Int.Data"))
        x <- DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Day3'),
                     group.by='Manuscript', pt.size=.2) + NoLegend() + ggtitle("Day3 Sciatic Nerve")
        x <- LabelClusters(x, id='Manuscript', size=4, repel=1, box.padding=1)
        x
    })
    output$Int.DimPlot.7d <- renderPlot({
        req(exists("Int.Data"))
        x <- DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Day7'),
                     group.by='Manuscript', pt.size=.2) + NoLegend() + ggtitle("Day7 Sciatic Nerve")
        x <- LabelClusters(x, id='Manuscript', size=4, repel=1, box.padding=1)
        x
    })
    output$Int.0d.feature <- renderPlot({
        req(exists("Int.Data"))
        Data <- Int.Data()
        req(input$IntGene %in% rownames(Data))
        FeaturePlot(Int.Data(), cells=which(Data$Time=='Naive'),
                    slot='data', features=input$IntGene, order=T, min.cutoff=1)
    })
    output$Int.b.feature <- renderPlot({
        req(exists("Int.Data"))
        Data <- Int.Data()
        req(input$IntGene %in% rownames(Data))
        FeaturePlot(Data, cells=which(Data$Time=='Blood'),
                    slot='data', features=input$IntGene, order=T, min.cutoff=1)
    })
    output$Int.1d.feature <- renderPlot({
        req(exists("Int.Data"))
        Data <- Int.Data()
        req(input$IntGene %in% rownames(Data))
        FeaturePlot(Data, cells=which(Data$Time=='Day1'),
                    slot='data', features=input$IntGene, order=T, min.cutoff=1)
    })
    output$Int.3d.feature <- renderPlot({
        req(exists("Int.Data"))
        Data <- Int.Data()
        req(input$IntGene %in% rownames(Data))
        FeaturePlot(Data, cells=which(Data$Time=='Day3'),
                    slot='data', features=input$IntGene, order=T, min.cutoff=1)
    })
    output$Int.7d.feature <- renderPlot({
        req(exists("Int.Data"))
        Data <- Int.Data()
        req(input$IntGene %in% rownames(Data))
        FeaturePlot(Data, cells=which(Data$Time=='Day7'),
                    slot='data', features=input$IntGene, order=T, min.cutoff=1)
    })
    output$Int.0d.seurat <- renderPlot({
        req(exists("Int.Data"))
        DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Naive'),
                group.by='seurat_clusters', pt.size=.2, label=T) +
            NoLegend() + ggtitle("Naive Sciatic Nerve")
    })
    output$Int.b.seurat <- renderPlot({
        req(exists("Int.Data"))
        DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Blood'),
                group.by='seurat_clusters', pt.size=.2, label=T) +
            NoLegend() + ggtitle("Blood")
    })
    output$Int.1d.seurat <- renderPlot({
        req(exists("Int.Data"))
        DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Day1'),
                group.by='seurat_clusters', pt.size=.2, label=T) +
            NoLegend() + ggtitle("Day1 Sciatic Nerve")
    })
    output$Int.3d.seurat <- renderPlot({
        req(exists("Int.Data"))
        DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Day3'),
                group.by='seurat_clusters', pt.size=.2, label=T) +
            NoLegend() + ggtitle("Day3 Sciatic Nerve")
    })
    output$Int.7d.seurat <- renderPlot({
        req(exists("Int.Data"))
        DimPlot(Int.Data(), cells=which(Int.Data()$Time=='Day7'),
                group.by='seurat_clusters', pt.size=.2, label=T) +
            NoLegend() + ggtitle("Day7 Sciatic Nerve")
    })
    int.seurat.markers <- reactive({
        top.int %>%
            filter(cluster==input$Int.Cluster) %>%
            filter(p_val_adj < .0001) %>%
            select(-p_val, -cluster) %>%
            mutate(across(where(is.numeric), round, 3)) %>%
            relocate('gene')
    })
    output$Int_Seurat_Markers <- renderDT(
        int.seurat.markers(),
        class="display nowrap compact",
        filter="top",
        rownames=FALSE,
        extensions = 'Buttons',
        options=list(order = list(list(4,'asc')),
                     dom = 'lfrtipB',
                     buttons = list(
                         'copy',
                         list(extend='csv',
                              filename=paste0('Positively_Regulated_Genes_For_Cluster_',
                                              input$Int.Cluster,
                                              '_in_Immune_Integrated_Data'))
                     )
                     )
    )
    
### Trajectory
    Traj.Data <- eventReactive(input$LoadTrajectory, {
        w.load.Traj$show()
        sds <- readRDS('sds.RDS')
        var.imp.1 <- readRDS('Variable_Importance_Mo_DC.RDS')
        var.imp.2 <- readRDS('Variable_Importance_Mo_Mac.RDS')
        w.load.CC$hide()
        list('sds'=sds, 'Mo_DC_var_imp'=var.imp.1, 'Mo_Mac_var_imp'=var.imp.2)
    })
    var_imp_table <- reactive({
        req(exists("Traj.Data"))
        mac <- Traj.Data()[['Mo_Mac_var_imp']]
        dc <- Traj.Data()[['Mo_DC_var_imp']]
        ix <- match(names(mac), names(dc))
        data.frame('Gene'=names(mac), 'Mo To Mac Importance Rank'=mac, 'Mo To DC Importance Rank'=dc[ix],
                   check.names=FALSE)
    })
    output$Var_Imp <- renderDT(
        var_imp_table(),
        class="display nowrap compact",
        filter="top",
        rownames=FALSE,
        extensions = 'Buttons',
        options=list(order = list(list(1,'asc')),
                     dom = 'lfrtipB',
                     buttons = list(
                         'copy',
                         list(extend='csv',
                              filename='Trajectory_Importance'))
                     )
    )
    output$Traj_PCA <- renderPlot({
        req(exists("Int.Data"))
        p <- DimPlot(Int.Data(), group.by=input$LabelBy, label=F, reduction='pca') + NoLegend()
        LabelClusters(p, id=input$LabelBy, size=4, repel=2, box.padding=.2, max.overlaps=20)
    })
    output$Traj_UMAP <- renderPlot({
        req(exists("Int.Data"))
        p <- DimPlot(Int.Data(), group.by=input$LabelBy, label=F) + NoLegend()
        LabelClusters(p, id=input$LabelBy, size=5, repel=1, box.padding=1)
    })
    output$Traj_PCA_pseudo <- renderPlot({
        req(exists("Int.Data"))
        req(exists("Traj.Data"))
        Data <- Int.Data()
        sds <- Traj.Data()[["sds"]]
        curve <- NULL
        if(input$Traj=='Mo_to_Mac') curve=9
        if(input$Traj=='Mo_to_DC') curve=6
        pseudotime.line <- slingPseudotime(SlingshotDataSet(sds))[,curve]
        pca.emb <- Embeddings(Data, "pca")[,1:2]
        pca.plot <- data.frame('PC_1'=pca.emb[,1], 'PC_2'=pca.emb[,2], 'pseudotime'=pseudotime.line)
        pca.layer1 <- pca.plot[which(is.na(pca.plot$pseudotime)),]
        pca.layer2 <- pca.plot[which(!is.na(pca.plot$pseudotime)),]
        curves <- slingCurves(sds, as.df = TRUE) %>% filter(Lineage==curve)
        ggplot() +
            geom_point(data=pca.layer1, aes(x=PC_1, y=PC_2), color='grey75') +
            geom_point(data=pca.layer2, aes(x=PC_1, y=PC_2, color=pseudotime)) + scale_color_viridis() +
            geom_path(data = curves %>% arrange(Order), aes(x=PC_1, y=PC_2, group = Lineage))
    })
    output$Traj_UMAP_pseudo <- renderPlot({
        req(exists("Int.Data"))
        req(exists("Traj.Data"))
        Data <- Int.Data()
        sds <- Traj.Data()[["sds"]]
        curve <- NULL
        if(input$Traj=='Mo_to_Mac') curve=9
        if(input$Traj=='Mo_to_DC') curve=6
        pseudotime.line <- slingPseudotime(SlingshotDataSet(sds))[,curve]
        umap <- Embeddings(Data, "umap")
        umap.plot <- data.frame('UMAP_1'=umap[,1], 'UMAP_2'=umap[,2], 'pseudotime'=pseudotime.line)
        umap.layer1 <- umap.plot[which(is.na(umap.plot$pseudotime)),]
        umap.layer2 <- umap.plot[which(!is.na(umap.plot$pseudotime)),]
        ggplot() +
            geom_point(data=umap.layer1, aes(x=UMAP_1, y=UMAP_2), color='grey75') +
            geom_point(data=umap.layer2, aes(x=UMAP_1, y=UMAP_2, color=pseudotime)) + scale_color_viridis()
    })
    output$Traj_PCA_Feature <- renderPlot({
        req(exists("Int.Data"))
        Data <- Int.Data()
        req(input$TrajGene %in% rownames(Data))
        FeaturePlot(Data, features=input$TrajGene, min.cutoff=TRUE, order=TRUE, reduction='pca')
    })
    output$Traj_UMAP_Feature <- renderPlot({
        req(exists("Int.Data"))
        Data <- Int.Data()
        req(input$TrajGene %in% rownames(Data))
        FeaturePlot(Data, features=input$TrajGene, min.cutoff=TRUE, order=TRUE)
    })
    
    ## CellSignal
    output$L_R_Day_0 <- renderTable({
        req(exists("CC.Data"))
        df.net.0d <- CC.Data()[['df.0d']]
        df.net.0d %>% filter(source == input$source.0d) %>% arrange(desc(prob)) %>%
            select(c(target, ligand, receptor, prob, pathway_name)) %>%
            mutate(across(where(is.numeric), round, 3)) %>% .[1:25,]
    })
    output$L_R_Day_1 <- renderTable({
        req(exists("CC.Data"))
        df.net.1d <- CC.Data()[['df.1d']]
        df.net.1d %>% filter(source == input$source.1d) %>% arrange(desc(prob)) %>%
            select(c(target, ligand, receptor, prob, pathway_name)) %>%
            mutate(across(where(is.numeric), round, 3)) %>% .[1:25,]
    })
    output$L_R_Day_3 <- renderTable({
        req(exists("CC.Data"))
        df.net.3d <- CC.Data()[['df.3d']]
        df.net.3d %>% filter(source == input$source.3d) %>% arrange(desc(prob)) %>%
            select(c(target, ligand, receptor, prob, pathway_name)) %>%
            mutate(across(where(is.numeric), round, 3)) %>% .[1:25,]
    })
    output$L_R_Day_7 <- renderTable({
        req(exists("CC.Data"))
        df.net.7d <- CC.Data()[['df.7d']]
        df.net.7d %>% filter(source == input$source.7d) %>% arrange(desc(prob)) %>%
            select(c(target, ligand, receptor, prob, pathway_name)) %>%
            mutate(across(where(is.numeric), round, 3)) %>% .[1:25,]
    })
    output$NetVisCircle.0d <- renderPlot({
        req(exists("CC.Data"))
        CellChat.0d <- CC.Data()[['CC.0d']]
        groupSize <- as.numeric(table(CellChat.0d@idents))
        mat <- CellChat.0d@net$weight
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        ix <- which(levels(CellChat.0d@idents) == input$source.0d)
        mat2[ix, ] <- mat[ix, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
                         edge.weight.max = max(mat), title.name = rownames(mat)[ix])
    })
    output$NetVisCircle.1d <- renderPlot({
        req(exists("CC.Data"))
        CellChat.1d <- CC.Data()[['CC.1d']]
        groupSize <- as.numeric(table(CellChat.1d@idents))
        mat <- CellChat.1d@net$weight
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        ix <- which(levels(CellChat.1d@idents) == input$source.1d)
        mat2[ix, ] <- mat[ix, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
                         edge.weight.max = max(mat), title.name = rownames(mat)[ix])
    })
    output$NetVisCircle.3d <- renderPlot({
        req(exists("CC.Data"))
        CellChat.3d <- CC.Data()[['CC.3d']]
        groupSize <- as.numeric(table(CellChat.3d@idents))
        mat <- CellChat.3d@net$weight
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        ix <- which(levels(CellChat.3d@idents) == input$source.3d)
        mat2[ix, ] <- mat[ix, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
                         edge.weight.max = max(mat), title.name = rownames(mat)[ix])
    })
    output$NetVisCircle.7d <- renderPlot({
        req(exists("CC.Data"))
        CellChat.7d <- CC.Data()[['CC.7d']]
        groupSize <- as.numeric(table(CellChat.7d@idents))
        mat <- CellChat.7d@net$weight
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        ix <- which(levels(CellChat.7d@idents) == input$source.7d)
        mat2[ix, ] <- mat[ix, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
                         edge.weight.max = max(mat), title.name = rownames(mat)[ix])
    })

    ## Error catching
    cc.err.0d <- reactive({
        req(exists("CC.Data"))
        CellChat.0d <- CC.Data()[['CC.0d']]
        temp <- tryCatch(netAnalysis_contribution(CellChat.0d, signaling = input$signal.path),
                         error = function(e){e$message})
        list(is.character(temp), temp)
    })
    cc.err.1d <- reactive({
        req(exists("CC.Data"))
        CellChat.1d <- CC.Data()[['CC.1d']]
        temp <- tryCatch(netAnalysis_contribution(CellChat.1d, signaling = input$signal.path),
                         error = function(e){e$message})
        list(is.character(temp), temp)
    })
    cc.err.3d <- reactive({
        req(exists("CC.Data"))
        CellChat.3d <- CC.Data()[['CC.3d']]
        temp <- tryCatch(netAnalysis_contribution(CellChat.3d, signaling = input$signal.path),
                         error = function(e){e$message})
        list(is.character(temp), temp)
    })
    cc.err.7d <- reactive({
        req(exists("CC.Data"))
        CellChat.7d <- CC.Data()[['CC.7d']]
        temp <- tryCatch(netAnalysis_contribution(CellChat.7d, signaling = input$signal.path),
                         error = function(e){e$message})
        list(is.character(temp), temp)
    })

    ## SigPath
    output$heatmap.path.0d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.0d()[[1]] == FALSE){
            pathways.0d <- CC.Data()[['path.0d']]
            CellChat.0d <- CC.Data()[['CC.0d']]
            print(netVisual_heatmap(CellChat.0d, signaling = input$signal.path, color.heatmap = "Reds"))
        } else {
            validate(
                need(!is.character(cc.err.0d()[[2]]), paste(cc.err.0d()[[2]]))
            )
            print(cc.err.0d()[[2]])
        }
    })
    output$heatmap.path.1d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.1d()[[1]] == FALSE){
            pathways.1d <- CC.Data()[['path.1d']]
            CellChat.1d <- CC.Data()[['CC.1d']]
            print(netVisual_heatmap(CellChat.1d, signaling = input$signal.path, color.heatmap = "Reds"))
        } else {
            validate(
                need(!is.character(cc.err.1d()[[2]]), paste(cc.err.1d()[[2]]))
            )
            print(cc.err.1d()[[2]])
        }
    })
    output$heatmap.path.3d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.3d()[[1]] == FALSE){
            pathways.3d <- CC.Data()[['path.3d']]
            CellChat.3d <- CC.Data()[['CC.3d']]
            print(netVisual_heatmap(CellChat.3d, signaling = input$signal.path, color.heatmap = "Reds"))
        } else {
            validate(
                need(!is.character(cc.err.3d()[[2]]), paste(cc.err.3d()[[2]]))
            )
            print(cc.err.3d()[[2]])
        }
    })
    output$heatmap.path.7d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.7d()[[1]] == FALSE){
            pathways.7d <- CC.Data()[['path.7d']]
            CellChat.7d <- CC.Data()[['CC.7d']]
            print(netVisual_heatmap(CellChat.7d, signaling = input$signal.path, color.heatmap = "Reds"))
        } else {
            validate(
                need(!is.character(cc.err.7d()[[2]]), paste(cc.err.7d()[[2]]))
            )
            print(cc.err.7d()[[2]])
        }
    })
    output$chord.path.0d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.0d()[[1]] == FALSE){
            CellChat.0d <- CC.Data()[['CC.0d']]
            group.cellType <- c(rep("MAC", 3), rep('IM',3), rep("Str", 7))
            names(group.cellType) <- levels(CellChat.0d@idents)
            print(netVisual_chord_cell(CellChat.0d, signaling = input$signal.path, group = group.cellType,
                                       title.name = paste0(input$signal.path, " signaling network")))
        } else {
            validate(
                need(!is.character(cc.err.0d()[[2]]), paste(cc.err.0d()[[2]]))
            )
            print(cc.err.0d()[[2]])
        }
    })
    output$chord.path.1d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.1d()[[1]] == FALSE){
            CellChat.1d <- CC.Data()[['CC.1d']]
            group.cellType <- c(rep("Mac", 6), rep("IM", 8), rep('Str', 7))
            names(group.cellType) <- levels(CellChat.1d@idents)
            print(netVisual_chord_cell(CellChat.1d, signaling = input$signal.path, group = group.cellType,
                                       title.name = paste0(input$signal.path, " signaling network")))
        } else {
            validate(
                need(!is.character(cc.err.1d()[[2]]), paste(cc.err.1d()[[2]]))
            )
            print(cc.err.1d()[[2]])
        }
    })
    output$chord.path.3d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.3d()[[1]] == FALSE){
            CellChat.3d <- CC.Data()[['CC.3d']]
            group.cellType <- c(rep("Mac", 5), rep("IM", 7), rep("Str", 7))
            names(group.cellType) <- levels(CellChat.3d@idents)
            print(netVisual_chord_cell(CellChat.3d, signaling = input$signal.path, group = group.cellType,
                                       title.name = paste0(input$signal.path, " signaling network")))
        } else {
            validate(
                need(!is.character(cc.err.3d()[[2]]), paste(cc.err.3d()[[2]]))
            )
            print(cc.err.3d()[[2]])
        }
    })
    output$chord.path.7d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.7d()[[1]] == FALSE){
            CellChat.7d <- CC.Data()[['CC.7d']]
            group.cellType <- c(rep("Mac", 6), rep("IM", 8), rep("Str", 8))
            names(group.cellType) <- levels(CellChat.7d@idents)
            print(netVisual_chord_cell(CellChat.7d, signaling = input$signal.path, group = group.cellType,
                                       title.name = paste0(input$signal.path, " signaling network")))
        } else {
            validate(
                need(!is.character(cc.err.7d()[[2]]), paste(cc.err.7d()[[2]]))
            )
            print(cc.err.7d()[[2]])
        }
    })
    output$contrib.path.0d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.0d()[[1]] == FALSE){
            print(cc.err.0d()[[2]])
        } else {
            validate(
                need(!is.character(cc.err.0d()[[2]]), paste(cc.err.0d()[[2]]))
            )
            print(cc.err.0d()[[2]])
        }
    })
    output$contrib.path.1d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.1d()[[1]] == FALSE){
            print(cc.err.1d()[[2]])
        } else {
            validate(
                need(!is.character(cc.err.1d()[[2]]), paste(cc.err.1d()[[2]]))
            )
            print(cc.err.1d()[[2]])
        }
    })
    output$contrib.path.3d <- renderPlot({
        req(exists("CC.Data"))
        if(cc.err.3d()[[1]] == FALSE){
            print(cc.err.3d()[[2]])
        } else {
            validate(
                need(!is.character(cc.err.3d()[[2]]), paste(cc.err.3d()[[2]]))
            )
            print(cc.err.3d()[[2]])
        }
    })
    output$contrib.path.7d <- renderPlot({
        if(cc.err.7d()[[1]] == FALSE){
            print(cc.err.7d()[[2]])
        } else {
            validate(
                need(!is.character(cc.err.7d()[[2]]), paste(cc.err.7d()[[2]]))
            )
            print(cc.err.7d()[[2]])
        }
    })

   values <- reactive({
        req(exists("CC.Data"))
        CC <- list(CC.Data()[['CC.0d']], CC.Data()[['CC.1d']],
                   CC.Data()[['CC.3d']], CC.Data()[['CC.7d']])
        ## lot of this from https://rdrr.io/github/sqjin/CellChat/src/R/visualization.R
        ## want to only give L-R pairs that pass Cell Chat filters for plot
        output.list <- list()
        for(i in 1:4){
            ## Get the pairs for the selected signaling pathway
            pairLR <- searchPair(signaling = input$signal.path,
                                 pairLR.use = CC[[i]]@LR$LRsig,
                                 key = "pathway_name", matching.exact = T, pair.only = F)
            pairLR.use.name <- dimnames(CC[[i]]@net$prob)[[3]]
            pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)

            ## Get the probabilities to determine which are significant enough to plot
            prob <- CC[[i]]@net$prob
            pval <- CC[[i]]@net$pval
            prob[pval > .05] <- 0

            ## this is the netVisual_individual code, if length of pairLR.name.use > 0 it will plot
            if (length(pairLR.name) > 1) {
                pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
            } else {
                pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
            }
            output.list[[i]] <- pairLR.name.use
        }
        output.list
    })
    observe({
        updateSelectizeInput(session, 'LR_Pair_0d', choices = values()[[1]], server = TRUE)
        updateSelectizeInput(session, 'LR_Pair_1d', choices = values()[[2]], server = TRUE)
        updateSelectizeInput(session, 'LR_Pair_3d', choices = values()[[3]], server = TRUE)
        updateSelectizeInput(session, 'LR_Pair_7d', choices = values()[[4]], server = TRUE)
    })
    check.pair <- function(object, signal){
        pairLR <- searchPair(signaling = input$signal.path, pairLR.use = object@LR$LRsig,
                             key = "pathway_name", matching.exact = T, pair.only = F)
        pairLR.use.name <- dimnames(object@net$prob)[[3]]
        pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
        pairLR.name <- intersect(pairLR.name, as.character(signal))
        length(pairLR.name)
    }
    output$LR_Pair.0d <- renderPlot({
        req(input$LR_Pair_0d)
        CellChat.0d <- CC.Data()[['CC.0d']]
        req(check.pair(CellChat.0d, input$LR_Pair_0d) > 0)
        group.cellType <- c(rep("MAC", 3), rep('IM',3), rep("Str", 7))
        names(group.cellType) <- levels(CellChat.0d@idents)
        netVisual_individual(CellChat.0d, signaling = input$signal.path, pairLR.use = input$LR_Pair_0d,
                             layout = "chord", vertex.weight=NULL, group = group.cellType,)
    })
    output$LR_Pair.1d <- renderPlot({
        CellChat.1d <- CC.Data()[['CC.1d']]
        req(check.pair(CellChat.1d, input$LR_Pair_1d) > 0)
        group.cellType <- c(rep("Mac", 6), rep("IM", 8), rep('Str', 7))
        names(group.cellType) <- levels(CellChat.1d@idents)
        netVisual_individual(CellChat.1d, signaling = input$signal.path, pairLR.use = input$LR_Pair_1d,
                             layout = "chord", vertex.weight=NULL, group = group.cellType,)
    })
    output$LR_Pair.3d <- renderPlot({
        req(exists("CC.Data"))
        CellChat.3d <- CC.Data()[['CC.3d']]
        req(check.pair(CellChat.3d, input$LR_Pair_3d) > 0)
        group.cellType <- c(rep("Mac", 5), rep("IM", 7), rep("Str", 7))
        names(group.cellType) <- levels(CellChat.3d@idents)
        netVisual_individual(CellChat.3d, signaling = input$signal.path, pairLR.use = input$LR_Pair_3d,
                             layout = "chord", vertex.weight=NULL, group = group.cellType,)
    })
    output$LR_Pair.7d <- renderPlot({
        CellChat.7d <- CC.Data()[['CC.7d']]
        req(check.pair(CellChat.7d, input$LR_Pair_7d) > 0)
        group.cellType <- c(rep("Mac", 6), rep("IM", 8), rep("Str", 8))
        names(group.cellType) <- levels(CellChat.7d@idents)
        netVisual_individual(CellChat.7d, signaling = input$signal.path, pairLR.use = input$LR_Pair_7d,
                             layout = "chord", vertex.weight=NULL, group = group.cellType,)
    })

    ## Spatial Distribution
    output$Location.UMAP <- renderPlot({
        l.umap <- DimPlot(s.obj.loc, group.by='CellType') +  NoLegend() + ggtitle("All Cells")
        suppressWarnings(print(LabelClusters(l.umap, id='CellType', size=5, repel=1, box.padding=2)))
    })
    output$Location.UMAP.Split <- renderPlot({
        l.umap.s <- DimPlot(s.obj.loc, group.by='CellType', split.by='Location') +  NoLegend() +
            ggtitle("Cells Split By Location")
        l.umap.s
    })
    output$Loc.Cluster.Size <- renderTable({
        dis.clus.size <- table(s.obj.loc$CellType[which(s.obj.loc$Location=='Distal')])
        inj.clus.size <- table(s.obj.loc$CellType[which(s.obj.loc$Location=='Injury')])
        tab <- data.frame(cbind(names(dis.clus.size), dis.clus.size, inj.clus.size), row.names=NULL)
        colnames(tab) <- c('Cell Type','Distal Cells','Injury Cells')
        tab
    })
    output$loc.markers <- renderTable({
        top.loc %>% filter(cluster == input$loc.clus) %>%
            select(-cluster) %>% `rownames<-`( NULL ) %>% .[1:35,]
    })
    output$DE.distal <- renderTable({
        distal.DE[[input$loc.clus]] %>% select(-p_val) %>% rownames_to_column(var='gene') %>%
            relocate(gene) %>% .[1:35,]
    })
    output$DE.injury <- renderTable({
        injury.DE[[input$loc.clus]] %>% select(-p_val) %>% rownames_to_column(var='gene') %>%
            relocate(gene) %>% .[1:35,]
    })
    output$Location.Feature <- renderPlot({
        req(loc.gene() %in% rownames(s.obj.loc))
        print(FeaturePlot(s.obj.loc, features=loc.gene(), order=T, min.cutoff=1))
    })
    output$Location.Feature.Split <- renderPlot({
        req(loc.gene() %in% rownames(s.obj.loc))
        print(FeaturePlot(s.obj.loc, features=loc.gene(), split.by='Location',
                          order=T, min.cutoff=1))
    })
    
    ## PBMC
    clusters.blood <- reactive({
        req(exists("PBMC.Data"))
        blood <- PBMC.Data()
        levels(blood$Manuscript)
    })
    observe({
        updateSelectizeInput(session, 'Cluster.blood', choices = clusters.blood(), server = TRUE)
    })
    genes.blood <- reactive({
        req(exists("PBMC.Data"))
        blood <- PBMC.Data()
        rownames(blood)
    })
    observe({
        updateSelectizeInput(session, 'Gene.pbmc', choices = genes.blood(), selected='Ccr2', server = TRUE)
        updateSelectizeInput(session, 'Gene2.pbmc', choices = genes.blood(), selected='Apoe', server = TRUE)
    })
    output$AllDimPlot.blood <- renderPlot({
        x <- DimPlot(PBMC.Data(), group.by='Manuscript', pt.size=.2) +
            NoLegend() + ggtitle("PBMC UMAP")
        LabelClusters(x, id='Manuscript', size=4, repel=1, box.padding=1)
    })
    output$AllFeaturePlot.blood <- renderPlot({
        req(input$Gene.pbmc)
        print(FeaturePlot(PBMC.Data(), slot='data', features=input$Gene.pbmc, order=T, min.cutoff=1))
    })
    output$Blend.blood <- renderPlot({
        req(input$Gene.pbmc)
        req(input$Gene2.pbmc)
        FeaturePlot(PBMC.Data(), min.cutoff=1, order=T, slot='data',
                    features=c(input$Gene.pbmc, input$Gene2.pbmc), blend=T, pt.size=0.2)
    })
    output$Fscatter.blood <- renderPlot({
        req(input$Gene.pbmc)
        req(input$Gene2.pbmc)
        FeatureScatter(PBMC.Data(), feature1=input$Gene.pbmc, feature2=input$Gene2.pbmc,
                       group.by='Manuscript', plot.cor=T, pt.size=0.2)
    })
    output$blood.markers <- renderTable({
        req(input$Cluster.blood)
        top.blood %>% filter(cluster == input$Cluster.blood) %>% select(-p_val) %>%
            relocate(gene) %>% .[1:50,]
    })
}

shinyApp(ui, server)

