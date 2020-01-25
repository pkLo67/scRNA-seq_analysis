#-------------------------------------------------------------------#
## Shiny App for Analysis of Single-Cell RNA-Seq Data ##
#-------------------------------------------------------------------#


library(shiny)
library(Seurat)
library(dplyr)
library(DT)

## Paste the File Path (e.g. "/Users/username/Desktop/folder1/folder2/folder3") into the App 
## before you click on "Start Analysis"
## The last folder should contain matrix, index and gene feature files.
## Enter "gene names" to display the Gene Feature Plot and Gene VlnPlot.
## Annotated plots can be displayed after cell cluster IDs are renamed.
## Be patient in waiting for the output results as some code execution need longer computing time.


# Define UI for application that generates analyzed scRNA-seq results
ui <- fluidPage(
  
  # Application title
  titlePanel("Single-Cell RNA-Seq Data Analysis"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      actionButton(
        inputId = "start", label = strong("Start Analysis")
      ),
      textInput(
        inputId = "file_path", label = strong("Paste the file path"), value = ""
      ),
      textInput(
        inputId = "project", label = strong("Enter the project name"), value = "scRNA-seq"
      ),
      selectInput(
        inputId = "plot", label = strong("Select a plot to display"),
        choices = c("VlnPlot", "FeatureScatter Plots", "VariableFeaturePlot", "VizDimLoadings Plots",
                    "ElbowPlot")
      ),
      sliderInput(
        inputId = "features_n", label = strong("Set the number of gene symbols showing in VariableFeaturePlot"),
        min = 1, max = 50, value = 10
      ),
      numericInput(
        inputId = "dims_n1", label = strong("Set the number of dims for VizDimLoadings Plots"),
        min = 1, max = NA, value = 2
      ),
      numericInput(
        inputId = "dims_n2", label = strong("Set the number of dims for DimHeatmap"),
        min = 1, max = NA, value = 6
      ),
      sliderInput(
        inputId = "pt_size", label = strong("Set the point size for PCA, t-SNE and UMAP plots"),
        min = 0.5, max = 1.5, value = 0.8
      ),
      numericInput(
        inputId = "dims_n3", label = strong("Set the number of dims for t-SNE and UMAP plots"),
        min = 3, max = NA, value = 15
      ),
      selectInput(
        inputId = "idn", label = strong("Select a cell cluster for gene marker identification"),
        choices = c(0:25), selected = 0
      ),
      numericInput(
        inputId = "marker_n", label = strong("Set the number of gene markers displayed in the table"),
        min = 10, max = NA, value = 20
      ),
      tags$hr(),
      textInput(
        inputId = "genes", label = strong("Enter gene names for the Gene Feature Plot (name1,name2,.., no space)"), value = ""
      ),
      actionButton(
        inputId = "show", label = strong("Show Plot")
      ),
      tags$hr(),
      textInput(
        inputId = "rename", label = strong("Rename cell cluster IDs (e.g. 0:CD4+,1:CD8+,.., no space)"), value = ""
      ),
      actionButton(
        inputId = "show2", label = strong("Show Annotated Plot")
      ),
      tags$hr(),
      numericInput(
        inputId = "label_size", label = strong("Set the size of cell cluster labels"),
        min = 4, max = NA, value = 8
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("QC & DIM plots",
                 plotOutput("plot")
        ),
        tabPanel("DimHeatmap",
                 plotOutput("DimHeatmap", height = "500px")
        ),
        tabPanel("PCA Plot",
                 plotOutput("pca", width = "70%", height = "500px")       
        ),
        tabPanel("t-SNE Plot",
                 plotOutput("tsne", width = "70%", height = "500px")
        ),
        tabPanel("UMAP Plot",
                 plotOutput("umap", width = "70%", height = "500px")
        ),
        tabPanel("Cluster-specific gene markers",
                 DT::DTOutput("marker_table")
        ),
        tabPanel("Gene Feature Plot",
                 plotOutput("gene_plot", width = "100%", height = "400px")
        ),
        tabPanel("Gene VlnPlot",
                 plotOutput("gene_vlnplot", width = "100%", height = "400px")
        ),
        tabPanel("Annotated t-SNE Plot",
                 plotOutput("tsne_an", width = "70%", height = "500px")
        ),
        tabPanel("Annotated UMAP Plot",
                 plotOutput("umap_an", width = "70%", height = "500px")
        ),
        tabPanel("Annotated Gene VlnPlot",
                 plotOutput("gene_vlnplot_an", width = "100%", height = "400px")
        )
      )
    )
  )
)

# Define server logic required to generate plots and tables
server <- function(input, output) {
  # Paste the file path into the entry block and click the "Start Analysis" bottom to start analysis.
  # The last folder should contain matrix, index and gene feature files.
  rawdata <- eventReactive(input$start, {
    Read10X(data.dir = input$file_path, gene.column = 2)
  })
  
  seuobj <- reactive({
    seuobj <- CreateSeuratObject(counts = rawdata(), project = input$project, min.cells = 3, min.features = 200)
    seuobj[["percent.mt"]] <- PercentageFeatureSet(seuobj, pattern = "^MT-")
    seuobj
  })   
  
  seuobj2 <- reactive({
    seuobj() %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  }) 
  
  seuobj3 <- reactive({
    seuobj2() %>% ScaleData() %>% 
      RunPCA(features = VariableFeatures(object = seuobj2()))
  })
  
  seuobj4 <- reactive({
    seuobj3() %>% FindNeighbors(dims = 1:input$dims_n3) %>%
      FindClusters(resolution = 0.5) %>%
      RunTSNE(dims = 1:input$dims_n3) %>% 
      RunTSNE(dims = 1:input$dims_n3) %>%
      RunUMAP(dims = 1:input$dims_n3)
  })  
  
  cluster_ids <- reactive({
    toString(levels(Idents(seuobj4())))  # cell cluster ids
  })
  
  cluster_n <- reactive({
    length(cluster_ids())   # The number of cell clusters
  })
  
  cidn <- reactive({
    if (input$idn <= cluster_n()) {
      input$idn
    } else {0}
  })
  
  cell_type_markers <- eventReactive(input$show, {
    unlist(strsplit(input$genes, "[,]"))
  })
  
  seuobj5 <- eventReactive(input$show2, {
    obj <- seuobj4()
    name <- unlist(strsplit(input$rename, "[:,]"))
    for (i in seq(1, length(name), 2)) {
      if (name[i]=='0') {obj <- RenameIdents(obj, '0'=name[i+1])}
      else if (name[i]=='1') {obj <- RenameIdents(obj, '1'=name[i+1])}
      else if (name[i]=='2') {obj <- RenameIdents(obj, '2'=name[i+1])}
      else if (name[i]=='3') {obj <- RenameIdents(obj, '3'=name[i+1])}
      else if (name[i]=='4') {obj <- RenameIdents(obj, '4'=name[i+1])}
      else if (name[i]=='5') {obj <- RenameIdents(obj, '5'=name[i+1])}
      else if (name[i]=='6') {obj <- RenameIdents(obj, '6'=name[i+1])}
      else if (name[i]=='7') {obj <- RenameIdents(obj, '7'=name[i+1])}
      else if (name[i]=='8') {obj <- RenameIdents(obj, '8'=name[i+1])}
      else if (name[i]=='9') {obj <- RenameIdents(obj, '9'=name[i+1])}
      else if (name[i]=='10') {obj <- RenameIdents(obj, '10'=name[i+1])}
      else if (name[i]=='11') {obj <- RenameIdents(obj, '11'=name[i+1])}
      else if (name[i]=='12') {obj <- RenameIdents(obj, '12'=name[i+1])}
      else if (name[i]=='13') {obj <- RenameIdents(obj, '13'=name[i+1])}
      else if (name[i]=='14') {obj <- RenameIdents(obj, '14'=name[i+1])}
      else if (name[i]=='15') {obj <- RenameIdents(obj, '15'=name[i+1])}
      else if (name[i]=='16') {obj <- RenameIdents(obj, '16'=name[i+1])}
      else if (name[i]=='17') {obj <- RenameIdents(obj, '17'=name[i+1])}
      else if (name[i]=='18') {obj <- RenameIdents(obj, '18'=name[i+1])}
    }
    obj
  })
  
  output$plot <- renderPlot({
    if (input$plot=="VlnPlot") {
      VlnPlot(seuobj(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0.3) ### VlnPlot
    } else if (input$plot=="FeatureScatter Plots") {
      plot1 <- FeatureScatter(seuobj(), feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(seuobj(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      CombinePlots(plots = list(plot1, plot2))  ### FeatureScatter Plots
    } else if (input$plot=="VariableFeaturePlot") {
      top10 <- head(VariableFeatures(seuobj2()), input$features_n)
      # Perform feature plot analysis of variable expressed genes
      plot1 <- VariableFeaturePlot(seuobj2())
      LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
    } else if (input$plot=="VizDimLoadings Plots") {
      VizDimLoadings(seuobj3(), dims = 1:input$dims_n1, reduction = "pca")   ### VizDimLoadings Plots
    } else {
      ElbowPlot(seuobj3())    ### ElbowPlot
    }
  })
  
  output$DimHeatmap <-renderPlot({
    DimHeatmap(seuobj3(), dims = 1:input$dims_n2, cells = 500, balanced = TRUE)   ### DimHeatmap
  })
  
  output$pca <- renderPlot({
    DimPlot(seuobj3(), reduction = "pca", pt.size = input$pt_size)   ### PCA Plot
  })
  
  output$tsne <- renderPlot({
    DimPlot(seuobj4(), reduction = "tsne", pt.size = input$pt_size)  ### t-SNE Plot
  })
  
  output$umap <- renderPlot({
    DimPlot(seuobj4(), reduction = "umap", pt.size = input$pt_size)    ### UMAP Plot
  })
  
  output$marker_table <- DT::renderDT({
    head(FindMarkers(seuobj4(), ident.1 = cidn(), min.pct = 0.25, logfc.threshold = 0.25), input$marker_n)
  })
  
  output$gene_plot <- renderPlot({
    FeaturePlot(seuobj4(), features = cell_type_markers(), pt.size = input$pt_size)  ### Gene Feature Plot
  })
  
  output$gene_vlnplot <- renderPlot({
    VlnPlot(seuobj4(), features = cell_type_markers(), slot = "counts", log = TRUE)  ### Gene VlnPlot
  })
  
  output$tsne_an <- renderPlot({
    DimPlot(seuobj5(), reduction = "tsne", pt.size = input$pt_size, label = TRUE, label.size = input$label_size)    ### Annotated t-SNE Plot
  })
  
  output$umap_an <- renderPlot({
    DimPlot(seuobj5(), reduction = "umap", pt.size = input$pt_size, label = TRUE, label.size = input$label_size)    ### Annotated UMAP Plot
  })
  
  output$gene_vlnplot_an <- renderPlot({
    VlnPlot(seuobj5(), features = cell_type_markers(), slot = "counts", log = TRUE)  ### Annotated Gene VlnPlot
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
