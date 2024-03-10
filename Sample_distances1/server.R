####library R packages####
library(shiny)
suppressMessages(library(DESeq2))
suppressMessages(library(DT))
suppressMessages(library(vsn))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(PoiClaClu))
suppressMessages(library(glmpca))

####File's size####
options(shiny.maxRequestSize = 90*1024^2)

####Server####
vals <- reactiveValues(count = 0)
shinyServer(function(input, output, session){
  #Increment the number of sessions when one is opened.
  #We use isolate() here to:
  # a) Provide a reactive context
  # b) Ensure that this expression doesn't take a reactive dependency on vals$count
  #    -- if it did, every time vals$count changed , this expression would run, leading to an infinite loop.
  isolate(vals$count <- vals$count + 1)
  
  #----RDS data input----#
  dds_upload <- reactive({
    file_in <- input$upload
    #Warning: while file is not unloaded.
    if(is.null(input$upload)){
      return(data.frame(x = "Click 'Browse...' to select or drop a file onto 'Browse' button"))
    } else {
      dds <- readRDS(file = file_in$datapath)
      #whether the input file is a DESeqDataSet.
      if("DESeqDataSet" %in% class(dds)){
        return(DESeq(dds))
      } else {
        return(data.frame(x = "The input data is not a 'DESeqDataset'"))
      }
    }
  })
  
  #----Display Uploaded Sample Table----#
  df_upload <- reactive({
    data <- dds_upload()
    if("DESeqDataSet" %in% class(data)){
      return(data.frame(colData(data)))
    } else {
      return(data)
    }
  })
  
  output$sampletable <- renderDataTable(
    df_upload(),
    options = list(pageLength = 20,
                   autoWidth = FALSE,
                   lengthMenu = c(20, 30, 40, 50)),
  )
  
  #----Pre-process Data----#
  #User select "log2(exp + 1)" or "VST" or "rlog".
  selectPreprocess <- reactive({
    return(input$preprocess)
  })
  #Pre-process
  preprocess <- reactive({
    dds <- dds_upload()
    preprocess_in <- selectPreprocess()
    #Warning: while file is not unloaded.
    #Whether the input file is a DESeqDataSet.
    if((!is.null(input$upload)) && ("DESeqDataSet" %in% class(dds))){
      #A count at least 10 for a minimal number of sample
      smallestGroupSize <- 4
      keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
      dds <- dds[keep,]
      #Normalize data
      if(preprocess_in == "log2(exp + 1)"){
        prodata <- log2(counts(dds, normalized = TRUE) + 1)
        nordata <- NULL
      } else if(preprocess_in == "VST"){
        nordata <- vst(dds, blind = FALSE)
        prodata <- assay(nordata)
      } else {
        nordata <- rlog(dds, blind = FALSE)
        prodata <- assay(nordata)
      }
      prolist <- list(prodata = prodata,
                      group = dds$treat,
                      pridata = nordata)
      return(prolist)
    }
  })
  #Display this result.
  output$Normalizetable <- renderDataTable(
    preprocess()$prodata,
    options = list(pageLength = 20,
                   autoWidth = FALSE,
                   lengthMenu = c(20, 30, 40, 50)),
  )
  
  #----Display Heatmap Analysis----#
  #User select "Euclidean" or "Poisson".
  selectDistance <- reactive({
    return(input$distancemeasure)
  })
  #Calculate the distance between samples.
  calculateDistance <- reactive({
    prolist_in <- preprocess()
    prodata_in <- prolist_in$prodata
    distance_in <- selectDistance()
    if(distance_in == "Euclidean"){
      sampleDists <- dist(t(prodata_in))
      sampleDistMatrix <- as.matrix(sampleDists)
    } else {
      sampleDists <- PoissonDistance(t(prodata_in))
      sampleDists <- sampleDists$dd
      sampleDistMatrix <- as.matrix(sampleDists)
      rownames(sampleDistMatrix) <- colnames(prodata_in)
      colnames(sampleDistMatrix) <- colnames(prodata_in)
    }
    distancelist <- list(sampleDistMatrix = sampleDistMatrix,
                         sampleDists = sampleDists,
                         group = prolist_in$group)
    return(distancelist)
  })
  # #Show distance matrix
  # output$Distancetable <- renderDataTable(
  #   calculateDistance()[[1]],
  #   options = list(pageLength = 20,
  #                  autoWidth = FALSE,
  #                  lengthMenu = c(20, 30, 40, 50)),
  # )
  #Heatmap the distance matrix.
  output$heatmap <- renderPlot({
    print(paste0("Heatmap Width: ", paste0(input$heatmapheight, "px")))
    print(paste0("Heatmap Height: ", paste0(input$heatmapheight, "px")))
    distancelist_in <- calculateDistance()
    sampleDistMatrix_in <- distancelist_in$sampleDistMatrix
    sampleDists_in <- distancelist_in$sampleDists
    group_in <- distancelist_in$group
    ann_colors <- brewer.pal(11,"PiYG")
    names(ann_colors) <- levels(group_in)
    ann_colors <- list(Group = ann_colors)
    ann_colors$Group <- ann_colors$Group[!is.na(names(ann_colors$Group))]
    ann_col <- data.frame(Group = as.character(group_in))
    ann_col <- HeatmapAnnotation(df = ann_col,
                                 col = ann_colors,
                                 show_legend = FALSE)
    ann_row <- data.frame(Group = as.character(group_in))
    ann_row <- rowAnnotation(df = ann_row,
                             col = ann_colors)
    Heatmap(sampleDistMatrix_in,
            clustering_distance_rows = sampleDists_in,
            clustering_distance_columns = sampleDists_in,
            top_annotation = ann_col,
            right_annotation = ann_row)
  }, width = reactive(as.numeric(input$heatmapwidth)), height = reactive(as.numeric(input$heatmapheight)))
  
  #----Display DR Plot Analysis----#
  #User select "PCA" or "GLM-PCA" or "MDS".
  selectDR <- reactive({
    return(input$drmeasure)
  })
  #Dimensionality reduction for samples
  dimensionalityReduction <- reactive({
    dr_in <- selectDR()
    prolist_in <- preprocess()
    group_in <- prolist_in$group
    proData_in <- prolist_in$prodata
    proData_in <- t(proData_in) %>% as.data.frame(.)
    if(dr_in == "PCA"){
      proPCA <- prcomp(proData_in, scale. = FALSE, center = FALSE)
      pcaData <- proPCA$x %>% as.data.frame(.)
      elemLoading <- (proPCA$sdev/sum(proPCA$sdev)) %>% round(., digits = 1)
      colnames(pcaData) <- paste0("Dim", 1:ncol(pcaData))
      pcaData$group <- group_in
      xLoading <- paste0("PC1:", elemLoading[1], "% variance")
      yLoading <- paste0("PC2:", elemLoading[2], "% variance")
      plotLoading <- list(xLoading = xLoading, yLoading = yLoading)
      drlist <- list(plotData = pcaData, dataloading = plotLoading)
    } else if(dr_in == "GLM-PCA"){
      dds <- dds_upload()
      gpca <- glmpca(counts(dds), L = 2)
      gpca.dat <- gpca$factors
      colnames(gpca.dat) <- paste0("Dim", 1:ncol(gpca.dat))
      gpca.dat$group <- dds$treat
      xLoading <- paste0("Dim1")
      yLoading <- paste0("Dim2")
      plotLoading <- list(xLoading = xLoading, yLoading = yLoading)
      drlist <- list(plotData = gpca.dat, dataloading = plotLoading)
    } else {
      sampleDistMatrix_in <- calculateDistance()$sampleDistMatrix
      mds <- cmdscale(sampleDistMatrix_in) %>% as.data.frame(.)
      colnames(mds) <- paste0("Dim",1:ncol(mds))
      mds$group <- group_in
      xLoading <- paste0("Dim1")
      yLoading <- paste0("Dim2")
      plotLoading <- list(xLoading = xLoading, yLoading = yLoading)
      drlist <- list(plotData = mds, dataloading = plotLoading)
    }
    return(drlist)
  })
  #DR plot the distance matrix
  output$drplot <- renderPlot({
    print(paste0("DR plot Width: ", paste0(input$drwidth, "px")))
    print(paste0("DR plot Height: ", paste0(input$drheight, "px")))
    drlist_in <- dimensionalityReduction()
    drData_in <- drlist_in$plotData
    drLoading_in <- drlist_in$plotLoading
    ggplot(drData_in, aes(x = Dim1, y = Dim2, color = group, shape = group)) +
      geom_point(size = 3) +
      scale_shape_manual(values = 0:11) +
      xlab(drLoading_in$xLoading) +
      ylab(drLoading_in$xLoading) +
      ggtitle("") +
      theme(panel.background = element_blank(),
            panel.grid.minor = element_line(color = "grey"),
            panel.grid.major = element_line(color = "grey"),
            panel.border = element_rect(color = "black", fill = NA))
  }, width = reactive(as.numeric(input$drwidth)), height = reactive(as.numeric(input$drheight)))
  
  #----Download----#
  #Normalized data Table
  output$dataDown <- downloadHandler(
    filename = function(){
      paste0(input$upload, ".", input$preprocess, ".csv")
    },
    content = function(file){
      data <- preprocess()$prodata
      write.csv(data, file = file, quote = FALSE,
                row.names = TRUE, col.names = TRUE)
    }
  )
  #Heatmap
  output$heatmapDown <- downloadHandler(
    filename = function(){
      paste0(input$upload,".Heatmap.",input$format)
    },
    content = function(file){
      if(input$format == "pdf"){
        pdf(file, width = input$heatmapwidth, height = input$heatmapheight)
      } else if(input$format == "png"){
        png(file, width = input$heatmapwidth, height = input$heatmapheight)
      } else {
        jpeg(file, width = input$heatmapwidth, height = input$heatmapheight)
      }
      distancelist_in <- calculateDistance()
      sampleDistMatrix_in <- distancelist_in$sampleDistMatrix
      sampleDists_in <- distancelist_in$sampleDists
      group_in <- distancelist_in$group
      ann_colors <- brewer.pal(11,"PiYG")
      names(ann_colors) <- levels(group_in)
      ann_colors <- list(Group = ann_colors)
      ann_colors$Group <- ann_colors$Group[!is.na(names(ann_colors$Group))]
      ann_col <- data.frame(Group = as.character(group_in))
      ann_col <- HeatmapAnnotation(df = ann_col,
                                   col = ann_colors,
                                   show_legend = FALSE)
      ann_row <- data.frame(Group = as.character(group_in))
      ann_row <- rowAnnotation(df = ann_row,
                               col = ann_colors)
      p <- Heatmap(sampleDistMatrix_in,
              clustering_distance_rows = sampleDists_in,
              clustering_distance_columns = sampleDists_in,
              top_annotation = ann_col,
              right_annotation = ann_row)
      print(p)
      dev.off()
    }
  )
  #PCA plot
  output$PCADown <- downloadHandler(
    filename = function(){
      paste0(input$upload,".PCA.",input$format)
    },
    content = function(file){
      if(input$format == "pdf"){
        pdf(file, width = input$drwidth, height = input$drheight)
      } else if(input$format == "png"){
        png(file, width = input$drwidth, height = input$drheight)
      } else {
        jpeg(file, width = input$drwidth, height = input$drheight)
      }
      drlist_in <- dimensionalityReduction()
      drData_in <- drlist_in$plotData
      drLoading_in <- drlist_in$plotLoading
      p <- ggplot(drData_in, aes(x = Dim1, y = Dim2, color = group, shape = group)) +
        geom_point(size = 3) +
        scale_shape_manual(values = 0:11) +
        xlab(drLoading_in$xLoading) +
        ylab(drLoading_in$xLoading) +
        ggtitle("") +
        theme(panel.background = element_blank(),
              panel.grid.minor = element_line(color = "grey"),
              panel.grid.major = element_line(color = "grey"),
              panel.border = element_rect(color = "black", fill = NA))
      print(p)
      dev.off()
    }
  )
  
  #When a session ends, decrement the counter.
  session$onSessionEnded(function(){
    #We use isolate() here for the same reasons as above.
    isolate(vals$count <- vals$count - 1)
  })

  
  #Reactively update the client.
  output$count <- renderText({vals$count})
})