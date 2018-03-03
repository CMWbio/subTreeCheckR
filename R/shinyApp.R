library(shiny)
library(shinydashboard)
library(leaflet)
library(shinyFiles)




popShiny <- function(fileName, contigs = "all"){

  header <- scanVcfHeader(TabixFile(fileName))
  contigMD <- as.data.frame(header@header$contig)
  if(contigs == "all") contigs <- rownames(contigMD)

  ui <- dashboardPage(
    dashboardHeader(title = "Population Statistics"),
    dashboardSidebar(
      sidebarMenu(
        menuItem(text = "", icon = icon("compressed", lib = "glyphicon"), tabName = "readVCF"),
        menuItem(text = "", icon = icon("tree"), tabName = "trees")
      ), width = "5%"
    ),
    dashboardBody(
      tabItem(tabName = "readVCF",
              h2("Read in VCF Data"),
              h5("Select VCF File"),
              # shinyFilesButton(id = "VCF", title = "Select VCF File", label = "Browse...", multiple = FALSE),
              textInput(inputId = "winSize", value = "10000", label = "Select Window Size")),
      selectizeInput("scaffoldNames", choices = NULL, label = "Select Contigs to Read in from VCF", multiple = TRUE),
      textInput("minSites", value = "1000", label = "Minimum sites in window to calculate from"),
      textInput("ploidy", value = "2", label = "Ploidy of samples"),
      actionButton("import", "Import Windows"),
      textOutput("nrows")

    )
  )

  server <- function(input, output, session){

    # volumes <- getVolumes()
    # observe({
    #   shinyFileChoose(input, "VCF", roots = volumes, session = session)
    # })

    values <- reactiveValues()

    updateSelectizeInput(session, "scaffoldNames", choices = c("all", contigs))

    observeEvent(input$import, {
      winSize <- as.numeric(input$winSize)
      percentage <- 0
      ploidy <- input$ploidy
      minSites <- input$minSites

      if(all(input$scaffoldNames == "all")){
        scaf <- contigs
      } else {
        scaf <- input$scaffoldNames
      }
      withProgress(
        values$dna <- lapply(scaf,  function(con){
          percentage <<- percentage + 1/length(scaf)*100
          incProgress(1/length(scaf), detail = paste0("Progress: ", round(percentage,2)))
          length <- as.integer(filter(contigMD, rownames(contigMD) == con)$length)
          if(length >= winSize){
            nWindows <- floor(length / winSize)

            scafDNA <- mclapply(seq(1, nWindows), mc.cores = nCores, function(winN){

              pos <- winN * winSize + 1
              start <- pos - winSize
              end <- pos

              p <- ScanVcfParam(which = GRanges(seqnames = scaf, ranges = IRanges(start = start, end = end)))

              nSites <- tryCatch(length(scanVcf(TabixFile(fileName), param = p)[[1]]$rowRanges),  error=function(e) 0)

              if(nSites >= minSites){
                #read in vcf
                dna <- vcfWindow(fileName = fileName, contig = scaf, param = p, ploidy = ploidy)
              } else dna <- c(NA)
              names(dna) <- paste0(scaf, ":", start, "-", end)
              dna
            })
            scafDNA
          } else scafDNA <- NA
        }) %>% unlist(recursive = FALSE),
        message = "Reading in Windows"

      )
    })



    output$nrows <- renderText({
      if(!is.null(values$dna)) length(values$dna)
    })


  }


  shinyApp(ui = ui, server = server)

}

