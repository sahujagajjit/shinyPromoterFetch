# Define CRAN and Bioconductor packages
cranPackages <- c("shiny", "DT")
biocPackages <- c(
  "Biostrings", "GenomicRanges", "BSgenome", 
  "BSgenome.Athaliana.TAIR.TAIR9", "BSgenome.Osativa.MSU.MSU7", 
  "TxDb.Athaliana.BioMart.plantsmart28", "biomaRt"
)

# Install CRAN packages if not present
installIfMissing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages if not present
installBiocIfMissing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# Apply installation
invisible(lapply(cranPackages, installIfMissing))
invisible(lapply(biocPackages, installBiocIfMissing))

# Load all packages
library(shiny)
library(DT)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(BSgenome.Osativa.MSU.MSU7)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(biomaRt)


# === Arabidopsis gene ranges ===
ensembl_ath <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")
genesBM_ath <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
  mart = ensembl_ath
)
athalianaGenes <- GRanges(
  seqnames = genesBM_ath$chromosome_name,
  ranges = IRanges(start = genesBM_ath$start_position, end = genesBM_ath$end_position),
  strand = ifelse(genesBM_ath$strand == 1, "+", "-"),
  gene_id = genesBM_ath$ensembl_gene_id
)
athalianaGenome <- BSgenome.Athaliana.TAIR.TAIR9

# === Oryza sativa gene ranges ===
ensembl_os <- useMart("plants_mart", dataset = "osativa_eg_gene", host = "https://plants.ensembl.org")
genesBM_os <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
  mart = ensembl_os
)
osativaGenes <- GRanges(
  seqnames = genesBM_os$chromosome_name,
  ranges = IRanges(start = genesBM_os$start_position, end = genesBM_os$end_position),
  strand = ifelse(genesBM_os$strand == 1, "+", "-"),
  gene_id = genesBM_os$ensembl_gene_id
)
osativaGenome <- BSgenome.Osativa.MSU.MSU7

# === UI ===
ui <- fluidPage(
  titlePanel("shinyPromoterFetch"),
  sidebarLayout(
    sidebarPanel(
      selectInput("organism", "Choose Organism:", choices = c("Arabidopsis thaliana", "Oryza sativa", "Cicer arietinum")),
      numericInput("promoterLength", "Upstream Length (bp):", value = 1000, min = 100, max = 3000, step = 100),
      numericInput("downstreamLength", "Downstream Length (bp):", value = 0, min = 0, max = 3000, step = 100),
      actionButton("submit", "Submit"),
      uiOutput("geneSelector")
    ),
    mainPanel(
      textOutput("statusText"),
      actionButton("selectIds", "Select from IDs"),
      DTOutput("geneTable"),
      br(),
      uiOutput("downloadUI") 
    )
  )
)

# === Server ===
server <- function(input, output, session) {
  selectedData <- reactiveVal(NULL)
  processed <- reactiveVal(FALSE)
  sequenceData <- reactiveVal(NULL)
  
  # Dynamic gene ranges
  getSelectedGenes <- reactive({
    if (input$organism == "Arabidopsis thaliana") {
      athalianaGenes
    } else if (input$organism == "Oryza sativa") {
      osativaGenes
    } else {
      NULL
    }
  })
  
  # Dynamic genome
  getSelectedGenome <- reactive({
    if (input$organism == "Arabidopsis thaliana") {
      athalianaGenome
    } else if (input$organism == "Oryza sativa") {
      osativaGenome
    } else {
      NULL
    }
  })
  
  # Status text
  output$statusText <- renderText({
    genes <- getSelectedGenes()
    if (is.null(genes)) {
      "Data for Cicer arietinum is not available yet."
    } else {
      paste("Total number of genes:", length(genes))
    }
  })
  
  # Gene table
  observeEvent(input$selectIds, {
    genes <- getSelectedGenes()
    req(!is.null(genes))
    geneDF <- data.frame(
      GeneID = genes$gene_id,
      GeneName = NA,
      Chromosome = as.character(seqnames(genes)),
      Start = start(genes),
      End = end(genes),
      Strand = as.character(strand(genes)),
      Length = width(genes)
    )
    selectedData(geneDF)
    processed(TRUE)
  })
  
  output$geneTable <- renderDT({
    req(processed())
    datatable(selectedData(), selection = 'multiple', options = list(scrollX = TRUE))
  })
  
  # General promoter function
  getPromoters <- function(gr, upstream, downstream, genome, organism) {
    if (organism == "Arabidopsis thaliana") {
      seqlevels(gr) <- paste0("Chr", seqlevels(gr))
      seqlevels(gr) <- sub("ChrPt", "ChrC", seqlevels(gr))
      seqlevels(gr) <- sub("ChrMt", "ChrM", seqlevels(gr))
    } else if (organism == "Oryza sativa") {
      # Adjust chromosome names for Oryza to match BSgenome
      seqlevels(gr) <- paste0("Chr", seqlevels(gr))
    }
    
    promoterRanges <- promoters(gr, upstream = upstream, downstream = downstream)
    commonSeqs <- intersect(seqlevels(promoterRanges), seqlevels(genome))
    promoterRanges <- keepSeqlevels(promoterRanges, commonSeqs, pruning.mode = "coarse")
    
    if (length(promoterRanges) == 0) return(DNAStringSet())
    
    seqlengths(promoterRanges) <- seqlengths(genome)[seqlevels(promoterRanges)]
    promoterRanges <- trim(promoterRanges)
    promoterSeqs <- getSeq(genome, promoterRanges)
    
    if (length(promoterSeqs) > 0) {
      names(promoterSeqs) <- paste0(
        promoterRanges$gene_id, "|",
        seqnames(promoterRanges), ":",
        strand(promoterRanges), "|",
        "Range:", start(promoterRanges), "_", end(promoterRanges),
        "|U:", upstream, "bp|D:", downstream, "bp"
      )
    }
    
    return(promoterSeqs)
  }
  
  # Handle promoter extraction
  observeEvent(input$submit, {
    output$downloadUI <- renderUI(NULL)
    
    if (input$organism == "Cicer arietinum") {
      sequenceData(NULL)
      showNotification("Promoter data for Cicer arietinum is not available.", type = "warning")
      return()
    }
    
    withProgress(message = "Extracting promoters...", value = 0, {
      incProgress(0.2, detail = "Preparing gene ranges...")
      genes <- getSelectedGenes()
      genome <- getSelectedGenome()
      
      req(!is.null(genes) && !is.null(genome))
      geneGR <- genes
      
      if (processed()) {
        sel <- input$geneTable_rows_selected
        if (!is.null(sel) && length(sel) > 0) {
          ids <- selectedData()$GeneID[sel]
          geneGR <- genes[genes$gene_id %in% ids]
        }
      }
      
      incProgress(0.5, detail = "Fetching sequences...")
      seqs <- getPromoters(
        gr = geneGR,
        upstream = input$promoterLength,
        downstream = input$downstreamLength,
        genome = genome,
        organism = input$organism
      )
      
      incProgress(0.3, detail = "Done.")
      sequenceData(seqs)
    })
    
    output$downloadUI <- renderUI({
      req(sequenceData())
      downloadButton("downloadBtn", "Download Promoter FASTA")
    })
  })
  
  # Download handler
  output$downloadBtn <- downloadHandler(
    filename = function() {
      paste0("promoters_", gsub(" ", "_", input$organism), ".fasta")
    },
    content = function(file) {
      req(sequenceData())
      writeXStringSet(sequenceData(), filepath = file)
    }
  )
}

# === Launch App ===
shinyApp(ui, server)
