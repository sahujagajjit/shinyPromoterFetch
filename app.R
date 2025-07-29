library(shiny)
library(DT)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(TxDb.Athaliana.BioMart.plantsmart28)

# === Load Arabidopsis gene ranges from RDS once ===
# rdsPath <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "athalianaGeneRanges.rds")
# if (!file.exists(rdsPath)) {
#   txdb <- TxDb.Athaliana.BioMart.plantsmart28
#   genes <- genes(txdb)
#   genes$gene_id <- names(genes)
#   saveRDS(genes, rdsPath)
# }
# athalianaGenes <- readRDS(rdsPath)
if (file.exists("athalianaGenes.rds")) {
  athalianaGenes <- readRDS("athalianaGenes.rds")
} else if (interactive()) {
  message("athalianaGenes.rds not found. Generating from TxDb...")
  txdb <- TxDb.Athaliana.BioMart.plantsmart28
  genes <- genes(txdb)
  genes$gene_id <- names(genes)
  saveRDS(genes, "athalianaGenes.rds")
  athalianaGenes <- genes
} else {
  stop("athalianaGenes.rds not found and cannot generate in non-interactive mode.")
}

athalianaGenome <- BSgenome.Athaliana.TAIR.TAIR9

# === Define UI ===
ui <- fluidPage(
  titlePanel("shinyPromoterFetch"),
  sidebarLayout(
    sidebarPanel(
      selectInput("organism", "Choose Organism:", choices = c("Arabidopsis thaliana", "Cicer arietinum")),
      numericInput("promoterLength", "Upstream Length (bp):", value = 1000, min = 100, max = 3000, step = 100),
      numericInput("downstreamLength", "Downstream Length (bp):", value = 0, min = 0, max = 3000, step = 100),  # 👈 New field
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
  
  # Text: Number of genes or unavailability message
  output$statusText <- renderText({
    if (input$organism == "Arabidopsis thaliana") {
      paste("Total number of genes:", length(athalianaGenes))
    } else {
      "Data for Cicer arietinum is not available yet."
    }
  })
  
  # Show gene table
  observeEvent(input$selectIds, {
    req(input$organism == "Arabidopsis thaliana")
    
    geneDF <- data.frame(
      GeneID = athalianaGenes$gene_id,
      GeneName = NA,
      Chromosome = as.character(seqnames(athalianaGenes)),
      Start = start(athalianaGenes),
      End = end(athalianaGenes),
      Strand = as.character(strand(athalianaGenes)),
      Length = width(athalianaGenes)
    )
    
    selectedData(geneDF)
    processed(TRUE)
  })
  
  output$geneTable <- renderDT({
    req(processed())
    datatable(selectedData(), selection = 'multiple', options = list(scrollX = TRUE))
  })
  
  # Function to extract promoter sequences
  getArabidopsisPromoters <- function(gr,
                                      upstream = 1000,
                                      downstream = 0,
                                      genome = BSgenome.Athaliana.TAIR.TAIR9,
                                      namesPromoterSeqs = TRUE,
                                      outputFasta = NULL) {
    # Step 1: Ensure chromosomes match between GR and genome
    seqlevels(gr) <- paste0("Chr", seqlevels(gr))  # Prefix "Chr" if not already
    
    # Step 2: Define TSS based on strand and build promoter ranges
    promoterRanges <- promoters(gr, upstream = upstream, downstream = downstream)
    
    # Step 3: Adjust unusual chromosome naming if needed
    seqlevels(promoterRanges) <- sub("ChrPt", "ChrC", seqlevels(promoterRanges))
    seqlevels(promoterRanges) <- sub("ChrMt", "ChrM", seqlevels(promoterRanges))
    
    # Step 4: Keep only valid chromosomes
    commonSeqs <- intersect(seqlevels(promoterRanges), seqlevels(genome))
    promoterRanges <- keepSeqlevels(promoterRanges, commonSeqs, pruning.mode = "coarse")
    
    # Step 5: Trim to genome bounds
    seqlengths(promoterRanges) <- seqlengths(genome)[seqlevels(promoterRanges)]
    promoterRanges <- trim(promoterRanges)
    
    # Step 6: Extract promoter sequences
    promoterSeqs <- getSeq(genome, promoterRanges)
    
    # Step 7: Add informative names
    if (namesPromoterSeqs) {
      names(promoterSeqs) <- paste0(
        promoterRanges$gene_id, "|",
        seqnames(promoterRanges), ":",
        strand(promoterRanges), "|",
        "Range:", start(promoterRanges), "_", end(promoterRanges),
        "|U:", input$promoterLength, "bp|D:", input$downstreamLength, "bp"
      )
    } else {
      names(promoterSeqs) <- promoterRanges$gene_id
    }
    
    # Step 8: Optionally write to FASTA
    if (!is.null(outputFasta)) {
      writeXStringSet(promoterSeqs, filepath = outputFasta)
      message("Promoter sequences written to ", outputFasta)
    }
    
    return(promoterSeqs)
  }
  
  # When Submit is clicked
  observeEvent(input$submit, {
    output$downloadUI <- renderUI(NULL)  # Disable download button initially
    
    if (input$organism == "Cicer arietinum") {
      sequenceData(NULL)
      showNotification("Promoter data for Cicer arietinum is not available.", type = "warning")
      return()
    }
    
    withProgress(message = "Extracting promoters...", value = 0, {
      incProgress(0.2, detail = "Preparing gene ranges...")
      geneGR <- athalianaGenes
      
      if (processed()) {
        sel <- input$geneTable_rows_selected
        if (!is.null(sel) && length(sel) > 0) {
          ids <- selectedData()$GeneID[sel]
          geneGR <- athalianaGenes[athalianaGenes$gene_id %in% ids]
        }
      }
      
      incProgress(0.5, detail = "Fetching sequences...")
      seqs <- getArabidopsisPromoters(
        gr = geneGR,
        upstream = input$promoterLength,
        downstream = input$downstreamLength,
        genome = athalianaGenome
      )
      
      incProgress(0.3, detail = "Done.")
      sequenceData(seqs)
    })
    
    # ✅ Show download button when ready
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

# === Launch ===
shinyApp(ui, server)
