library(shiny)
library(DT)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(TxDb.Athaliana.BioMart.plantsmart28)

# === Load Arabidopsis gene ranges from RDS once ===
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    !is.null(rstudioapi::getActiveDocumentContext()$path)) {
  rdsPath <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "athalianaGeneRanges.rds")
} else {
  rdsPath <- "athalianaGeneRanges.rds"
}

if (!file.exists(rdsPath)) {
  txdb <- TxDb.Athaliana.BioMart.plantsmart28
  genes <- genes(txdb)
  genes$gene_id <- names(genes)
  saveRDS(genes, rdsPath)
}
athalianaGenes <- readRDS(rdsPath)
athalianaGenome <- BSgenome.Athaliana.TAIR.TAIR9

# === Define UI ===
ui <- fluidPage(
  titlePanel("shinyPromoterFetch"),
  sidebarLayout(
    sidebarPanel(
      selectInput("organism", "Choose Organism:", choices = c("Arabidopsis thaliana", "Cicer arietinum")),
      numericInput("promoterLength", "Upstream Length:", value = 1000, min = 100, max = 3000, step = 100),
      actionButton("submit", "Submit"),
      br(),
      downloadButton("downloadBtn", "Download Promoter FASTA")
    ),
    mainPanel(
      textOutput("statusText"),
      actionButton("selectIds", "Select from IDs"),
      DTOutput("geneTable")
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
                                      genome = BSgenome.Athaliana.TAIR.TAIR9,
                                      namesPromoterSeqs = TRUE,
                                      outputFasta = NULL) {
    # Step 1: Ensure chromosomes match between GR and genome
    seqlevels(gr) <- paste0("Chr", seqlevels(gr))  # Prefix "Chr" if not already
    
    # Step 2: Define TSS based on strand and build promoter ranges
    promoterRanges <- promoters(gr, upstream = upstream, downstream = 0)
    
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
        start(promoterRanges), "_", end(promoterRanges)
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
    if (input$organism == "Cicer arietinum") {
      sequenceData(NULL)
      showNotification("Promoter data for Cicer arietinum is not available.", type = "warning")
      return()
    }
    
    # Start with all genes
    geneGR <- athalianaGenes
    
    # If user selected specific genes
    if (processed()) {
      sel <- input$geneTable_rows_selected
      if (!is.null(sel) && length(sel) > 0) {
        ids <- selectedData()$GeneID[sel]
        geneGR <- athalianaGenes[athalianaGenes$gene_id %in% ids]
      }
    }
    
    # Extract promoters
    seqs <- getArabidopsisPromoters(
      gr = geneGR,
      upstream = input$promoterLength,
      genome = athalianaGenome
    )
    
    sequenceData(seqs)
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
