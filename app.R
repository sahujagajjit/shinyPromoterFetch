library(shiny)
library(DT)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(biomaRt)

# === Generate Arabidopsis gene ranges directly using biomaRt ===
ensembl <- useMart("plants_mart", 
                   dataset = "athaliana_eg_gene", 
                   host = "https://plants.ensembl.org")

genesBM <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", 
                 "start_position", "end_position", "strand"),
  mart = ensembl
)

athalianaGenes <- GRanges(
  seqnames = genesBM$chromosome_name,
  ranges = IRanges(start = genesBM$start_position,
                   end = genesBM$end_position),
  strand = ifelse(genesBM$strand == 1, "+", "-"),
  gene_id = genesBM$ensembl_gene_id
)

athalianaGenome <- BSgenome.Athaliana.TAIR.TAIR9

# === Define UI ===
ui <- fluidPage(
  titlePanel("shinyPromoterFetch"),
  sidebarLayout(
    sidebarPanel(
      selectInput("organism", "Choose Organism:", choices = c("Arabidopsis thaliana", "Cicer arietinum")),
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
  
  output$statusText <- renderText({
    if (input$organism == "Arabidopsis thaliana") {
      paste("Total number of genes:", length(athalianaGenes))
    } else {
      "Data for Cicer arietinum is not available yet."
    }
  })
  
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
  
  getArabidopsisPromoters <- function(gr, upstream = 1000, downstream = 0,
                                      genome = BSgenome.Athaliana.TAIR.TAIR9,
                                      namesPromoterSeqs = TRUE) {
    seqlevels(gr) <- paste0("Chr", seqlevels(gr))
    promoterRanges <- promoters(gr, upstream = upstream, downstream = downstream)
    seqlevels(promoterRanges) <- sub("ChrPt", "ChrC", seqlevels(promoterRanges))
    seqlevels(promoterRanges) <- sub("ChrMt", "ChrM", seqlevels(promoterRanges))
    commonSeqs <- intersect(seqlevels(promoterRanges), seqlevels(genome))
    promoterRanges <- keepSeqlevels(promoterRanges, commonSeqs, pruning.mode = "coarse")
    seqlengths(promoterRanges) <- seqlengths(genome)[seqlevels(promoterRanges)]
    promoterRanges <- trim(promoterRanges)
    promoterSeqs <- getSeq(genome, promoterRanges)
    
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
    
    return(promoterSeqs)
  }
  
  observeEvent(input$submit, {
    output$downloadUI <- renderUI(NULL)
    
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
    
    output$downloadUI <- renderUI({
      req(sequenceData())
      downloadButton("downloadBtn", "Download Promoter FASTA")
    })
  })
  
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
