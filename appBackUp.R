library(shiny)
library(DT)
library(Biostrings)
library(GenomicRanges)
# library(BSgenome)
# library(BSgenome.Athaliana.TAIR.TAIR9)
# library(TxDb.Athaliana.BioMart.plantsmart28)
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
  
  library(httr)
  library(Biostrings)
  
  getArabidopsisPromoters <- function(gr, upstream = 1000, downstream = 0,
                                      namesPromoterSeqs = TRUE) {
    # Get coordinates using biomaRt
    mart <- useMart("plants_mart", dataset = "athaliana_eg_gene",
                    host = "https://plants.ensembl.org")
    
    coords <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "strand", 
                                   "start_position", "end_position"),
                    filters = "ensembl_gene_id",
                    values = gr$gene_id,
                    mart = mart)
    
    promoterSeqs <- DNAStringSet()
    
    for (i in seq_len(nrow(coords))) {
      gene <- coords[i, ]
      if (gene$strand == 1) {
        start <- max(1, gene$start_position - upstream)
        end <- gene$start_position + downstream - 1
      } else {
        start <- gene$end_position - downstream + 1
        end <- gene$end_position + upstream
      }
      
      url <- paste0("https://rest.ensembl.org/sequence/region/arabidopsis_thaliana/",
                    gene$chromosome_name, ":", start, "..", end,
                    ifelse(gene$strand == 1, ":", ":-1"),
                    "?content-type=text/plain")
      
      res <- httr::GET(url)
      if (status_code(res) == 200) {
        seqText <- httr::content(res, "text")
        promoterSeqs <- append(promoterSeqs, DNAStringSet(seqText))
        if (namesPromoterSeqs) {
          names(promoterSeqs)[length(promoterSeqs)] <- paste0(
            gene$ensembl_gene_id, "|", gene$chromosome_name, ":", start, "-", end,
            "|strand:", gene$strand
          )
        }
      } else {
        warning(paste("Failed for gene:", gene$ensembl_gene_id))
      }
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
      incProgress(0.2, detail = "Preparing gene list...")
      geneGR <- athalianaGenes
      
      if (processed()) {
        sel <- input$geneTable_rows_selected
        if (!is.null(sel) && length(sel) > 0) {
          ids <- selectedData()$GeneID[sel]
          geneGR <- athalianaGenes[athalianaGenes$gene_id %in% ids]
        }
      }
      
      incProgress(0.5, detail = "Fetching sequences using biomaRt...")
      seqs <- getArabidopsisPromoters(
        gr = geneGR,
        upstream = input$promoterLength,
        downstream = input$downstreamLength
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
