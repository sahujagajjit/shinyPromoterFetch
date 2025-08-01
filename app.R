# ProFetShi
# The idea of this tool is to serve as a resource which will provide user-friendly access to promoter sequences of various species (initially plant species)
# Author: JeeT
# Date Created: 29 July, 2025
# Date last modified: 30 July, 2025

# To Do:
## Modify to have more dynamic codes in terms of species
## Add and fit options for accessing species and their datasets on biomart directly from the tool
## Identify packages or write codes to provide promoter analysis options (may be similar to plantcare)
## Check with possibility of parallelization to make some processes a bit faster, for example, downloading all sequences using biomart
## Check if there is any scope for accomodating the species for which neither of BSgenome and biomart is available

# Define CRAN and Bioconductor packages
cranPackages <- c("shiny", "DT", "httr")
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

# Develop dynamic codes for biomart use for those species for which BSgenome is not available
## The first part i.e. Gene can be dynamic for any species for which biomart is available and if no BSgenome then in the next phase while using
### getPromoters4mUsemart function they can be handled and the Dynamic genome part would be saved to the Gene

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

# === Setaria viridis gene ranges ===
ensembl_sviridis <- useMart("plants_mart", dataset = "sviridis_eg_gene", host = "https://plants.ensembl.org")

genesBM_sviridis <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
  mart = ensembl_sviridis
)

sviridisGenes <- GRanges(
  seqnames = genesBM_sviridis$chromosome_name,
  ranges = IRanges(start = genesBM_sviridis$start_position, end = genesBM_sviridis$end_position),
  strand = ifelse(genesBM_sviridis$strand == 1, "+", "-"),
  gene_id = genesBM_sviridis$ensembl_gene_id
)



# === UI ===
ui <- fluidPage(
  titlePanel("shinyPromoterFetch"), # It can be renamed as Promoter Fetching Shiny; in shorts ProFetShi
  sidebarLayout(
    sidebarPanel(
      selectInput("organism", "Choose Organism:", choices = c("Arabidopsis thaliana", "Oryza sativa", "Setaria viridis", "Cicer arietinum")),
      numericInput("promoterLength", "Upstream Length (bp):", value = 1000, min = 100, max = 3000, step = 100),
      numericInput("downstreamLength", "Downstream Length (bp):", value = 0, min = 0, max = 3000, step = 100),
      checkboxInput("namesPromoterSeqs", "SeqIdAddons", value = FALSE),
      br(),
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
    } else if (input$organism == "Setaria viridis") {
      sviridisGenes
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
    } else if (input$organism == "Setaria viridis") {
      sviridisGenes
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
      if(input$namesPromoterSeqs)
      {
        names(promoterSeqs) <- paste0(
          promoterRanges$gene_id, "|",
          seqnames(promoterRanges), ":",
          strand(promoterRanges), "|",
          "Range:", start(promoterRanges), "_", end(promoterRanges),
          "|U:", upstream, "bp|D:", downstream, "bp"
        )
      }else{
        names(promoterSeqs) <- promoterRanges$gene_id
      }
      
    }
    
    return(promoterSeqs)
  }
  
  getPromoters4mUsemart <- function(gr, upstream = 1000, downstream = 0) {
    mart <- useMart("plants_mart", dataset = "sviridis_eg_gene",
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
      
      url <- paste0("https://rest.ensembl.org/sequence/region/Setaria_viridis/",
                    gene$chromosome_name, ":", start, "..", end,
                    ifelse(gene$strand == 1, ":", ":-1"),
                    "?content-type=text/plain")
      
      res <- httr::GET(url)
      if (httr::status_code(res) == 200) {
        seqText <- httr::content(res, "text")
        newSeq <- DNAStringSet(seqText)
        
        if (input$namesPromoterSeqs) {
          names(newSeq) <- paste0(
            gene$ensembl_gene_id, "|", gene$chromosome_name, ":", start, "-", end,
            "|strand:", gene$strand
          )
        } else {
          names(newSeq) <- gene$ensembl_gene_id
        }
        
        promoterSeqs <- append(promoterSeqs, newSeq)
      } else {
        warning(paste("Failed for gene:", gene$ensembl_gene_id))
      }
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
      # Here need to modify to have an option for species for which BSgenome libraries are available if not the sequence can be downloaded using biomart itself
      if (input$organism %in% c("Arabidopsis thaliana", "Oryza sativa")){
        seqs <- getPromoters(
          gr = geneGR,
          upstream = input$promoterLength,
          downstream = input$downstreamLength,
          genome = genome,
          organism = input$organism
        )
      } else if (input$organism == "Setaria viridis"){
        seqs <- getPromoters4mUsemart(
          gr = geneGR,
          upstream = input$promoterLength,
          downstream = input$downstreamLength
          # organism = input$organism
        )
      }
      
      
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
