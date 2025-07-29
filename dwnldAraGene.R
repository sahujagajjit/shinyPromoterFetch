# dwnldAraGene.R
# Purpose: Download Arabidopsis gene metadata and full sequences, and save to an RDS file in the script's directory

# Load required libraries
suppressMessages({
  library(biomaRt)
  library(BSgenome.Athaliana.TAIR.TAIR9)
  library(GenomicRanges)
  library(Biostrings)
  library(this.path)  # helps get the script's own directory
})

# Get the directory of the script itself
scriptDir <- this.path::this.dir()

# Define output path in the same directory
outputPath <- file.path(scriptDir, "athalianaGenes.rds")

# Connect to Ensembl Plants via biomaRt
ensembl <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")

# Retrieve gene metadata
geneTable <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                  "chromosome_name", "start_position", "end_position", "strand"),
                   mart = ensembl)

# Filter only valid chromosomes (TAIR9: Chr1 to Chr5)
geneTable <- subset(geneTable, chromosome_name %in% c("1", "2", "3", "4", "5"))

# Create GRanges object for genes
geneGR <- GRanges(seqnames = paste0("Chr", geneTable$chromosome_name),
                  ranges = IRanges(start = geneTable$start_position, end = geneTable$end_position),
                  strand = ifelse(geneTable$strand == 1, "+", "-"),
                  gene_id = geneTable$ensembl_gene_id,
                  gene_name = geneTable$external_gene_name)

# Fetch gene sequences from BSgenome
geneSeqs <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, geneGR)

# Build final data frame
geneDF <- data.frame(
  GeneID = mcols(geneGR)$gene_id,
  GeneName = mcols(geneGR)$gene_name,
  Chromosome = as.character(seqnames(geneGR)),
  Start = start(geneGR),
  End = end(geneGR),
  Strand = as.character(strand(geneGR)),
  Length = width(geneGR),
  stringsAsFactors = FALSE
)

# Add sequence as character column
geneDF$Sequence <- as.character(geneSeqs)

# Save to .rds file
saveRDS(geneDF, file = outputPath)

# Confirmation
message("Arabidopsis gene data saved to: ", outputPath)
