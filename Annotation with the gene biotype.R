library(biomaRt)
library(stringr)

# Read raw counts
counts <- read.csv("C:/Users/raghu/Downloads/Meta-analysis MSMF/RNA_Primary/Counts_primary.csv",
                   row.names = 1,
                   check.names = FALSE)

# Remove Ensembl version numbers
ensembl_ids <- str_remove(rownames(counts), "\\..*$")

# Connect to Ensembl
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# Get annotations
annotation <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "description",
    "entrezgene_id",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  filters = "ensembl_gene_id",
  values = unique(ensembl_ids),
  mart = mart
)

# Keep original order
annotation <- annotation[
  match(ensembl_ids, annotation$ensembl_gene_id),
]

# Combine with counts
annotated_counts <- cbind(annotation, counts)

# Save
write.csv(
  annotated_counts,
  "Annotated_Raw_Counts_1.csv",
  row.names = FALSE
)

