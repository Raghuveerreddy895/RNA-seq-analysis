#Load libraries#

library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(plotly)
library(biomaRt)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(htmlwidgets)


#1. Read Counts#

counts_raw <- read.csv(
  "<input-path>/Counts.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(counts_raw) <- counts_raw[1, ]
counts_raw <- counts_raw[-1, ]

rownames(counts_raw) <- counts_raw$Geneid

counts_matrix <- counts_raw[, -1]

counts_matrix <- as.matrix(
  sapply(counts_matrix, as.numeric)
)

rownames(counts_matrix) <- rownames(counts_raw)

mode(counts_matrix) <- "integer"

#Remove duplicate genes if present
counts_matrix <- counts_matrix[!duplicated(rownames(counts_matrix)), ]


#2. Read Metadata


metadata <- read.csv(
  "<input-path>/Metadata.csv",
  stringsAsFactors = FALSE
)

metadata <- metadata[
  metadata$ID != "" &
    !is.na(metadata$ID),
]

metadata <- metadata[
  metadata$ID %in% colnames(counts_matrix),
]

counts_matrix <- counts_matrix[, metadata$ID]

metadata$Survival <- factor(
  metadata$Survival,
  levels = c("Long term","Short term") #Update according to the condition#
)

stopifnot(all(colnames(counts_matrix) == metadata$ID))


#3. DESeq2


dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = metadata,
  design = ~ Survival
)

keep <- rowSums(counts(dds) >= 10) >= 5

dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "Survival",
    "Short term",
    "Long term"                   #Update according to the condition#
  )
)


#4. Convert Results


res_df <- as.data.frame(res)

res_df <- tibble::rownames_to_column(
  res_df,
  var = "Geneid"
)

res_df$ensembl_clean <- sub(
  "\\..*",
  "",
  res_df$Geneid
)


#5. Connect to Ensembl


cat("Connecting to Ensembl...\n")

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

annotations <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description"
  ),
  filters = "ensembl_gene_id",
  values = unique(res_df$ensembl_clean),
  mart = ensembl
)


#6. Merge Annotation


res_annotated <-
  dplyr::left_join(
    res_df,
    annotations,
    by = c(
      "ensembl_clean" = "ensembl_gene_id"
    )
  )

res_annotated <-
  dplyr::mutate(
    res_annotated,
    
    gene_symbol = ifelse(
      is.na(external_gene_name) |
        external_gene_name == "",
      Geneid,
      external_gene_name
    ),
    
    Significance =
      dplyr::case_when(
        
        padj < 0.05 &
          log2FoldChange > 1
        ~ "Upregulated ",
        
        padj < 0.05 &
          log2FoldChange < -1
        ~ "Downregulated ",
        
        TRUE
        ~ "Not Significant"
        
      )
  )

res_annotated <-
  dplyr::select(
    res_annotated,
    
    Geneid,
    ensembl_clean,
    gene_symbol,
    log2FoldChange,
    pvalue,
    padj,
    Significance,
    description
  )

res_annotated <-
  dplyr::arrange(
    res_annotated,
    padj
  )

write.csv(
  res_annotated,
  "DEG_Annotated_Results_1.csv",
  row.names = FALSE
)

cat("Annotation completed.\n")


#7. Volcano Plot


plot_data <- res_annotated[
  !is.na(res_annotated$padj),
]

volcano <-
  ggplot(
    plot_data,
    
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      
      color = Significance,
      
      text = paste0(
        
        "Gene: ", gene_symbol,
        
        "<br>Ensembl: ", ensembl_clean,
        
        "<br>log2FC: ",
        round(log2FoldChange,2),
        
        "<br>padj: ",
        signif(padj,3)
        
      )
    )
  ) +
  
  geom_point(alpha = 0.7,size = 2) +
  
  scale_color_manual(
    
    values = c(
      
      "Upregulated "="red",
      
      "Downregulated "="blue",
      
      "Not Significant"="grey70"
      
    )
    
  ) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    
    title="Volcano Plot",
    
    x="Log2 Fold Change",
    
    y="-Log10 Adjusted P-value"
    
  )

interactive_volcano <- ggplotly(
  volcano,
  tooltip="text"
)

htmlwidgets::saveWidget(
  interactive_volcano,
  "Interactive_Volcano.html",
  selfcontained = TRUE
)


#8. Heatmap of Top 50 DEGs


vsd <- vst(dds)

top50 <-
  rownames(
    head(
      res[
        order(res$padj),
      ],
      50
    )
  )

mat <- assay(vsd)[top50, ]

Heatmap(
  mat,
  name="Expression",
  show_row_names=FALSE,
  column_title="Top 50 Differentially Expressed Genes"
)

cat("\nFinished Successfully!\n")

#9. Export Upregulated and Downregulated Genes


#Significant upregulated genes (Short term > Long term)
upregulated_genes <- res_annotated %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.05,
    log2FoldChange > 1
  ) %>%
  dplyr::arrange(padj, dplyr::desc(log2FoldChange))

#Significant downregulated genes (Short term < Long term)
downregulated_genes <- res_annotated %>%
  dplyr::filter(
    !is.na(padj),
    padj < 0.05,
    log2FoldChange < -1
  ) %>%
  dplyr::arrange(padj, log2FoldChange)

#Save CSV files
write.csv(
  upregulated_genes,
  "Upregulated_Genes_ShortTerm_vs_LongTerm.csv",
  row.names = FALSE
)

write.csv(
  downregulated_genes,
  "Downregulated_Genes_ShortTerm_vs_LongTerm.csv",
  row.names = FALSE
)

library(dplyr)

#1. First, define the Status column properly in your main dataset
res_annotated <- res_annotated %>%
  mutate(
    Status = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

#2. Now run your Top 50 Upregulated extraction
top_50_upregulated <- res_annotated %>%
  dplyr::filter(Status == "Upregulated") %>%
  dplyr::arrange(padj) %>%
  head(50)

#3. Pull the Top 50 Downregulated extraction
top_50_downregulated <- res_annotated %>%
  dplyr::filter(Status == "Downregulated") %>%
  dplyr::arrange(padj) %>%
  head(50)

#4. Save your top 50 lists safely to CSV files
write.csv(top_50_upregulated, "Top_50_Upregulated_Genes.csv", row.names = FALSE)
write.csv(top_50_downregulated, "Top_50_Downregulated_Genes.csv", row.names = FALSE)
