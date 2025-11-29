library(ChIPseeker)
library(DESeq2)
library(stringr)
library(limma)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
library(Gviz)
library(AnnotationHub)


get_top_genes_for_mark <- function(histoneMark,
                                   save_dir = "./") {

  data_dir     <- str_c(save_dir, "data/promoter/", histoneMark, "/")
  analysis_dir <- str_c(save_dir, "analysis/promoter/", histoneMark, "/")

  # List all up/down CSV files
  files <- list.files(data_dir, pattern = "_genes.csv$", full.names = TRUE)

  # Extract contrasts (GAL1_Ctrl, G3BP_Ctrl, POSTN_Ctrl)
  contrasts <- unique(gsub("_UP_genes.csv|_DOWN_genes.csv", "", basename(files)))

  result_list <- list()

  for (cont in contrasts) {

    up_file   <- str_c(data_dir, cont, "_UP_genes.csv")
    down_file <- str_c(data_dir, cont, "_DOWN_genes.csv")

    most_up <- NA
    most_down <- NA

    # UP genes: choose highest logFC
    if (file.exists(up_file)) {
      df_up <- read.csv(up_file)
      if (nrow(df_up) > 0) {
        # highest logFC
        best <- df_up[which.max(df_up$logFC), ]
        most_up <- best$SYMBOL
      }
    }

    # DOWN genes: choose lowest logFC
    if (file.exists(down_file)) {
      df_down <- read.csv(down_file)
      if (nrow(df_down) > 0) {
        # lowest logFC
        best <- df_down[which.min(df_down$logFC), ]
        most_down <- best$SYMBOL
      }
    }

    result_list[[cont]] <- data.frame(
      contrast     = cont,
      top_up_gene  = most_up,
      top_down_gene = most_down,
      stringsAsFactors = FALSE
    )
  }

  # Combine results
  result_df <- do.call(rbind, result_list)

  # Save summary
  out_file <- str_c(analysis_dir, "top_genes_summary_", histoneMark, ".csv")
  write.csv(result_df, out_file, row.names = FALSE)

  message("Saved: ", out_file)

  return(result_df)
}

top_A1 <- get_top_genes_for_mark("A1")
top_A2 <- get_top_genes_for_mark("A2")
top_P  <- get_top_genes_for_mark("P")
top_R1 <- get_top_genes_for_mark("R1")
top_R2 <- get_top_genes_for_mark("R2")

###############################################################
# Function: plot_gene_gviz
#
# histoneMark:      "A1", "A2", "C", "P", "R1", "R2"
# gene:             gene symbol
# sample_type:      "G3BP", "GAL1", "POSTN", etc.
# replacement_file: CSV file containing columns: x, replacement
# geneModels_file:  CSV with gene annotation
# bw_dir:           directory containing bigwig files
# save_dir:         where to store PDF output
###############################################################

geneModels = read.csv("./ref/geneModels.csv",header = T,row.names = 1)

plot_gene_gviz <- function(histoneMark,
                           gene,
                           sample_type,
                           replacement_file,
                           geneModels,
                           bw_dir = "./bw/",
                           save_dir = "./analysis/promoter/",
                           promoter_extend = 3000,
                           suffix = "") {
  # Skip if NA gene
  if (is.na(gene) || is.null(gene) || gene == "") {
    message("Skipping plot: gene_to_plot is NA for ", 
            sample_type, " in ", histoneMark)
    return(NULL)
  }
  message("====== Plotting Gviz for ", gene, " in ", histoneMark, " =====")

  # ---------------------------------------------------------
  # 1. Load mapping file (your replacement_file)
  # ---------------------------------------------------------
  map <- read.csv(replacement_file, header = TRUE)

  # The screenshot shows:
  # map$x            (original filename, no extension)
  # map$replacement  (standardized name)
  # Filter: only samples of this histone mark
  map_sub <- map[grepl(str_c("_", histoneMark, "_"), map$replacement), ]

  # Control = begins with "Ctrl"
  ctrl_samples <- map_sub[grepl("^Ctrl", map_sub$replacement), ]
  print(ctrl_samples)

  # Treatment = begins with chosen sample_type
  treat_samples <- map_sub[grepl(paste0("^", sample_type), map_sub$replacement), ]
  print(treat_samples)

  # ---------------------------------------------------------
  # 2. Load BigWigs
  # bigwig path = x + ".bw" (from your screenshot)
  # ---------------------------------------------------------
  load_bw <- function(x) {
    bw_file <- str_c(bw_dir, x, ".bw")
    if (!file.exists(bw_file))
      stop("File not found: ", bw_file)
    import.bw(bw_file)
  }

  ctrl_bw  <- lapply(ctrl_samples$x, load_bw)
  treat_bw <- lapply(treat_samples$x, load_bw)


  # ---------------------------------------------------------
  # 3. Gene Model
  # ---------------------------------------------------------
  gm <- geneModels[geneModels$symbol == gene, ]
  
  if (nrow(gm) == 0) {
      message("Skipping plot for ", gene,
              ": gene not found in geneModels_file.")
      return(NULL)
  }
  gm <- gm[!gm$exon %in% 1, ]

  chr     <- as.character(unique(gm$chromosome))
  start_g <- min(gm$start)
  end_g   <- max(gm$end)

  # Query region (±3000 bp)
  q <- GRanges(seqnames = chr,
               ranges = IRanges(
                 start = start_g - promoter_extend,
                 end   = end_g   + promoter_extend),
               strand = "*")


  # ---------------------------------------------------------
  # 4. Subset bigwig regions to gene region
  # ---------------------------------------------------------
  subset_bw <- function(bw_list) {
    lapply(bw_list, function(x) subsetByOverlaps(x, q))
  }

  ctrl_bw_sub  <- subset_bw(ctrl_bw)
  treat_bw_sub <- subset_bw(treat_bw)

  # Determine max y-scale
  all_scores <- unlist(lapply(c(ctrl_bw_sub, treat_bw_sub), function(x) x$score))
  ymax <- max(all_scores, na.rm = TRUE)


  # ---------------------------------------------------------
  # 5. Build Gviz Tracks
  # ---------------------------------------------------------
  itrack  <- IdeogramTrack(genome = "hg38", chromosome = chr)
  gtrack  <- GenomeAxisTrack()
  grtrack <- GeneRegionTrack(gm, genome = "hg38", chromosome = chr,
                             name = str_c("Gene: ", gene))

  make_dt <- function(gr, name, color) {
    DataTrack(range = gr, name = name, genome = "hg38",
              type = "polygon", col = color, ylim = c(0, ymax))
  }

  dtracks <- list()

  # Ctrl tracks
  for (i in seq_along(ctrl_bw_sub)) {
    nm <- ctrl_samples$replacement[i]
    dtracks[[length(dtracks)+1]] <- make_dt(ctrl_bw_sub[[i]], nm, "red")
  }

  # Treatment tracks
  for (i in seq_along(treat_bw_sub)) {
    nm <- treat_samples$replacement[i]
    dtracks[[length(dtracks)+1]] <- make_dt(treat_bw_sub[[i]], nm, "purple")
  }


  # ---------------------------------------------------------
  # 6. Save PDF
  # ---------------------------------------------------------
  out_pdf <- str_c(save_dir, histoneMark, "/Gviz_", gene, "_", sample_type, "_", histoneMark, "_", suffix, ".pdf")
  dir.create(str_c(save_dir, histoneMark, "/"), recursive = TRUE, showWarnings = FALSE)

  pdf(out_pdf, width = 12, height = 15)
  plotTracks(
    c(list(itrack, gtrack, grtrack), dtracks),
    transcriptAnnotation = "transcript"
  )
  dev.off()

  message("Saved: ", out_pdf)
  return(out_pdf)
}

for (i in c("A1", "A2", "P", "R1", "R2")) {
  top_variable <- get(str_c("top_", i))
  gene_to_plot <- top_variable[top_variable$contrast == "G3BP_Ctrl", "top_up_gene"]
  plot_gene_gviz(
    histoneMark = i,
    gene = gene_to_plot,
    sample_type = "G3BP",
    replacement_file = "./data/sample_name_replacement.csv",
    geneModels = geneModels,
    bw_dir = "./data/bw/",
    suffix = "up"
  )

  gene_to_plot <- top_variable[top_variable$contrast == "GAL1_Ctrl", "top_up_gene"]
  plot_gene_gviz(
    histoneMark = i,
    gene = gene_to_plot,
    sample_type = "GAL1",
    replacement_file = "./data/sample_name_replacement.csv",
    geneModels = geneModels,
    bw_dir = "./data/bw/",
    suffix = "up"
  )

  gene_to_plot <- top_variable[top_variable$contrast == "POSTN_Ctrl", "top_up_gene"]
  plot_gene_gviz(
    histoneMark = i,
    gene = gene_to_plot,
    sample_type = "POSTN",
    replacement_file = "./data/sample_name_replacement.csv",
    geneModels = geneModels,
    bw_dir = "./data/bw/",
    suffix = "up"
  )
}

for (i in c("A1", "A2", "P", "R1", "R2")) {
  top_variable <- get(str_c("top_", i))
  gene_to_plot <- top_variable[top_variable$contrast == "G3BP_Ctrl", "top_down_gene"]
  plot_gene_gviz(
    histoneMark = i,
    gene = gene_to_plot,
    sample_type = "G3BP",
    replacement_file = "./data/sample_name_replacement.csv",
    geneModels = geneModels,
    bw_dir = "./data/bw/",
    suffix = "down"
  )

  gene_to_plot <- top_variable[top_variable$contrast == "GAL1_Ctrl", "top_down_gene"]
  plot_gene_gviz(
    histoneMark = i,
    gene = gene_to_plot,
    sample_type = "GAL1",
    replacement_file = "./data/sample_name_replacement.csv",
    geneModels = geneModels,
    bw_dir = "./data/bw/",
    suffix = "down"
  )

  gene_to_plot <- top_variable[top_variable$contrast == "POSTN_Ctrl", "top_down_gene"]
  plot_gene_gviz(
    histoneMark = i,
    gene = gene_to_plot,
    sample_type = "POSTN",
    replacement_file = "./data/sample_name_replacement.csv",
    geneModels = geneModels,
    bw_dir = "./data/bw/",
    suffix = "down"
  )
}