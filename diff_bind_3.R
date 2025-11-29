library(DESeq2)
library(stringr)
library(limma)

df <- read.table("./data/hg38_promoters_chr1-22_X_3kb.tsv", header = T, sep = "\t", check.name = F,               
                 comment.char = "",
                 quote = "'")
rownames(df) <- str_c("pr", seq(1, nrow(df)), sep = "_")
tmp <- read.csv("./data/sample_name_replacement.csv", row.names = 1)
identical(tmp$x, colnames(df))
colnames(df) <- tmp$replacement # if previous code is TRUE
df_meta <- df[,c(1:3)]

write.csv(df_meta, file = "./data/pr_metadata.csv")
write.csv(df, file = "./data/hg38_promoters_chr1-22_X_3kb_clean.tsv.csv")

run_promoter_limma <- function(df,
                               histoneMark = "A1",
                               save_dir = "./") {

  message("===== Running promoter LIMMA for histone mark ", histoneMark, " =====")

  # 1. Folders ---------------------------------------------------------------
  data_save     <- str_c(save_dir, "data/promoter/", histoneMark, "/")
  analysis_save <- str_c(save_dir, "analysis/promoter/", histoneMark, "/")
  dir.create(data_save, recursive = TRUE, showWarnings = FALSE)
  dir.create(analysis_save, recursive = TRUE, showWarnings = FALSE)

  # 2. Select samples with this histone mark ---------------------------------
  pattern <- str_c("_", histoneMark, "_")
  df_sel <- df[, grep(pattern, colnames(df))]
  if (ncol(df_sel) == 0)
    stop("No samples found for histone mark: ", histoneMark)

  message("Selected ", ncol(df_sel), " samples.")

  # 3. Parse sample name: group_histone_rep ----------------------------------
  # Example: Ctrl_A1_r1 → group=Ctrl, histone=A1, rep=r1
  parsed <- str_match(colnames(df_sel), "^(.+?)_(.+?)_(r\\d+)$")

  coldata <- data.frame(
    sample  = colnames(df_sel),
    group   = parsed[, 2],  # Ctrl, GAL1, POSTN, G3BP
    histone = parsed[, 3],  # A1, A2, C, P, R1, R2
    rep     = parsed[, 4]   # r1, r2, r3, ...
  )

  # Ensure factor levels: Ctrl always first if present
  if ("Ctrl" %in% coldata$group) {
    coldata$group <- factor(coldata$group,
                            levels = c("Ctrl", setdiff(sort(unique(coldata$group)), "Ctrl")))
  } else {
    coldata$group <- factor(coldata$group)
  }

  rownames(coldata) <- coldata$sample

  # 4. Filter low expression --------------------------------------------------
  keep <- rowMeans(df_sel) > 2
  df_f <- df_sel[keep, ]

  pdf(str_c(analysis_save, "hist_raw_values.pdf"))
  hist(as.numeric(as.matrix(df_f)), breaks = 50,
       main = "Raw promoter signal", xlab = "Raw signal")
  dev.off()

  # 5. log10 normalization ----------------------------------------------------
  df_log <- log10(df_f + 1)

  pdf(str_c(analysis_save, "hist_log10_values.pdf"))
  hist(as.numeric(as.matrix(df_log)), breaks = 50,
       main = "Log10 promoter signal", xlab = "log10(x+1)")
  dev.off()

  # 6. LIMMA design -----------------------------------------------------------
  design <- model.matrix(~ 0 + group, data = coldata)
  colnames(design) = levels(coldata$group)

  fit <- lmFit(df_log, design)

  # 7. Contrasts: each group vs Ctrl -----------------------------------------
  ctrl_group <- levels(coldata$group)[1]
  contrasts <- setdiff(levels(coldata$group), ctrl_group)

  contrast.matrix <- makeContrasts(
    contrasts = sapply(contrasts,
                       function(g) str_c(g, "-", ctrl_group)),
    levels = design
  )

  fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))

  # 8. Save results -----------------------------------------------------------
  for (coef_name in colnames(contrast.matrix)) {
    out <- str_c(data_save, "bw_limma_", coef_name, ".csv")
    write.csv(topTable(fit2, coef = coef_name, number = Inf),
              out, row.names = TRUE)
    message("Saved: ", out)
  }

  message("===== DONE: ", histoneMark, " =====")

  return(list(
    coldata = coldata,
    design = design,
    contrasts = contrast.matrix,
    fit = fit2
  ))
}

res_A1 = run_promoter_limma(df, histoneMark = "A1")
res_A2 = run_promoter_limma(df, "A2")
res_P  = run_promoter_limma(df, "P")
res_R1 = run_promoter_limma(df, "R1")
res_R2 = run_promoter_limma(df, "R2")

########################################################
# limma_result:   output from run_promoter_limma()
# histoneMark:    "A1", "A2", "C", "P", "R1", "R2"
# save_dir:       project root (default "./")
# pr_metadata:    shared promoter metadata file
# ref_bed:        promoter BED annotation
########################################################
annotate_all_contrasts <- function(limma_result,
                                   histoneMark  = "A1",
                                   save_dir     = "./",
                                   pr_metadata  = "./data/pr_metadata.csv",
                                   ref_bed      = "./ref/hg38_promoters_chr1-22_X_3kb.bed",
                                   logFC_min    = 0,
                                   p_cutoff     = 0.05) {

  message("===== Annotating all contrasts for histone mark ", histoneMark, " =====")

  # folders
  data_save     <- str_c(save_dir, "data/promoter/", histoneMark, "/")
  analysis_save <- str_c(save_dir, "analysis/promoter/", histoneMark, "/")
  dir.create(analysis_save, recursive = TRUE, showWarnings = FALSE)

  # 1. Load promoter BED annotation -----------------------------------------
  tmp <- read.table(ref_bed, header = FALSE)
  colnames(tmp) <- c("chr", "start", "end", "SYMBOL", "score", "strand")

  # 2. Load *shared* promoter metadata --------------------------------------
  df_meta <- read.csv(pr_metadata, header = TRUE)
  colnames(df_meta) <- c("name", "chr", "start", "end")

  # 3. Merge BED + promoter region metadata ---------------------------------
  annot <- merge(tmp, df_meta, by = c("chr", "start", "end"))

  # 4. Get all contrasts defined in LIMMA result -----------------------------
  contrast_names <- colnames(limma_result$contrasts)

  summary_list <- list()

  for (coef_name in contrast_names) {

    message("Processing contrast: ", coef_name)

    # 4a. full differential table
    df_limma <- topTable(limma_result$fit,
                         coef   = coef_name,
                         number = Inf)

    # 4b. Identify up/down regulated promoters
    df_up <- df_limma[df_limma$logFC >  logFC_min & df_limma$P.Value < p_cutoff, ]
    df_dn <- df_limma[df_limma$logFC < -logFC_min & df_limma$P.Value < p_cutoff, ]

    df_up$region_id <- rownames(df_up)
    df_dn$region_id <- rownames(df_dn)

    # 4c. Annotate with gene symbols using "name"
    up_annot  <- merge(annot, df_up, by.x = "name", by.y = "region_id")
    dn_annot  <- merge(annot, df_dn, by.x = "name", by.y = "region_id")

    # ensure filename safe
    safe_name <- gsub("-", "_", coef_name)

    # 4d. Save annotated results
    up_file <- str_c(data_save, safe_name, "_UP_genes.csv")
    dn_file <- str_c(data_save, safe_name, "_DOWN_genes.csv")

    write.csv(up_annot, up_file, row.names = FALSE)
    write.csv(dn_annot, dn_file, row.names = FALSE)

    message("  UP:   ", nrow(up_annot),  " promoters → ", up_file)
    message("  DOWN: ", nrow(dn_annot), " promoters → ", dn_file)

    # 4e. Return summary
    summary_list[[coef_name]] <- list(
      up_n   = nrow(up_annot),
      down_n = nrow(dn_annot),
      up_file = up_file,
      down_file = dn_file
    )
  }

  message("===== Done annotating for ", histoneMark, " =====")
  invisible(summary_list)
}

summary_A1 <- annotate_all_contrasts(res_A1,
                                     histoneMark = "A1")
summary_A2 <- annotate_all_contrasts(res_A2,
                                     histoneMark = "A2")
summary_P <- annotate_all_contrasts(res_P,
                                     histoneMark = "P")
summary_R1 <- annotate_all_contrasts(res_R1,
                                     histoneMark = "R1")
summary_R2 <- annotate_all_contrasts(res_R2,
                                     histoneMark = "R2")

