library(Seurat)
library(SeuratObject)
library(DESeq2)
library(fgsea)
library(data.table)
library(ggplot2)
library(msigdbr)
library(plyr)
library(dplyr)

# differential expression analysis with pseudobulk approach followed by GSEA

## preparing pathways
## retrieve gene sets
cell_sets <- msigdbr(
  species = "Homo sapiens",
  collection = "C8"
) %>%
  dplyr::filter(grepl(
    "^(VANGURP_PANCREATIC_|MURARO_PANCREAS_)", gs_name
  )) %>%
  dplyr::transmute(
    gs_name, gene_symbol,
    Collection = "Islet cell type"
  )

pathway_names <- c(
  "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING",
  "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
  "GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE"
)

UPR_pathways <- msigdbr(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "GO:BP"
) %>%
  dplyr::filter(gs_name %in% pathway_names) %>%
  dplyr::transmute(
    gs_name, gene_symbol,
    Collection = "GO UPR/ISR"
  )

reactome_ids <- c(
  "R-HSA-381070",   # IRE1alpha activates chaperones (XBP1s targets)    50 genes
  "R-HSA-381042",   # PERK regulates gene expression                     32 genes
  "R-HSA-380994",   # ATF4 activates genes in response to ER stress      27 genes
  "R-HSA-381033"   # ATF6 (ATF6-alpha) activates chaperones             12 genes
)

reactome_UPR <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
) %>%
  dplyr::filter(gs_exact_source %in% reactome_ids)

reactome_UPR <- reactome_UPR %>%
  dplyr::transmute(
    gs_name, gene_symbol,
    Collection = "Reactome UPR"
  )


msig_h <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
) %>%
  dplyr::transmute(
    gs_name, gene_symbol,
    Collection = "Hallmark"
  )

pathway_sets <- dplyr::bind_rows(
  cell_sets, UPR_pathways, reactome_UPR, msig_h
) %>%
  dplyr::distinct()

# List for fgsea
pathways <- split(
  pathway_sets$gene_symbol,
  pathway_sets$gs_name
)
pathways <- lapply(pathways, unique)


## differential expression analysis
## =============================================================================
## Pseudobulk DE and GSEA in every cell type: within-sex time responses, then the sex difference
##   One DESeq2 model per cell type:  ~ sex + sex:donor_n + sex:condition
##     (donor nested in sex; every sample is compared with its own donor's 0 h)
##   Within sex:  female / male, 12 h and 48 h vs 0 h
##   Sex contrast: male response - female response (sex x time interaction)
## =============================================================================

ERstress.integrated <- readRDS("ERstress.integrated.rds")
pseudo_bulk <- AggregateExpression(ERstress.integrated, assays = "RNA", return.seurat = T, 
                                   group.by = c("cell.type.final","sex","condition",  "donor"))

skip      <- c("doublets", "epsilon")     # doublets; epsilon pseudobulks too small (median 0.06 M reads)
out_dir   <- "DEG/results"; dir.create(out_dir, showWarnings = FALSE)
pb  <- pseudo_bulk
cnt <- LayerData(pb, assay = "RNA", layer = "counts")
md  <- pb[[]]
gs  <- pathways
cell_types <- setdiff(sort(unique(md$cell.type.final)), skip)
args <- commandArgs(trailingOnly = TRUE); if (length(args)) cell_types <- strsplit(args[1], ",")[[1]]   # optional: run a subset

contrasts <- function(dds) {
  out <- list()
  for (tm in c("12h", "48h")) { tp <- paste0(tm, ".Tg")
  out[[paste0("female_", tm)]] <- results(dds, name = paste0("sexfemale.condition", tp))
  out[[paste0("male_", tm)]]   <- results(dds, name = paste0("sexmale.condition", tp))
  out[[paste0("sex_x_", tm)]]  <- results(dds, contrast = list(paste0("sexmale.condition", tp), paste0("sexfemale.condition", tp)))
  }
  out
}

de_all <- list(); gsea_all <- list(); summ <- list()
for (ct in cell_types) {
  t0 <- Sys.time()
  s <- md$cell.type.final == ct
  m <- droplevels(md[s, ])
  m$sex       <- factor(m$sex, levels = c("female", "male"))
  m$condition <- factor(m$condition, levels = c("0h-Tg", "12h-Tg", "48h-Tg"))
  m$donor_n   <- factor(ave(as.integer(factor(m$donor)), m$sex, FUN = function(v) as.integer(factor(v))))
  dds <- DESeqDataSetFromMatrix(round(as.matrix(cnt[, s])), m, design = ~ sex + sex:donor_n + sex:condition)
  dds <- dds[rowSums(counts(dds) >= 10) >= 1, ]
  dds <- suppressMessages(DESeq(dds, quiet = TRUE))
  for (k in names(cs <- contrasts(dds))) {
    r <- cs[[k]]
    de_all[[paste(ct, k)]] <- data.frame(celltype = ct, contrast = k, gene = rownames(r), log2FC = r$log2FoldChange,
                                         lfcSE = r$lfcSE, stat = r$stat, pvalue = r$pvalue, padj = r$padj)
    st <- sort(setNames(r$stat, rownames(r))[is.finite(r$stat)], decreasing = TRUE)
    set.seed(42)
    f <- suppressWarnings(fgsea(gs, st, minSize = 5, maxSize = 2000))
    f$leadingEdge <- vapply(f$leadingEdge, paste, "", collapse = ",")
    gsea_all[[paste(ct, k)]] <- data.frame(celltype = ct, contrast = k, f)
    summ[[paste(ct, k)]] <- data.frame(celltype = ct, contrast = k, genes = nrow(r), DE_FDR0.05 = sum(r$padj < 0.05, na.rm = TRUE),
                                       up = sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE),
                                       down = sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE),
                                       pct_p0.05 = round(100 * mean(r$pvalue < 0.05, na.rm = TRUE), 1))
  }
  saveRDS(list(de = de_all[grepl(paste0("^", ct, " "), names(de_all))], gsea = gsea_all[grepl(paste0("^", ct, " "), names(gsea_all))],
               summ = summ[grepl(paste0("^", ct, " "), names(summ))]), file.path(out_dir, paste0("ct_", ct, ".rds")))
  message(sprintf("%-12s %d genes, %.1f min", ct, nrow(dds), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
## combine every cell type finished so far
parts <- lapply(list.files(out_dir, "^ct_.*\\.rds$", full.names = TRUE), readRDS)
write.csv(do.call(rbind, unlist(lapply(parts, `[[`, "de"), recursive = FALSE)),   file.path(out_dir, "DE_all_celltypes.csv"), row.names = FALSE)
write.csv(do.call(rbind, unlist(lapply(parts, `[[`, "gsea"), recursive = FALSE)), file.path(out_dir, "GSEA_all_celltypes.csv"), row.names = FALSE)
write.csv(do.call(rbind, unlist(lapply(parts, `[[`, "summ"), recursive = FALSE)), file.path(out_dir, "DE_summary.csv"), row.names = FALSE)




## =============================================================================
## GSEA heatmaps, all gene sets, all cell types
##   Fig A  within-sex response (each sex and timepoint vs its own 0 h)        red = induced, blue = repressed
##   Fig B  sex contrast (male - female response), bulk islets as reference      purple = male > female, green = female > male
##   tile colour = NES; * = FDR below the cut-off; grey = not tested
##   rows = every gene set (top_n = Inf) or the top_n gene sets of each analysis (lowest FDR), grouped by theme
##          (Hallmark, islet cell types, GO UPR, Reactome UPR)
## Input: GSEA_all_celltypes.csv from the DESeq2 + fgsea loop
## Requires: ggplot2, patchwork
## =============================================================================
library(ggplot2)
library(patchwork)

## ---- settings ----------------------------------------------------------------------------
top_n      <- 5                                  # top gene sets per analysis; Inf = all gene sets
upr_always <- TRUE                               # sex contrast: always show the UPR sets (GO and Reactome)
gsea_file  <- "DEG/results/GSEA_all_celltypes.csv"
bulk_file  <- "fgsea_sex_all_times.csv"          # bulk sex x time GSEA; set to NA to leave bulk out
out_dir    <- "DEG/figures"
fdr_within <- 0.05                               # * in the within-sex heatmap
fdr_sex    <- 0.1                                # * in the sex-contrast heatmap
nes_cap    <- 2.5                                # colour scale saturates at +/- this NES
row_h      <- 0.14                               # inches per gene set
cell_order <- c("alpha", "beta", "delta", "PP", "acinar", "ductal", "endothelial", "fibroblast", "immune")
ct_cols    <- c(bulk = "#525252", alpha = "#88BDFF", beta = "#FF8000", delta = "#FFFA0A", PP = "#A65628", acinar = "#F781BF",
                ductal = "#E41A1C", endothelial = "#4DAF4A", fibroblast = "#984EA3", immune = "#BDBDBD")
ann_cols   <- list(condition = c(Tg_12h = "#A6CEE3", Tg_48h = "#1F78B4"), sex = c(female = "#7FC97F", male = "#BEAED4"))
## heatmap colours: the same hues as the dot plots, darker so the tiles are readable
col_up <- "#B2182B"; col_down <- "#2166AC"            # within-sex
col_male <- "#BEAED4"; col_female <- "#7FC97F"        # sex contrast (dot plot)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- data ------------------------------------------------------------------------------------
g <- read.csv(gsea_file)
g$time <- sub(".*_", "", g$contrast)                                  # "12h" / "48h"
g$sex  <- ifelse(grepl("^sex_x", g$contrast), NA, sub("_.*", "", g$contrast))

## ---- readable names ----------------------------------------------------------------------------
pretty_name <- function(p) {
  src <- ifelse(grepl("^MURARO_", p), " (Muraro)", ifelse(grepl("^VANGURP_", p), " (van Gurp)",
                                                          ifelse(grepl("^GOBP_", p), " (GO)", ifelse(grepl("^REACTOME_", p), " (Reactome)", ""))))
  s <- sub("^(HALLMARK|GOBP|REACTOME|MURARO_PANCREAS|VANGURP_PANCREATIC)_", "", p)
  s <- sub("^IRE1ALPHA_ACTIVATES_CHAPERONES$", "IRE1a_activates_chaperones", s)
  s <- sub("^ATF4_ACTIVATES_GENES_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS$", "ATF4_activates_ER-stress_genes", s)
  s <- sub("^ATF6_ATF6_ALPHA_", "ATF6_", s)
  s <- sub("^UNFOLDED_PROTEIN_RESPONSE_UPR$", "UNFOLDED_PROTEIN_RESPONSE", s)
  s <- tolower(gsub("_", " ", s)); s <- paste0(toupper(substring(s, 1, 1)), substring(s, 2))
  fix <- c(mtorc1 = "mTORC1", tnfa = "TNF-alpha", nfkb = "NF-kB", myc = "MYC", v1 = "V1", v2 = "V2", ire1 = "IRE1",
           perk = "PERK", atf6 = "ATF6", e2f = "E2F", g2m = "G2M", il6 = "IL-6", jak = "JAK", stat3 = "STAT3", il2 = "IL-2",
           stat5 = "STAT5", kras = "KRAS", uv = "UV", dna = "DNA", tgf = "TGF", pi3k = "PI3K", akt = "AKT", mtor = "mTOR",
           ire1a = "IRE1a", atf4 = "ATF4", creb3 = "CREB3", `er-stress` = "ER-stress")
  for (f in names(fix)) s <- gsub(paste0("\\b", f, "\\b"), fix[[f]], s, ignore.case = TRUE)
  s <- sub(" mediated unfolded protein response", "-mediated UPR", s)
  s <- sub("Pancreatic polypeptide cell", "PP cell", s)
  s <- sub("IL-6 JAK STAT3", "IL-6/JAK/STAT3", s); s <- sub("IL-2 STAT5", "IL-2/STAT5", s); s <- sub("PI3K AKT mTOR", "PI3K/AKT/mTOR", s)
  s <- sub(" dn$", " down", s); s <- sub("Epithelial mesenchymal", "Epithelial-mesenchymal", s); s <- sub("Wnt beta catenin", "Wnt/beta-catenin", s)
  s <- sub("TGF beta", "TGF-beta", s)
  paste0(s, src)
}

## ---- themes (checked from top to bottom) and row order -------------------------------------------
theme_rules <- c(
  "ER stress"               = "UNFOLDED_PROTEIN|INTEGRATED_STRESS|PROTEIN_SECRETION|IRE1|PERK|ATF4_ACTIVATES|ATF6|CREB3",
  "Inflammation"            = "TNFA|INTERFERON|IL6|IL2|INFLAMMATORY|COMPLEMENT|ALLOGRAFT",
  "Growth signaling"        = "MTORC1|MYC_TARGETS|PI3K",
  "Cell cycle"              = "E2F|G2M|MITOTIC",
  "Cell stress & death"     = "P53|APOPTOSIS|HYPOXIA|UV_RESPONSE|REACTIVE_OXYGEN|DNA_REPAIR",
  "Metabolism"              = "OXIDATIVE|GLYCOLYSIS|FATTY_ACID|CHOLESTEROL|HEME|ADIPOGENESIS|BILE_ACID|XENOBIOTIC|PEROXISOME",
  "Endocrine cell identity" = "ALPHA_CELL|BETA_CELL|DELTA_CELL|POLYPEPTIDE|EPSILON|GAMMA_CELL",
  "Exocrine cell identity"  = "ACINAR|DUCTAL",
  "Endothelial & ECM"       = "ENDOTHELIAL|STROMAL|TGF_BETA|MESENCHYMAL_TRANSITION|ANGIOGENESIS|APICAL_JUNCTION|COAGULATION")
## ER-stress rows in UPR-branch order: whole UPR, IRE1, PERK / ISR, ATF6, other
er_order <- c("Unfolded protein response", "Unfolded protein response (Reactome)", "Cellular response to unfolded protein (GO)",
              "IRE1a activates chaperones (Reactome)", "IRE1-mediated UPR (GO)",
              "PERK regulates gene expression (Reactome)", "ATF4 activates ER-stress genes (Reactome)",
              "PERK-mediated UPR (GO)", "Integrated stress response signaling (GO)",
              "ATF6 activates chaperones (Reactome)", "ATF6 activates chaperone genes (Reactome)", "ATF6-mediated UPR (GO)",
              "CREB3 factors activate genes (Reactome)", "Protein secretion")
## Endocrine cell-identity rows by cell type: alpha, beta, delta, epsilon, PP (gamma = PP cells)
celltype_order <- c("Alpha cell (Muraro)", "Alpha cell (van Gurp)",
                    "Beta cell (Muraro)", "Beta cell (van Gurp)", "Pancreas beta cells",
                    "Delta cell (Muraro)", "Delta cell (van Gurp)",
                    "Epsilon cell (Muraro)",
                    "PP cell (Muraro)", "Gamma cell (van Gurp)")
## themes whose rows follow a fixed order (all other themes: by mean NES)
fixed_orders <- list("ER stress" = er_order, "Endocrine cell identity" = celltype_order)

prepare_rows <- function(d, order_by) {
  d$theme <- "Other pathways"
  for (t in rev(names(theme_rules))) d$theme[grepl(theme_rules[[t]], d$pathway)] <- t
  d$theme <- factor(d$theme, c(names(theme_rules), "Other pathways"))
  d$label <- pretty_name(d$pathway)
  ro <- aggregate(NES ~ label + theme, d[order_by & !is.na(d$NES), ], mean)       # by theme, then mean NES
  ro <- ro[order(ro$theme, -ro$NES), ]
  for (th in names(fixed_orders)) {                                      # fixed order where defined; rows not listed go last
    i  <- which(ro$theme == th); if (!length(i)) next
    rk <- match(ro$label[i], fixed_orders[[th]]); rk[is.na(rk)] <- 99 + seq_along(rk[is.na(rk)])
    ro[i, ] <- ro[i, ][order(rk), ]
  }
  d$label <- factor(d$label, rev(ro$label)); d$theme <- droplevels(d$theme)
  d
}

## top gene sets: the top_n of each analysis by FDR (ties broken alphabetically), pooled across analyses
top_sets <- function(d, n) {
  d <- d[!is.na(d$padj), ]; d <- d[order(d$padj, d$pathway), ]
  unique(unlist(lapply(split(d$pathway, paste(d$celltype, d$contrast)), head, n)))
}
suffix   <- if (is.finite(top_n)) paste0("top", top_n) else "all_sets"

## ---- heatmap with annotation bars ---------------------------------------------------------------
heatmap_figure <- function(d, cols, ann_vars, low, high, low_lab, high_lab, fdr_cut, title, subtitle, file, width) {
  x_lim <- range(cols$x) + c(-0.5, 0.5)
  bar <- function(var, order) {
    vals <- if (var == "celltype") ct_cols[levels(cols$celltype)] else ann_cols[[var]]
    ggplot(cols, aes(x, var, fill = .data[[var]])) +
      geom_tile(width = 1, height = 1, colour = "white", linewidth = 0.3) +
      scale_fill_manual(values = vals, name = var, guide = guide_legend(order = order, nrow = 1, title.position = "left")) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) + scale_y_discrete(position = "right", expand = c(0, 0)) +
      theme_void(base_size = 10) +
      theme(axis.text.y.right = element_text(face = "bold", size = 7.5, hjust = 0, margin = margin(l = 3)), plot.margin = margin(0, 4, 1, 4))
  }
  d$star <- ifelse(!is.na(d$padj) & d$padj < fdr_cut, "*", "")
  d$star_col <- ifelse(!is.na(d$NES) & abs(d$NES) > 0.45 * nes_cap, "white", "black")      # white on dark tiles
  hm <- ggplot(d, aes(x, label)) +
    geom_tile(aes(fill = pmax(pmin(NES, nes_cap), -nes_cap)), width = 1, height = 1, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = star, colour = star_col), size = 2.6, vjust = 0.75) + scale_colour_identity() +
    facet_grid(theme ~ ., scales = "free_y", space = "free_y", labeller = labeller(theme = label_wrap_gen(14))) +
    scale_x_continuous(limits = x_lim, expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(low = low, mid = "white", high = high, limits = c(-nes_cap, nes_cap), na.value = "grey85",
                         breaks = c(-nes_cap, 0, nes_cap), labels = c(low_lab, "0", high_lab), name = "NES",
                         guide = guide_colourbar(order = 4, direction = "horizontal", barwidth = unit(5, "cm"), barheight = unit(0.3, "cm"),
                                                 title.position = "left", title.vjust = 1)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
          axis.text.y = element_text(size = 6.8), panel.border = element_rect(colour = "grey40", linewidth = 0.3),
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 7.2), strip.background = element_rect(fill = "grey92"),
          panel.spacing.y = unit(0.2, "lines"), plot.margin = margin(0, 4, 6, 4))
  bars <- lapply(seq_along(ann_vars), function(i) bar(ann_vars[i], i))
  fig <- wrap_plots(c(bars, list(hm)), ncol = 1) +
    plot_layout(heights = c(rep(0.15, length(bars)), row_h * nlevels(d$label) + 0.05 * nlevels(d$theme)), guides = "collect") +
    plot_annotation(title = title, subtitle = subtitle,
                    theme = theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 7.6))) &
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.justification = "left",
          legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 7.3), legend.title = element_text(size = 7.8, face = "bold"),
          legend.spacing.y = unit(0.02, "cm"))
  h <- 2.6 + row_h * nlevels(d$label) + 0.05 * nlevels(d$theme) + 0.15 * length(bars)
  ggsave(file.path(out_dir, paste0(file, ".pdf")), fig, width = width, height = h, device = if (capabilities("cairo")) cairo_pdf else "pdf")
  ggsave(file.path(out_dir, paste0(file, ".png")), fig, width = width, height = h, dpi = 300, bg = "white")
}





## ---- dot plot with annotation bars (used for Fig B) ------------------------------------------------
dot_figure <- function(d, cols, ann_vars, low, high, low_lab, high_lab, fdr_cut, title, subtitle, file, width, row_h = 0.2) {
  x_lim <- range(cols$x) + c(-0.55, 0.55)
  bar <- function(var, order) {
    vals <- if (var == "celltype") ct_cols[levels(cols$celltype)] else ann_cols[[var]]
    ggplot(cols, aes(x, var, fill = .data[[var]])) + geom_tile(width = 1, height = 1, colour = "white", linewidth = 0.4) +
      scale_fill_manual(values = vals, name = var, guide = guide_legend(order = order, nrow = 1, title.position = "left")) +
      scale_x_continuous(limits = x_lim, expand = c(0, 0)) + scale_y_discrete(position = "right", expand = c(0, 0)) + theme_void(base_size = 10) +
      theme(axis.text.y.right = element_text(face = "bold", size = 7.5, hjust = 0, margin = margin(l = 3)), plot.margin = margin(0, 4, 1, 4))
  }
  d <- d[!is.na(d$NES), ]                                              # gene sets without a result get no dot
  d$sig <- !is.na(d$padj) & d$padj < fdr_cut
  dots <- ggplot(d, aes(x, label)) +
    geom_point(aes(size = pmin(-log10(padj), 6), fill = pmax(pmin(NES, nes_cap), -nes_cap), colour = sig), shape = 21, stroke = 0.5) +
    facet_grid(theme ~ ., scales = "free_y", space = "free_y", labeller = labeller(theme = label_wrap_gen(14))) +
    scale_x_continuous(breaks = cols$x, labels = NULL, limits = x_lim, expand = c(0, 0)) +
    scale_fill_gradient2(low = low, mid = "white", high = high, limits = c(-nes_cap, nes_cap), breaks = c(-nes_cap, 0, nes_cap),
                         labels = c(low_lab, "0", high_lab), name = "NES",
                         guide = guide_colourbar(order = 4, direction = "horizontal", barwidth = unit(4.5, "cm"), barheight = unit(0.3, "cm"),
                                                 title.position = "left", title.vjust = 1)) +
    scale_size_continuous(range = c(0.4, 3.3), limits = c(0, 6), breaks = c(1, 2, 4, 6), labels = c("1", "2", "4", ">=6"),
                          name = expression(-log[10]~FDR), guide = guide_legend(order = 5, nrow = 1, title.position = "left")) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), breaks = c("TRUE", "FALSE"), labels = c(paste("FDR <", fdr_cut), "n.s."),
                        name = NULL, guide = guide_legend(order = 6, nrow = 1, override.aes = list(size = 3, fill = "white"))) +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey93"), axis.ticks.x = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 7.2), strip.background = element_rect(fill = "grey92"),
          axis.text.y = element_text(size = 7.4), panel.spacing.y = unit(0.25, "lines"), plot.margin = margin(0, 4, 6, 4))
  bars <- lapply(seq_along(ann_vars), function(i) bar(ann_vars[i], i))
  fig <- wrap_plots(c(bars, list(dots)), ncol = 1) +
    plot_layout(heights = c(rep(0.16, length(bars)), row_h * nlevels(d$label) + 0.06 * nlevels(d$theme)), guides = "collect") +
    plot_annotation(title = title, subtitle = subtitle,
                    theme = theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 7.6))) &
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left", legend.justification = "left",
          legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 7.3), legend.title = element_text(size = 7.8, face = "bold"),
          legend.spacing.y = unit(0.02, "cm"))
  h <- 3.6 + row_h * nlevels(d$label) + 0.06 * nlevels(d$theme) + 0.16 * length(bars)
  ggsave(file.path(out_dir, paste0(file, ".pdf")), fig, width = width, height = h, device = if (capabilities("cairo")) cairo_pdf else "pdf")
  ggsave(file.path(out_dir, paste0(file, ".png")), fig, width = width, height = h, dpi = 300, bg = "white")
}



## ---------------------------------------------------------------------------
## Fig A. Within-sex response
## ---------------------------------------------------------------------------
a <- g[!grepl("^sex_x", g$contrast) & g$celltype %in% cell_order, ]
if (is.finite(top_n)) a <- a[a$pathway %in% top_sets(a, top_n), ]
a <- prepare_rows(a, rep(TRUE, nrow(a)))
gap_sex <- 0.3; gap_cell <- 0.9                         # narrow gap between sexes, wide gap between cell types
a$x <- (match(a$celltype, cell_order) - 1) * (4 + gap_sex + gap_cell) + (a$sex == "male") * (2 + gap_sex) + (a$time == "48h") + 1
ca <- unique(a[, c("celltype", "sex", "time", "x")]); ca$celltype <- factor(ca$celltype, cell_order); ca$condition <- paste0("Tg_", ca$time)
heatmap_figure(a, ca, c("celltype", "condition", "sex"), col_down, col_up, "repressed", "induced", fdr_within,
               paste0("Thapsigargin response within each sex: ", if (is.finite(top_n)) sprintf("top %d gene sets", top_n) else "all gene sets", ", all islet cell types"),
               sprintf("%sNES from fgsea of each sex and timepoint vs its own 0 h (DESeq2 per cell type, donor nested in sex; Wald statistic). * FDR < %s; grey = not tested.\nFemale 48 h: one donor was not sampled at this timepoint.",
                       if (is.finite(top_n)) sprintf("Rows: the top %d gene sets of each analysis (lowest FDR), %d in total. ", top_n, nlevels(a$label)) else "", fdr_within),
               paste0("heatmap_within_sex_", suffix), width = 11)

## ---------------------------------------------------------------------------
## Fig B. Sex contrast (male - female response)
## ---------------------------------------------------------------------------

b <- g[grepl("^sex_x", g$contrast) & g$celltype %in% cell_order, ]
if (is.finite(top_n)) {
  keep <- top_sets(b, top_n)
  if (upr_always) keep <- unique(c(keep, b$pathway[grepl("^(GOBP|REACTOME)_|^HALLMARK_UNFOLDED_PROTEIN", b$pathway)]))
  b <- b[b$pathway %in% keep, ]
}
b <- b[, c("celltype", "time", "pathway", "NES", "padj")]
groups <- cell_order

b <- prepare_rows(b, b$celltype != "bulk")
gap <- 0.8
b$x <- (match(b$celltype, groups) - 1) * (2 + gap) + (b$time == "48h") + 1
cb <- unique(b[, c("celltype", "time", "x")]); cb$celltype <- factor(cb$celltype, groups); cb$condition <- paste0("Tg_", cb$time)
dot_figure(b, cb, c("celltype", "condition"), col_female, col_male, "female > male", "male > female", fdr_sex,
           paste0("Sex differences in the thapsigargin response: ", if (is.finite(top_n)) sprintf("top %d gene sets", top_n) else "all gene sets", ", all islet cell types"),
           paste0("fgsea of the sex x time interaction (male - female change from each donor's 0 h); DESeq2 per cell type, donor nested in sex.\n",
                  if (is.finite(top_n)) sprintf("Rows: the top %d gene sets of each cell type and timepoint%s, %d in total. ", top_n,
                                                if (upr_always) ", plus all UPR sets (ordered by UPR branch)" else "", nlevels(b$label)) else "",
                  sprintf("Black outline = FDR < %s.", fdr_sex),
                  if ("bulk" %in% groups) "\nBulk islets shown for reference; gene sets not tested in bulk (Reactome UPR sets, GO parent term) have no bulk dots." else ""),
           paste0("dotplot_sex_contrast_", suffix), width = 9.2)

