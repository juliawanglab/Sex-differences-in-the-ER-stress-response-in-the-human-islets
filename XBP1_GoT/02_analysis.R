# Clean up the dataframes
HP23166_0h <- read.table("HP23166_0h_S1_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(HP23166_0h) <- sub(";.*$", "", HP23166_0h$BC)
HP23166_0h_1 <- HP23166_0h[(HP23166_0h$WT.calls!=0)|(HP23166_0h$MUT.calls!=0),]
HP23166_0h_1 <- HP23166_0h_1[,c("WT.calls","MUT.calls")]
write.table(HP23166_0h_1,"XBP1_GoT/HP23166_0h.csv",quote=F,sep=",")

HP23166_12h <- read.table("HP23166_12h_S2_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(HP23166_12h) <- sub(";.*$", "", HP23166_12h$BC)
HP23166_12h_1 <- HP23166_12h[(HP23166_12h$WT.calls!=0)|(HP23166_12h$MUT.calls!=0),]
HP23166_12h_1 <- HP23166_12h_1[,c("WT.calls","MUT.calls")]
write.table(HP23166_12h_1,"XBP1_GoT/HP23166_12h.csv",quote=F,sep=",")

SAMN35848421_0h <- read.table("SAMN35848421_0h_S3_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN35848421_0h) <- sub(";.*$", "", SAMN35848421_0h$BC)
SAMN35848421_0h_1 <- SAMN35848421_0h[(SAMN35848421_0h$WT.calls!=0)|(SAMN35848421_0h$MUT.calls!=0),]
SAMN35848421_0h_1 <- SAMN35848421_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN35848421_0h_1,"XBP1_GoT/SAMN35848421_0h.csv",quote=F,sep=",")

SAMN35848421_12h <- read.table("SAMN35848421_12h_S4_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN35848421_12h) <- sub(";.*$", "", SAMN35848421_12h$BC)
SAMN35848421_12h_1 <- SAMN35848421_12h[(SAMN35848421_12h$WT.calls!=0)|(SAMN35848421_12h$MUT.calls!=0),]
SAMN35848421_12h_1 <- SAMN35848421_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN35848421_12h_1,"XBP1_GoT/SAMN35848421_12h.csv",quote=F,sep=",")

SAMN35848421_48h <- read.table("SAMN35848421_48h_S5_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN35848421_48h) <- sub(";.*$", "", SAMN35848421_48h$BC)
SAMN35848421_48h_1 <- SAMN35848421_48h[(SAMN35848421_48h$WT.calls!=0)|(SAMN35848421_48h$MUT.calls!=0),]
SAMN35848421_48h_1 <- SAMN35848421_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN35848421_48h_1,"XBP1_GoT/SAMN35848421_48h.csv",quote=F,sep=",")

SAMN36705973_0h <- read.table("SAMN36705973_0h_S6_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36705973_0h) <- sub(";.*$", "", SAMN36705973_0h$BC)
SAMN36705973_0h_1 <- SAMN36705973_0h[(SAMN36705973_0h$WT.calls!=0)|(SAMN36705973_0h$MUT.calls!=0),]
SAMN36705973_0h_1 <- SAMN36705973_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36705973_0h_1,"XBP1_GoT/SAMN36705973_0h.csv",quote=F,sep=",")

SAMN36705973_12h <- read.table("SAMN36705973_12h_S7_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36705973_12h) <- sub(";.*$", "", SAMN36705973_12h$BC)
SAMN36705973_12h_1 <- SAMN36705973_12h[(SAMN36705973_12h$WT.calls!=0)|(SAMN36705973_12h$MUT.calls!=0),]
SAMN36705973_12h_1 <- SAMN36705973_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36705973_12h_1,"XBP1_GoT/SAMN36705973_12h.csv",quote=F,sep=",")

SAMN36705973_48h <- read.table("SAMN36705973_48h_S8_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36705973_48h) <- sub(";.*$", "", SAMN36705973_48h$BC)
SAMN36705973_48h_1 <- SAMN36705973_48h[(SAMN36705973_48h$WT.calls!=0)|(SAMN36705973_48h$MUT.calls!=0),]
SAMN36705973_48h_1 <- SAMN36705973_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36705973_48h_1,"XBP1_GoT/SAMN36705973_48h.csv",quote=F,sep=",")

SAMN36823227_0h <- read.table("SAMN36823227_0h_S9_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36823227_0h) <- sub(";.*$", "", SAMN36823227_0h$BC)
SAMN36823227_0h_1 <- SAMN36823227_0h[(SAMN36823227_0h$WT.calls!=0)|(SAMN36823227_0h$MUT.calls!=0),]
SAMN36823227_0h_1 <- SAMN36823227_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36823227_0h_1,"XBP1_GoT/SAMN36823227_0h.csv",quote=F,sep=",")

SAMN36823227_12h <- read.table("SAMN36823227_12h_S10_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36823227_12h) <- sub(";.*$", "", SAMN36823227_12h$BC)
SAMN36823227_12h_1 <- SAMN36823227_12h[(SAMN36823227_12h$WT.calls!=0)|(SAMN36823227_12h$MUT.calls!=0),]
SAMN36823227_12h_1 <- SAMN36823227_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36823227_12h_1,"XBP1_GoT/SAMN36823227_12h.csv",quote=F,sep=",")

SAMN36823227_48h <- read.table("SAMN36823227_48h_S11_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN36823227_48h) <- sub(";.*$", "", SAMN36823227_48h$BC)
SAMN36823227_48h_1 <- SAMN36823227_48h[(SAMN36823227_48h$WT.calls!=0)|(SAMN36823227_48h$MUT.calls!=0),]
SAMN36823227_48h_1 <- SAMN36823227_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN36823227_48h_1,"XBP1_GoT/SAMN36823227_48h.csv",quote=F,sep=",")

SAMN37871873_0h <- read.table("SAMN37871873_0h_S12_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN37871873_0h) <- sub(";.*$", "", SAMN37871873_0h$BC)
SAMN37871873_0h_1 <- SAMN37871873_0h[(SAMN37871873_0h$WT.calls!=0)|(SAMN37871873_0h$MUT.calls!=0),]
SAMN37871873_0h_1 <- SAMN37871873_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN37871873_0h_1,"XBP1_GoT/SAMN37871873_0h.csv",quote=F,sep=",")

SAMN37871873_12h <- read.table("SAMN37871873_12h_S13_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN37871873_12h) <- sub(";.*$", "", SAMN37871873_12h$BC)
SAMN37871873_12h_1 <- SAMN37871873_12h[(SAMN37871873_12h$WT.calls!=0)|(SAMN37871873_12h$MUT.calls!=0),]
SAMN37871873_12h_1 <- SAMN37871873_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN37871873_12h_1,"XBP1_GoT/SAMN37871873_12h.csv",quote=F,sep=",")

SAMN37871873_48h <- read.table("SAMN37871873_48h_S14_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN37871873_48h) <- sub(";.*$", "", SAMN37871873_48h$BC)
SAMN37871873_48h_1 <- SAMN37871873_48h[(SAMN37871873_48h$WT.calls!=0)|(SAMN37871873_48h$MUT.calls!=0),]
SAMN37871873_48h_1 <- SAMN37871873_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN37871873_48h_1,"XBP1_GoT/SAMN37871873_48h.csv",quote=F,sep=",")

SAMN39523303_0h <- read.table("SAMN39523303_0h_S15_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN39523303_0h) <- sub(";.*$", "", SAMN39523303_0h$BC)
SAMN39523303_0h_1 <- SAMN39523303_0h[(SAMN39523303_0h$WT.calls!=0)|(SAMN39523303_0h$MUT.calls!=0),]
SAMN39523303_0h_1 <- SAMN39523303_0h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN39523303_0h_1,"XBP1_GoT/SAMN39523303_0h.csv",quote=F,sep=",")

SAMN39523303_12h <- read.table("SAMN39523303_12h_S16_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN39523303_12h) <- sub(";.*$", "", SAMN39523303_12h$BC)
SAMN39523303_12h_1 <- SAMN39523303_12h[(SAMN39523303_12h$WT.calls!=0)|(SAMN39523303_12h$MUT.calls!=0),]
SAMN39523303_12h_1 <- SAMN39523303_12h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN39523303_12h_1,"XBP1_GoT/SAMN39523303_12h.csv",quote=F,sep=",")

SAMN39523303_48h <- read.table("SAMN39523303_48h_S17_myGoT.summTable.concat.umi_collapsed.txt",header=T,sep="\t")
rownames(SAMN39523303_48h) <- sub(";.*$", "", SAMN39523303_48h$BC)
SAMN39523303_48h_1 <- SAMN39523303_48h[(SAMN39523303_48h$WT.calls!=0)|(SAMN39523303_48h$MUT.calls!=0),]
SAMN39523303_48h_1 <- SAMN39523303_48h_1[,c("WT.calls","MUT.calls")]
write.table(SAMN39523303_48h_1,"XBP1_GoT/SAMN39523303_48h.csv",quote=F,sep=",")


# exporting metadata from the scRNA-seq data
library(Seurat)
ERstress.integrated <- readRDS("ERstress.integrated.rds")
meta.data <- ERstress.integrated@meta.data

# prepare the individual GoT results for merging with meta.data
## Rename rownames based on filename to incorporate sample name and write out new CSVs

# Directory containing the input CSVs
dir_path <- "XBP1_GoT"

# Identify all CSV files
files <- list.files(path = dir_path, pattern = "\\.csv$", full.names = TRUE)

for (f in files) {
  
  # Extract base filename without extension
  base <- sub("\\.csv$", "", basename(f))
  
  # Construct prefix: e.g. "SAMN39523303_48h_Tg_"
  prefix <- paste0(base, "_Tg_")
  
  # Read CSV with rownames in first column
  df <- read.csv(
    f,
    header = TRUE,
    row.names = 1  
    )
  
  # Modify rownames: PREFIX + old_rowname + "-1"
  rownames(df) <- paste0(prefix, rownames(df), "-1")
  
  # Output filename
  out_file <- file.path(
    dir_path,
    paste0(base, "_renamed.csv")
  )
  
  # Write CSV
  write.csv(
    df,
    file = out_file,
    quote = FALSE
  )
}

## append all the files together 
dir_path <- "XBP1_GoT"

# Identify all CSV files
files <- list.files(path = dir_path, pattern = "\\_renamed.csv$", full.names = TRUE)
combined <- data.frame()
for (f in files) {
  df <- read.csv(
    f,
    header = TRUE,
    row.names = 1  
  )
  combined <- rbind(combined,df)
}
write.table(combined,file.path(dir_path,"combined.csv"),quote=F,sep=",")

## put everything together
# Left outer join: keeps all rows from df1 and matching rows from df2 (fills NAs for non-matches)
meta.data$cell <- rownames(meta.data)
combined$cell <- rownames(combined)
merged <- merge(meta.data, combined, by = "cell", all.x = TRUE) # by=0 is shorthand for by="row.names"
rownames(merged) = merged$cell

# compute normalized counts for uXBP1 and sXBP1, using the LogNormalize method as in Seurat
## WT call corresponding to uXBP1, MUT call corresponding to sXBP1
merged$norm.WT.calls <- log1p(merged$WT.calls/merged$nCount_RNA * 10000)
merged$norm.MUT.calls <- log1p(merged$MUT.calls/merged$nCount_RNA * 10000)

merged_1 <- na.omit(merged)
write.table(merged,"meta.data.GoT.merged.csv",quote=F,sep=",")

ERstress.integrated[["uXBP1"]] = merged$WT.calls
ERstress.integrated[["sXBP1"]] = merged$MUT.calls
ERstress.integrated[["norm.uXBP1"]] = merged$norm.WT.calls
ERstress.integrated[["norm.sXBP1"]] = merged$norm.MUT.calls

## filter the Seurat data to retain only those with corresponding to GoT data
ERstress.integrated[["CellName"]] = colnames(ERstress.integrated)
ERstress.integrated.sub <- subset(ERstress.integrated, subset = CellName %in% rownames(merged_1) )


# Merge with pseudotime data
beta <- subset(ERstress.integrated.sub,idents="beta")
beta_sce <- readRDS("beta_sub_sce_filtered_MNN.rds")
pseudotime <- data.frame(colnames(beta_sce),beta_sce$slingPseudotime_1)
common_cells <- intersect(beta$CellName, pseudotime$colnames.beta_sce.)
beta.sub <- subset(beta, subset = CellName %in% common_cells)
rownames(pseudotime) = pseudotime$colnames.beta_sce.
pseudotime.sub <- pseudotime[colnames(beta.sub),]
beta.sub$pseudotime <- pseudotime.sub$beta_sce.slingPseudotime_1

# Extracting INS transcript info
INS_UMI <- FetchData(beta.sub, vars="INS",layer = "counts")
INS_data <- FetchData(beta.sub,vars="INS",layer="data")
INS <- cbind(INS_UMI,INS_data)
colnames(INS) = c("INS","norm.INS")
beta.sub[["INS"]] = INS$INS
beta.sub[["norm.INS"]] = INS$norm.INS



# plotting
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
meta.beta <- beta.sub@meta.data

## QC
meta.beta$GoT_sum <- meta.beta$uXBP1+meta.beta$sXBP1
ggplot(meta.beta, aes(XBP1, GoT_sum)) +
  geom_point( color =  "#1f78b4") +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "firebrick",
  ) +
  stat_cor(aes(label = ..r.label..), color = "black", size = 4) +
  stat_regline_equation(aes(label = ..eq.label..), color = "black", size = 4, label.y.npc = 0.9) +
  labs(
    x = "XBP1 expression",
    y = "GoT total counts"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

## boxplot
scales <- brewer.pal(7,"Accent")[1:2]
ggplot(meta.beta, aes(x = condition, y = norm.sXBP1, fill=sex,color=sex)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.8,
    alpha = 0.75) +
  scale_fill_manual(values = scales)+
  scale_color_manual(values=scales)+
  theme_classic() +
  theme(
    axis.title=element_text(face="bold",size=20),
    axis.text=element_text(face="bold",size=20)
  )


## Plotting regression
scales <- brewer.pal(7,"Accent")[1:2]

ggplot(meta.beta, aes(x = pseudotime, y = norm.uXBP1, color=sex)) +
  geom_smooth() +
  scale_color_manual(values=scales)+
  ylim(0, 2) +
  theme_classic() 

ggplot(meta.beta, aes(x = pseudotime, y = norm.sXBP1, color=sex)) +
  geom_smooth() +
  scale_color_manual(values=scales)+
  ylim(0, 2) +
  theme_classic()

ggplot(meta.beta, aes(x = norm.sXBP1, y = norm.INS, color=sex)) +
  geom_smooth() +
  scale_color_manual(values=scales)+
  theme_classic() 

## Statistical test of transcriptome-wide association with sXBP1,  donor effect as random effect, restricted to cells with higher than median expression of sXBP1

# model:
#gene_expression ~ norm.sXBP1 * sex +
#                  condition * sex +
#                  (1 | donor)


library(Seurat)
library(Matrix)
library(dplyr)
library(lme4)
library(lmerTest)

# ============================================================
# SET OBJECT
# ============================================================

seu <- Beta


# ============================================================
# MIXED-EFFECT MODEL
#
# Model:
#
# gene_expression ~ norm.sXBP1 * sex +
#                   condition * sex +
#                   (1 | donor)
#
# Female = reference sex
#
# quantile_cutoff:
#   0.50 = upper 50% (above pooled median)
#   0.60 = upper 40%
#   0.70 = upper 30%
#
# Expression:
#   Seurat RNA "data" layer
#   log-normalized expression
#
# NO rank transformation
# NO scaling of genes
# NO scaling of sXBP1
# ============================================================

run_high_sXBP1_mixed <- function(
    seu,
    quantile_cutoff = 0.50,
    assay = "RNA",
    layer = "data",
    min_pct = 0.05,
    genes = NULL,
    n_cores = 1,
    verbose_every = 250
) {
  
  
  # ==========================================================
  # 1. Metadata
  # ==========================================================
  
  meta <- seu@meta.data
  
  
  # ==========================================================
  # 2. Define high-sXBP1 cutoff
  # ==========================================================
  
  sx_cutoff <- quantile(
    meta$norm.sXBP1,
    probs = quantile_cutoff,
    na.rm = TRUE,
    names = FALSE
  )
  
  
  # ==========================================================
  # 3. Keep cells ABOVE cutoff
  # ==========================================================
  
  cells_high <- rownames(meta)[
    meta$norm.sXBP1 > sx_cutoff
  ]
  
  meta_high <- meta[
    cells_high,
    ,
    drop = FALSE
  ]
  
  
  # ==========================================================
  # 4. Prepare factors
  # ==========================================================
  
  # ---- Sex ----
  
  meta_high$sex_model <- factor(
    meta_high$sex
  )
  
  # Female is reference
  meta_high$sex_model <- relevel(
    meta_high$sex_model,
    ref = "female"
  )
  
  
  # ---- Donor ----
  
  meta_high$donor_model <- factor(
    meta_high$donor
  )
  
  
  # ---- Condition ----
  
  meta_high$condition_model <- factor(
    meta_high$condition
  )
  
  
  # 0h_Tg the reference condition
  
  meta_high$condition_model <- relevel(
    meta_high$condition_model,
    ref = "0h_Tg"
  )
  
  # ==========================================================
  # 5. Get normalized RNA expression
  # ==========================================================
  
  expr <- LayerData(
    seu,
    assay = assay,
    layer = layer
  )
  
  expr <- expr[
    ,
    rownames(meta_high),
    drop = FALSE
  ]
  
  # ==========================================================
  # 6. Gene filtering
  # ==========================================================
  
  pct_expressed <- Matrix::rowMeans(
    expr > 0
  )
  
  genes_keep <- names(pct_expressed)[
    pct_expressed >= min_pct
  ]
  
  
  # Optional specific gene list
  if (!is.null(genes)) {
    
    genes_keep <- intersect(
      genes_keep,
      genes
    )
  }
  
  
  expr <- expr[
    genes_keep,
    ,
    drop = FALSE
  ]
  
  
  # ==========================================================
  # 7. FIT ONE GENE
  # ==========================================================
  
  fit_one_gene <- function(i) {
    
    gene <- genes_keep[i]
    
    
    # --------------------------------------------------------
    # Expression vector
    # --------------------------------------------------------
    
    y <- as.numeric(
      expr[
        gene,
        ,
        drop = TRUE
      ]
    )
    
    
    # --------------------------------------------------------
    # Skip invariant genes
    # --------------------------------------------------------
    
    if (
      length(unique(y)) < 2 ||
      is.na(sd(y)) ||
      sd(y) == 0
    ) {
      
      return(
        data.frame(
          
          gene = gene,
          
          female_slope = NA_real_,
          female_SE = NA_real_,
          female_p = NA_real_,
          
          male_slope = NA_real_,
          male_SE = NA_real_,
          male_p = NA_real_,
          
          male_minus_female = NA_real_,
          interaction_SE = NA_real_,
          interaction_p = NA_real_,
          
          pct_expressed =
            as.numeric(pct_expressed[gene]),
          
          singular = NA,
          convergence_warning = NA,
          model_failed = TRUE
        )
      )
    }
    
    
    # --------------------------------------------------------
    # Model dataframe
    # --------------------------------------------------------
    
    model_df <- meta_high
    
    # IMPORTANT:
    # original Seurat log-normalized expression
    model_df$gene_expression <- y
    
    
    # --------------------------------------------------------
    # Fit mixed model
    # --------------------------------------------------------
    
    fit <- tryCatch(
      
      lmerTest::lmer(
        
        gene_expression ~
          
          norm.sXBP1 * sex_model +
          
          condition_model * sex_model +
          
          (1 | donor_model),
        
        data = model_df,
        
        REML = TRUE,
        
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(
            maxfun = 100000
          )
        )
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    
    # --------------------------------------------------------
    # Model failure
    # --------------------------------------------------------
    
    if (is.null(fit)) {
      
      return(
        data.frame(
          
          gene = gene,
          
          female_slope = NA_real_,
          female_SE = NA_real_,
          female_p = NA_real_,
          
          male_slope = NA_real_,
          male_SE = NA_real_,
          male_p = NA_real_,
          
          male_minus_female = NA_real_,
          interaction_SE = NA_real_,
          interaction_p = NA_real_,
          
          pct_expressed =
            as.numeric(pct_expressed[gene]),
          
          singular = NA,
          convergence_warning = TRUE,
          model_failed = TRUE
        )
      )
    }
    
    
    # ========================================================
    # 8. Fixed effects
    # ========================================================
    
    beta <- lme4::fixef(fit)
    
    beta_names <- names(beta)
    
    
    # --------------------------------------------------------
    # Female slope
    #
    # Because female is reference:
    #
    # female slope = coefficient for norm.sXBP1
    # --------------------------------------------------------
    
    female_slope <- unname(
      beta["norm.sXBP1"]
    )
    
    
    # --------------------------------------------------------
    # Find sex interaction term
    # --------------------------------------------------------
    
    interaction_name <- grep(
      "norm.sXBP1:sex_modelmale|sex_modelmale:norm.sXBP1",
      beta_names,
      value = TRUE
    )
    
    
    if (length(interaction_name) != 1) {
      
      return(
        data.frame(
          
          gene = gene,
          
          female_slope =
            female_slope,
          
          female_SE = NA_real_,
          female_p = NA_real_,
          
          male_slope = NA_real_,
          male_SE = NA_real_,
          male_p = NA_real_,
          
          male_minus_female = NA_real_,
          interaction_SE = NA_real_,
          interaction_p = NA_real_,
          
          pct_expressed =
            as.numeric(pct_expressed[gene]),
          
          singular =
            lme4::isSingular(
              fit,
              tol = 1e-4
            ),
          
          convergence_warning = TRUE,
          model_failed = TRUE
        )
      )
    }
    
    
    # --------------------------------------------------------
    # Male - female difference
    # --------------------------------------------------------
    
    delta <- unname(
      beta[
        interaction_name
      ]
    )
    
    
    # --------------------------------------------------------
    # Male slope
    #
    # male = female + interaction
    # --------------------------------------------------------
    
    male_slope <-
      female_slope +
      delta
    
    
    # ========================================================
    # 9. FEMALE SLOPE TEST
    # ========================================================
    
    L_female <- rep(
      0,
      length(beta)
    )
    
    names(L_female) <- beta_names
    
    L_female[
      "norm.sXBP1"
    ] <- 1
    
    
    female_test <- tryCatch(
      
      lmerTest::contest1D(
        fit,
        L = L_female,
        ddf = "Satterthwaite"
      ),
      
      error = function(e) NULL
    )
    
    
    # ========================================================
    # 10. MALE SLOPE TEST
    #
    # male slope =
    # norm.sXBP1 +
    # norm.sXBP1:sexmale
    # ========================================================
    
    L_male <- rep(
      0,
      length(beta)
    )
    
    names(L_male) <- beta_names
    
    L_male[
      "norm.sXBP1"
    ] <- 1
    
    L_male[
      interaction_name
    ] <- 1
    
    
    male_test <- tryCatch(
      
      lmerTest::contest1D(
        fit,
        L = L_male,
        ddf = "Satterthwaite"
      ),
      
      error = function(e) NULL
    )
    
    
    # ========================================================
    # 11. SEX DIFFERENCE TEST
    #
    # H0:
    # male slope - female slope = 0
    # ========================================================
    
    L_difference <- rep(
      0,
      length(beta)
    )
    
    names(L_difference) <- beta_names
    
    L_difference[
      interaction_name
    ] <- 1
    
    
    difference_test <- tryCatch(
      
      lmerTest::contest1D(
        fit,
        L = L_difference,
        ddf = "Satterthwaite"
      ),
      
      error = function(e) NULL
    )
    
    
    # ========================================================
    # 12. Extract statistics
    # ========================================================
    
    female_SE <- if (!is.null(female_test)) {
      female_test[["Std. Error"]][1]
    } else {
      NA_real_
    }
    
    
    female_p <- if (!is.null(female_test)) {
      female_test[["Pr(>|t|)"]][1]
    } else {
      NA_real_
    }
    
    
    male_SE <- if (!is.null(male_test)) {
      male_test[["Std. Error"]][1]
    } else {
      NA_real_
    }
    
    
    male_p <- if (!is.null(male_test)) {
      male_test[["Pr(>|t|)"]][1]
    } else {
      NA_real_
    }
    
    
    interaction_SE <- if (!is.null(difference_test)) {
      difference_test[["Std. Error"]][1]
    } else {
      NA_real_
    }
    
    
    interaction_p <- if (!is.null(difference_test)) {
      difference_test[["Pr(>|t|)"]][1]
    } else {
      NA_real_
    }
    
    
    # ========================================================
    # 13. Diagnostics
    # ========================================================
    
    singular <- lme4::isSingular(
      fit,
      tol = 1e-4
    )
    
    
    convergence_warning <-
      !is.null(
        fit@optinfo$conv$lme4$messages
      )
    
    
    # ========================================================
    # 14. Return
    # ========================================================
    
    data.frame(
      
      gene = gene,
      
      female_slope =
        female_slope,
      
      female_SE =
        female_SE,
      
      female_p =
        female_p,
      
      male_slope =
        male_slope,
      
      male_SE =
        male_SE,
      
      male_p =
        male_p,
      
      male_minus_female =
        delta,
      
      interaction_SE =
        interaction_SE,
      
      interaction_p =
        interaction_p,
      
      pct_expressed =
        as.numeric(
          pct_expressed[gene]
        ),
      
      singular =
        singular,
      
      convergence_warning =
        convergence_warning,
      
      model_failed =
        FALSE
    )
  }
  
  
  # ==========================================================
  # 15. RUN ALL GENES
  # ==========================================================
  
  gene_indices <- seq_along(
    genes_keep
  )
  
  
  # ----------------------------------------------------------
  # Parallel
  # ----------------------------------------------------------
  
  if (n_cores > 1) {
    
    result_list <- parallel::mclapply(
      
      gene_indices,
      
      fit_one_gene,
      
      mc.cores = n_cores
    )
    
    
  } else {
    
    # --------------------------------------------------------
    # Serial
    # --------------------------------------------------------
    
    result_list <- vector(
      "list",
      length(gene_indices)
    )
    
    
    for (i in gene_indices) {
      
      if (
        i == 1 ||
        i %% verbose_every == 0
      ) {
        
        cat(
          "Processing gene",
          i,
          "of",
          length(gene_indices),
          "\n"
        )
      }
      
      
      result_list[[i]] <-
        fit_one_gene(i)
    }
  }
  
  
  # ==========================================================
  # 16. Combine
  # ==========================================================
  
  results <- bind_rows(
    result_list
  )
  
  
  # ==========================================================
  # 19. FDR correction
  # ==========================================================
  
  results <- results |>
    mutate(
      
      female_FDR =
        p.adjust(
          female_p,
          method = "BH"
        ),
      
      male_FDR =
        p.adjust(
          male_p,
          method = "BH"
        ),
      
      interaction_FDR =
        p.adjust(
          interaction_p,
          method = "BH"
        ),
      
      quantile_cutoff =
        quantile_cutoff,
      
      sXBP1_cutoff =
        sx_cutoff,
      
      n_cells =
        nrow(meta_high),
      
      n_donors =
        nlevels(
          meta_high$donor_model
        )
    )
  
  
  return(results)
}		


# run transcriptome-wide


result_q50 <- run_high_sXBP1_mixed(
  seu = Beta,
  quantile_cutoff = 0.50,
  min_pct = 0.05,
  n_cores = 8
)



## plotting, scatter plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggnewscale)

# ============================================================
# 1. Prepare plotting data
# ============================================================

plot_df <- result_q50 |>
  filter(
    !is.na(gene),
    !is.na(female_slope),
    !is.na(male_slope),
    !is.na(male_minus_female),
    !is.na(interaction_FDR)
  ) |>
  mutate(
    direction = case_when(
      interaction_FDR < 0.05 & male_minus_female < 0 ~ "Female",
      interaction_FDR < 0.05 & male_minus_female > 0 ~ "Male",
      TRUE ~ "Not significant"
    )
  )


# ============================================================
# 2. Significant genes only for top-difference labels
# ============================================================

sig_df <- plot_df |>
  filter(direction != "Not significant")


# ============================================================
# 3. Top 10 genes with Female > Male slope difference
# ============================================================

top_female_diff <- sig_df |>
  filter(direction == "Female") |>
  arrange(male_minus_female) |>
  slice_head(n = 10)


# ============================================================
# 4. Top 10 genes with Male > Female slope difference
# ============================================================

top_male_diff <- sig_df |>
  filter(direction == "Male") |>
  arrange(desc(male_minus_female)) |>
  slice_head(n = 10)


# ============================================================
# 5. Gene with strongest POSITIVE slope overall
#
# Selects the gene with the largest positive slope
# in either females or males
# ============================================================

top_positive_gene <- plot_df |>
  mutate(
    max_positive_slope = pmax(
      female_slope,
      male_slope,
      na.rm = TRUE
    )
  ) |>
  arrange(desc(max_positive_slope)) |>
  slice(1)

top_positive_gene


# ============================================================
# 6. Combine colored labels
#
# Remove top_positive_gene so it is not labeled twice
# ============================================================

label_genes <- bind_rows(
  top_female_diff,
  top_male_diff
) |>
  distinct(gene, .keep_all = TRUE) |>
  filter(gene != top_positive_gene$gene)


# ============================================================
# 7. Set identical x/y limits
# ============================================================

lims <- range(
  c(
    plot_df$female_slope,
    plot_df$male_slope
  ),
  na.rm = TRUE
)

# add 5% padding
pad <- diff(lims) * 0.05

lims <- c(
  lims[1] - pad,
  lims[2] + pad
)


# ============================================================
# 8. Plot
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = female_slope,
    y = male_slope
  )
) +
  
  # ----------------------------------------------------------
# All genes
# ----------------------------------------------------------
geom_point(
  aes(color = direction),
  size = 1.6,
  alpha = 0.75
) +
  
  scale_color_manual(
    values = c(
      "Female" = "#7FC97F",
      "Male" = "#BEAED4",
      "Not significant" = "grey80"
    ),
    name = NULL
  ) +
  
  # ----------------------------------------------------------
# Bigger dots for top 10 sex-difference genes
# ----------------------------------------------------------
geom_point(
  data = label_genes,
  aes(
    x = female_slope,
    y = male_slope,
    color = direction
  ),
  size = 3.2,
  alpha = 1
) +
  
  # ----------------------------------------------------------
# Strongest positive gene:
# BLACK dot
# ----------------------------------------------------------
geom_point(
  data = top_positive_gene,
  aes(
    x = female_slope,
    y = male_slope
  ),
  color = "black",
  size = 4,
  alpha = 1
) +
  
  # ----------------------------------------------------------
# Reference lines
# ----------------------------------------------------------
geom_abline(
  slope = 1,
  intercept = 0,
  linetype = "dashed",
  linewidth = 0.5
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  
  # ==========================================================
# New color scale for text labels
# ==========================================================
ggnewscale::new_scale_color() +
  
  # ----------------------------------------------------------
# Colored labels for top 10 in each direction
# ----------------------------------------------------------
geom_text_repel(
  data = label_genes,
  aes(
    x = female_slope,
    y = male_slope,
    label = gene,
    color = direction
  ),
  size = 6,
  fontface = "bold",
  max.overlaps = Inf,
  box.padding = 0.5,
  point.padding = 0.4,
  min.segment.length = 0,
  show.legend = FALSE
) +
  
  scale_color_manual(
    values = c(
      "Female" = "#1B9E77",
      "Male" = "#6A3D9A"
    )
  ) +
  
  # ----------------------------------------------------------
# BLACK label for strongest positive gene
# ----------------------------------------------------------
geom_text_repel(
  data = top_positive_gene,
  aes(
    x = female_slope,
    y = male_slope,
    label = gene
  ),
  color = "black",
  size = 6,
  fontface = "bold",
  box.padding = 0.7,
  point.padding = 0.5,
  min.segment.length = 0,
  max.overlaps = Inf,
  show.legend = FALSE
) +
  
  # ----------------------------------------------------------
# Same X/Y limits and equal aspect ratio
# ----------------------------------------------------------
scale_x_continuous(
  limits = lims
) +
  
  scale_y_continuous(
    limits = lims
  ) +
  
  coord_equal() +

  
  labs(
    x = "Female slope",
    y = "Male slope"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 20,face="bold"),
    axis.text = element_text(size = 20, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    legend.position = "top",
    legend.justification = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.margin = margin(6, 6, 6, 6)
  )
p

		  


