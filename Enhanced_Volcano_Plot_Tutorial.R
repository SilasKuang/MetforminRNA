# This tutorial is adapted from:
# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

if(!require(EnhancedVolcano))
  BiocManager::install('EnhancedVolcano')
if(!require(airway))
  BiocManager::install("airway")
if(!require(magrittr))
  BiocManager::install("magrittr")
if(!require(org.Hs.eg.db))
  BiocManager::install("org.Hs.eg.db")
if(!require(DESeq2))
  BiocManager::install("DESeq2")

library(EnhancedVolcano)

library(airway)
library(magrittr)

data('airway')
airway$dex %<>% relevel('untrt')

ens <- rownames(airway)

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway <- airway[keep,]

library(DESeq2)

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds,
               contrast = c('dex','trt','untrt'))
res <- lfcShrink(dds,
                 contrast = c('dex','trt','untrt'), res=res, type = 'normal')

# Basic plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

# Modify cut-offs for log2FC and P value; specify title; adjust point and label size 
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

# Modify cut-offs for log2FC and P value; specify title; adjust point and label size
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

# Adjust shape of plotted points
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                pointSize = 4.0,
                labSize = 6.0,
                shape = 8,
                colAlpha = 1)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                shape = c(1, 4, 23, 25),
                colAlpha = 1)

# Adjust cut-off lines and add extra threshold lines
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-6, 6),
                title = 'N061011 versus N61311',
                pCutoff = 10e-12,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 6.0,
                colAlpha = 1,
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.8,
                hline = c(10e-20,
                          10e-20 * 10e-30,
                          10e-20 * 10e-60,
                          10e-20 * 10e-90),
                hlineCol = c('pink', 'hotpink', 'purple', 'black'),
                hlineType = c('solid', 'longdash', 'dotdash', 'dotted'),
                hlineWidth = c(1.0, 1.5, 2.0, 2.5),
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

# Adjust legend position, size, and text
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-12,
                FCcutoff = 1.5,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)

# Fit more labels by adding connectors
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-32,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)

# Only label key variables
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('TMEM176B','ADH1A'),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                shape = c(4, 35, 17, 18),
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5.0)

# Draw labels in boxes
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('VCAM1','KCTD12','ADAM12',
                              'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA'),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

# Italicise labels and flip volcano on it’s side
lab_italics <- paste0("italic('", rownames(res), "')")
selectLab_italics = paste0(
  "italic('",
  c('VCAM1','KCTD12','ADAM12', 'CXCL12','CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA'),
  "')")

EnhancedVolcano(res,
                lab = lab_italics,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = selectLab_italics,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('black', 'pink', 'purple', 'red3'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()

# Over-ride colouring scheme with custom key-value pairs
# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  res$log2FoldChange < -2.5, 'royalblue',
  ifelse(res$log2FoldChange > 2.5, 'gold',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Custom colour over-ride',
                pCutoff = 10e-14,
                FCcutoff = 1.0,
                pointSize = 3.5,
                labSize = 4.5,
                shape = c(6, 4, 2, 11),
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'left',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                arrowheads = FALSE,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black')

# Over-ride colour and/or shape scheme with custom key-value pairs
# define different cell-types that will be shaded
celltype1 <- c('VCAM1','KCTD12','ADAM12','CXCL12')
celltype2 <- c('CACNB2','SPARCL1','DUSP1','SAMHD1','MAOA')

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  rownames(res) %in% celltype1, 17,
  ifelse(rownames(res) %in% celltype2, 64,
         3))
keyvals.shape[is.na(keyvals.shape)] <- 3
names(keyvals.shape)[keyvals.shape == 3] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Cell-type 1'
names(keyvals.shape)[keyvals.shape == 64] <- 'Cell-type 2'

p1 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      selectLab = rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                      xlab = bquote(~Log[2]~ 'fold change'),
                      title = 'Custom shape over-ride',
                      pCutoff = 10e-14,
                      FCcutoff = 1.0,
                      pointSize = 4.5,
                      labSize = 4.5,
                      shapeCustom = keyvals.shape,
                      colCustom = NULL,
                      colAlpha = 1,
                      legendLabSize = 15,
                      legendPosition = 'left',
                      legendIconSize = 5.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      colConnectors = 'grey50',
                      gridlines.major = TRUE,
                      gridlines.minor = FALSE,
                      border = 'partial',
                      borderWidth = 1.5,
                      borderColour = 'black')

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals.colour <- ifelse(
  res$log2FoldChange < -2.5, 'royalblue',
  ifelse(res$log2FoldChange > 2.5, 'gold',
         'black'))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'gold'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'low'

p2 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      selectLab = rownames(res)[which(names(keyvals) %in% c('High', 'Low'))],
                      xlab = bquote(~Log[2]~ 'fold change'),
                      title = 'Custom shape & colour over-ride',
                      pCutoff = 10e-14,
                      FCcutoff = 1.0,
                      pointSize = 5.5,
                      labSize = 0.0,
                      shapeCustom = keyvals.shape,
                      colCustom = keyvals.colour,
                      colAlpha = 1,
                      legendPosition = 'right',
                      legendLabSize = 15,
                      legendIconSize = 5.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      colConnectors = 'grey50',
                      gridlines.major = TRUE,
                      gridlines.minor = FALSE,
                      border = 'full',
                      borderWidth = 1.0,
                      borderColour = 'black')

library(gridExtra)
library(grid)
grid.arrange(p1, p2,
             ncol=2,
             top = textGrob('EnhancedVolcano',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))

# Encircle / highlight certain variables
# define different cell-types that will be shaded
celltype1 <- c('VCAM1','CXCL12')
celltype2 <- c('SORT1', 'KLF15')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c(celltype1, celltype2),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'Shading cell-type 1|2',
                pCutoff = 10e-14,
                FCcutoff = 1.0,
                pointSize = 8.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                shape = 42,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 20,
                legendIconSize = 20.0,
                # encircle
                encircle = celltype1,
                encircleCol = 'black',
                encircleSize = 2.5,
                encircleFill = 'pink',
                encircleAlpha = 1/2,
                # shade
                shade = celltype2,
                shadeAlpha = 1/2,
                shadeFill = 'skyblue',
                shadeSize = 1,
                shadeBins = 5,
                drawConnectors = TRUE,
                widthConnectors = 2.0,
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 5,
                borderColour = 'black')

# Highlighting key variables via custom point sizes
library("pasilla")
pasCts <- system.file("extdata", "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
rownames(coldata) <- sub("fb", "", rownames(coldata))
cts <- cts[, rownames(coldata)]
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
dds <- DESeq(dds)
res <- results(dds)

p1 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = "log2FoldChange",
                      y = "pvalue",
                      pCutoff = 10e-4,
                      FCcutoff = 2,
                      ylim = c(0, -log10(10e-12)),
                      pointSize = c(ifelse(res$log2FoldChange>2, 8, 1)),
                      labSize = 6.0,
                      shape = c(6, 6, 19, 16),
                      title = "DESeq2 results",
                      subtitle = "Differential expression",
                      caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      col = c("grey30", "forestgreen", "royalblue", "red2"),
                      colAlpha = 0.9,
                      drawConnectors = TRUE,
                      hline = c(10e-8),
                      widthConnectors = 0.5)

p1

# Change to continuous colour scheme
p1 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = "log2FoldChange",
                      y = "pvalue",
                      pCutoff = 10e-4,
                      FCcutoff = 2,
                      ylim = c(0, -log10(10e-12)),
                      pointSize = c(ifelse(res$log2FoldChange>2, 8, 1)),
                      labSize = 6.0,
                      shape = c(6, 6, 19, 16),
                      title = "DESeq2 results",
                      subtitle = "Differential expression",
                      caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                      legendPosition = "right",
                      legendLabSize = 14,
                      colAlpha = 0.9,
                      colGradient = c('red3', 'royalblue'),
                      drawConnectors = TRUE,
                      hline = c(10e-8),
                      widthConnectors = 0.5)

p1

# Custom axis tick marks
p1 +
  ggplot2::coord_cartesian(xlim=c(-6, 6)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1))


