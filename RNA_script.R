if(!require(tidyverse))
  install.packages("tidyverse")
if(!require(ggplot2))
  install.packages("ggplot2")
if(!require(clusterProfiler))
  BiocManager::install("clusterProfiler")
if(!require(org.Mm.eg.db))
  BiocManager::install("org.Mm.eg.db")
if(!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

library(readxl)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(EnhancedVolcano)

# Tutorial for Enhanced Volcano Plot:
# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

# Note: Since in this project, we have not tested for MET only
# So, MET means AET+MET unless otherwise specified
import_dataset <- function(filename){
  rawDGE <- read_excel(filename)
  curatedDGE <- rawDGE[,c(8,2,5,6,7)]
  colnames(curatedDGE) <- c("Symbol", "logFC", "pval", "FDR", "Group")
  
  # Add a column named direction, to show whether this gene is UP-regulated
  # or DOWN-regulated after the treatment
  curatedDGE <- curatedDGE %>% 
    mutate("logFDR" = -log(FDR)) %>%
    mutate(direction = case_when(FDR < 0.05 & logFC > 0 ~ "UP",
                                 FDR < 0.05 & logFC < 0 ~ "DOWN",
                                 FDR >= 0.05 ~ "NS")) 
  return(curatedDGE)
}

curatedDGE_AET_vs_MET <- import_dataset("Genes_AET_vs_AET+MET.xlsx")
curatedDGE_MET_vs_SED <- import_dataset("Genes_AET+MET_vs_SED.xlsx")
curatedDGE_SED_vs_AET <- import_dataset("Genes_SED_vs_AET.xlsx")

category_plot <- function(curatedDGE, title){
  totExpressed <- curatedDGE %>% 
    tally() %>% 
    mutate(direction = "all trans")
  sigDirection <- curatedDGE %>% 
    group_by(direction) %>% 
    tally() 
  allSig <- curatedDGE %>%
    filter(direction != "NS") %>%
    tally() %>% 
    mutate(direction = "normal Sig")
  restrictSig <- curatedDGE %>% 
    filter(FDR <= 0.01) %>% 
    tally() %>% 
    mutate(direction = "restrictive Sig")
  
  # Plot transcript categories
  transcriptBreakdown <- bind_rows(
    totExpressed, sigDirection, allSig, restrictSig)
  transcriptNumsPlot <- transcriptBreakdown %>% 
    arrange(match(direction, c("all trans", "DOWN", "UP", 
                               "normal Sig", "NS", 
                               "restrictive Sig")), 
            desc(direction)) %>% 
    mutate(direction = factor(direction, levels = direction)) %>% 
    ggplot(aes(x = direction, y = n, 
               fill = direction)) + geom_col(color = "black", 
                                             size = 0.25) +
    xlab("Transcript Category") + ylab("count") + theme_bw() + 
    theme(axis.line = element_line(color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") + 
    # 10000 is the number of rows in curatedDGE
    scale_y_continuous(expand = c(0,0), limits = c(0, 10000)) + 
    geom_text(aes(label = n), vjust = -0.5) + 
    scale_fill_manual(values = c("black", "#808080", 
                                 "#524fa1", "#fdb913", 
                                 "red", "cyan"))
  ggsave(paste(title,"transcriptCategories.pdf",sep=" "), 
         transcriptNumsPlot, 
         height = 4, width = 4)
}

category_plot(curatedDGE_AET_vs_MET, "AET_vs_MET")
category_plot(curatedDGE_MET_vs_SED, "MET_vs_SED")
category_plot(curatedDGE_SED_vs_AET, "SED_vs_AET")

volcano_plot <- function(curatedDGE, title){
  curatedDGE %>% 
    ggplot(aes(logFC, logFDR, color = direction,label = curatedDGE$Symbol)) + geom_text()+
    geom_point(size = 2) + theme_bw() + 
    theme(axis.line = element_line(color = "black"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("#662d91", 
                                  "#999999", "#ed1c24")) + 
    xlab(expression(paste("Log"[2],"(FC)"))) + 
    ylab(expression(paste("-log"[10],"(FDR)"))) +
    ggtitle("Volcano plot of AET+MET vs AET Control")
  ggsave(paste(title,"Volcano.pdf",sep=" "),
         volcanoPlot, height = 4, width = 4) 
}

volcano_plot(curatedDGE_AET_vs_MET, "AET_vs_MET")
volcano_plot(curatedDGE_MET_vs_SED, "MET_vs_SED")
volcano_plot(curatedDGE_SED_vs_AET, "SED_vs_AET")


sig_gene <- function(curatedDGE, title){
  # Extract restrictively significant genes
  sigGenes <- curatedDGE[,c("Symbol", "logFC","FDR")] %>% 
    filter(FDR < 0.01) %>% 
    # Calculate the non-log fold change and filter out rows with values smaller than 50%
    # log2(1.5) is about 0.5849, we take 0.58 here
    filter(abs(logFC)>= 0.58)
  
  write.csv(sigGenes, file = paste(title,"sigGenes.csv",sep=" "))
  return(sigGenes)
}

sigGenes_AET_vs_MET <- sig_gene(curatedDGE_AET_vs_MET, "AET_vs_MET")
sigGenes_MET_vs_SED <- sig_gene(curatedDGE_MET_vs_SED, "MET_vs_SED")
sigGenes_SED_vs_AET <- sig_gene(curatedDGE_SED_vs_AET, "SED_vs_AET")

# Run the enrichedKEGG GSEA Pathway analysis
path_generate <- function(geneset, geneType){
  # Map the EMSEMBL IDs to their ENTREZID
  sameGenesEntrez <- na.exclude(mapIds(org.Mm.eg.db, 
                                       keys = as.character(geneset$Symbol),
                                       keytype = "SYMBOL", column = "ENTREZID"))
  
  # Function for enrichKEGG analysis and pathway dataset generation
  # By default, this works on mouse (mmu) only
  allPaths <- enrichKEGG(gene = sameGenesEntrez, organism = "mmu")
  write.csv(allPaths, file = paste("allPaths_",geneType,".csv",sep = ""))
  
  allPaths <- as.data.frame(allPaths) %>%
    mutate(bkgdSize = 
             as.numeric(substring(BgRatio, 
                                  regexpr("/", BgRatio) + 1))) %>%
    mutate(pathBkgd = 
             as.numeric(substring(BgRatio, 1,
                                  regexpr("/", BgRatio)-1))) %>%
    mutate(bkgdPerc = pathBkgd/bkgdSize) %>%
    mutate(GeneRatTotal = 
             as.numeric(substring(GeneRatio, 
                                  regexpr("/", GeneRatio) + 1))) %>%
    mutate(percPath = Count/GeneRatTotal) %>%
    mutate(Enrichment = percPath/bkgdPerc)
  return(allPaths)
} 

stringent_path_plot <- function(pathway, title){
  pdf(paste(title,"KEGG pathways p less than 0.01.pdf",sep=" "))  
  plot(pathway %>% 
        filter(p.adjust < 0.01) %>% 
        ggplot(aes(x = Enrichment, y = Description, 
                   color = p.adjust, size = Count)) + 
        geom_point() + expand_limits(x = 0) + 
        labs(x = "Enrichment", y = "KEGG pathway", 
             color = "FDR", size = "Count") +
        theme_bw() + scale_color_gradient(low = "#B72668", 
                                          high = "#dba3b2") + 
        ggtitle(paste(title,"KEGG pathways p < 0.01", sep = " ")))
  dev.off()
}

relaxed_path_plot <- function(pathway, title){
  pdf(paste(title,"KEGG p between 0.05 and 0.01.pdf",sep=" "))  
  plot(pathway %>% 
         filter(p.adjust < 0.05 & p.adjust >= 0.01) %>% 
         ggplot(aes(x = Enrichment, y = Description, 
                    color = p.adjust, size = Count)) + 
         geom_point() + expand_limits(x = 0) + 
         labs(x = "Enrichment", y = "KEGG pathway",
              color = "FDR", size = "Count") +
         theme_bw() + scale_color_gradient(low = "#B72668", 
                                           high = "#dba3b2") + 
         ggtitle(paste(title,"KEGG pathways 0.05 < p < 0.01", sep = " ")))
  dev.off()
}

allPaths_AET_vs_MET <- path_generate(sigGenes_AET_vs_MET, "AET_vs_MET")
allPaths_MET_vs_SED <- path_generate(sigGenes_MET_vs_SED, "MET_vs_SED")
allPaths_SED_vs_AET <- path_generate(sigGenes_SED_vs_AET, "SED_vs_AET")

stringent_path_plot(allPaths_AET_vs_MET, "AET_vs_MET")
relaxed_path_plot(allPaths_AET_vs_MET, "AET_vs_MET")

stringent_path_plot(allPaths_MET_vs_SED, "MET_vs_SED")
relaxed_path_plot(allPaths_MET_vs_SED, "MET_vs_SED")

stringent_path_plot(allPaths_SED_vs_AET, "SED_vs_AET")
relaxed_path_plot(allPaths_SED_vs_AET, "SED_vs_AET")