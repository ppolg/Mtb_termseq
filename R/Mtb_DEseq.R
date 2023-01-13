###########################################################################

# Additional data analysis for D'Halluin et al., 2022
#   • Note - most of this is copy-paste and not solidly defined functions,
#   • due to lack of time (and energy). Apologies in advance if you check it!

#   • DEseq initial analysis on RhoDUC RNAseq
#   • DEseq analysis of term-seq (to select actually relevant Rho-dependents?)

###########################################################################

########
# Libs #
########

packages <- c("hash","here","tidyverse","ape","data.table","ggthemes","ggplot2","ggridges","Biostrings","forcats","reshape2","DESeq2","vsn","RColorBrewer","pheatmap")
invisible(lapply(packages, require, character.only = TRUE))

########
# Func #
########

# nin - because it just looks clean and neat!
`%nin%` = Negate(`%in%`)

########
# Init #
########

colours <- hash::hash(
  red = "#A3280A",
  orange = "#E3812B",
  brown = "#8A5122",
  yellow = "#E0D253",
  grey = "#858482",
  green = "#195928",
  lime = "#83e649",
  lightlime = "#c5f081",
  lightlime2 = "#bed698",
  blue = "#4A749E",
  darkblue = "#18099c",
  azure = "#092780",
  purple = "#612882",
  pink = "#d64daf",
  darkgold = "#755807",
  alex_R1 = "#F0BFC1",
  alex_R2 = "#D46859",
  alex_R3 = "#871E12",
  alex_B1 = "#D7E0F7",
  alex_B2 = "#87A9D9",
  alex_B3 = "#5577AA",
  alex_B1_edge = "#BBD0FA"
)

# ggplot2 theme for figures
theme_alex <- function(base_size=10) {
  (theme_foundation(base_size=base_size, base_family="Arial")
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(2.8), hjust = 0.5),
           text = element_text(color = "black"),
           panel.background = element_rect(colour = NA, fill = "white"),
           plot.background = element_rect(colour = NA, fill = "white"),
           panel.border = element_blank(),
           axis.title = element_text(face = "bold",size = rel(1.4)),
           axis.title.y = element_text(angle=90,vjust = 3, size = rel(1)),
           axis.title.x = element_text(vjust = -0.1, size = rel(1)),
           axis.text = element_text(face="bold",size = rel(1.3)), 
           axis.line = element_line(colour="black",size=1),
           axis.ticks = element_line(colour="black",size = 1),
           axis.ticks.length = unit(.35, "cm"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(color = NA, fill = "white"),
           legend.position = "top",
           legend.background= element_rect(color = NA, fill = "white"),
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.spacing.x = unit(0.3, "cm"),
           legend.text = element_text(color = "black", size = rel(1.2)),
           legend.title = element_text(face="italic", color = "black"),
           plot.margin=unit(c(12,6,6,6),"mm"),
           strip.background=element_rect(colour="grey90",fill="grey70"),
           strip.text = element_text(face="bold")
   ))
}

########
# Main #
########

# Define directory to get count files from
directory <- paste(here::here("NUGA_Data"))

# Load data

sample_files <- grep(".count",list.files(directory),value=TRUE)
sample_times <- str_extract(sample_files, "[[:digit:]]+h*")
sample_types <- str_extract(sample_files, "rna|term")
sample_replicates <-  str_extract(sample_files, "r[123]")

sample_table <- data.frame(sampleName = sample_files,
                       fileName = sample_files,
                       time = sample_times,
                       type = sample_types,
                       rep = sample_replicates)

rna_table <- sample_table %>% filter(type == "rna")
term_table <- sample_table %>% filter(type == "term")

no6 <- rna_table %>%
  filter(time != '6h')

dds <- DESeqDataSetFromHTSeqCount(sampleTable = rna_table,
                                      directory = directory,
                                      design = formula(~ time))

dds_term <- DESeqDataSetFromHTSeqCount(sampleTable = term_table,
                                  directory = directory,
                                  design = formula(~ time))


dds_no6 <- DESeqDataSetFromHTSeqCount(sampleTable = no6,
                                  directory = directory,
                                  design = formula(~ time))

# Pre-filter
filter_limit <- 10
keep <- rowSums(counts(dds)) >= filter_limit
dds <- dds[keep,]

keep_term <- rowSums(counts(dds_term)) >= filter_limit
dds_term <- dds_term[keep_term,]

keep_no6 <- rowSums(counts(dds_no6)) >= filter_limit
dds_no6 <- dds_no6[keep,]

# Differential expression
dds <- DESeq(dds)
dds_term <- DESeq(dds_term)
dds_no6 <- DESeq(dds_no6)

# MA

cap_limit <- 4

# 0h to 3h
ma3 <- as.data.frame(results(dds, contrast = c("time","0h","3h"), cooksCutoff = F))

ma3 <- ma3 %>%
  dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit, 
                                             sign(log2FoldChange) * cap_limit, log2FoldChange))

ggplot(ma3) +
  geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
             size = 0.8,
             color = ifelse(ma3$padj < 0.05, colours$alex_R2, "grey60"),
             fill = ifelse(ma3$padj < 0.05, colours$alex_R2, "grey60"),
             shape = ifelse(ma3$log2FoldChange >= cap_limit, 24, 
                            ifelse(ma3$log2FoldChange <= -1*cap_limit, 25, 21)),
             alpha = 0.8) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6")) +
  scale_y_continuous(limits=c(-4,4)) +
  theme_alex() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
  xlab("Mean of counts") +
  ylab("log2-fold change") +
  ggtitle("MA plot between 0h and 3h")

# 0h to 4.5h
ma45 <- as.data.frame(results(dds, contrast = c("time","0h","45h"), cooksCutoff = F))

ma45 <- ma45 %>%
  dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit, 
                               sign(log2FoldChange) * cap_limit, log2FoldChange))

ggplot(ma45) +
  geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
             size = 0.8,
             color = ifelse(ma45$padj < 0.05, colours$alex_R2, "grey60"),
             fill = ifelse(ma45$padj < 0.05, colours$alex_R2, "grey60"),
             shape = ifelse(ma45$log2FoldChange >= cap_limit, 24, 
                            ifelse(ma45$log2FoldChange <= -1*cap_limit, 25, 21)),
             alpha = 0.8) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6")) +
  scale_y_continuous(limits=c(-4,4)) +
  theme_alex() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
  xlab("Mean of counts") +
  ylab("log2-fold change") +
  ggtitle("MA plot between 0h and 4.5h")

# 0h to 6h
ma6 <- as.data.frame(results(dds, contrast = c("time","0h","6h"), cooksCutoff = F))

ma6 <- ma6 %>%
  dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit, 
                                             sign(log2FoldChange) * cap_limit, log2FoldChange))

ggplot(ma6) +
  geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
             size = 0.8,
             color = ifelse(ma6$padj < 0.05, colours$alex_R2, "grey60"),
             fill = ifelse(ma6$padj < 0.05, colours$alex_R2, "grey60"),
             shape = ifelse(ma6$log2FoldChange >= cap_limit, 24, 
                            ifelse(ma6$log2FoldChange <= -1*cap_limit, 25, 21)),
             alpha = 0.8) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6")) +
  scale_y_continuous(limits=c(-4,4)) +
  theme_alex() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
  xlab("Mean of counts") +
  ylab("log2-fold change") +
  ggtitle("MA plot between 0h and 6h")

# MA - term
# 0h to 3h
ma3_term <- as.data.frame(results(dds_term, contrast = c("time","0h","3h"), cooksCutoff = F))

ma3_term <- ma3_term %>%
  dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit, 
                                             sign(log2FoldChange) * cap_limit, log2FoldChange))

ggplot(ma3_term) +
  geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
             size = 0.8,
             color = ifelse(ma3_term$padj < 0.05, colours$alex_R2, "grey60"),
             fill = ifelse(ma3_term$padj < 0.05, colours$alex_R2, "grey60"),
             shape = ifelse(ma3_term$log2FoldChange >= cap_limit, 24, 
                            ifelse(ma3_term$log2FoldChange <= -1*cap_limit, 25, 21)),
             alpha = 0.8) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6")) +
  scale_y_continuous(limits=c(-4,4)) +
  theme_alex() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
  xlab("Mean of counts") +
  ylab("log2-fold change") +
  ggtitle("MA plot between 0h and 3h (termseq)")

# 0h to 4.5h
ma45_term <- as.data.frame(results(dds_term, contrast = c("time","0h","45h"), cooksCutoff = F))

ma45_term <- ma45_term %>%
  dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit, 
                                             sign(log2FoldChange) * cap_limit, log2FoldChange))

ggplot(ma45_term) +
  geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
             size = 0.8,
             color = ifelse(ma45_term$padj < 0.05, colours$alex_R2, "grey60"),
             fill = ifelse(ma45_term$padj < 0.05, colours$alex_R2, "grey60"),
             shape = ifelse(ma45_term$log2FoldChange >= cap_limit, 24, 
                            ifelse(ma45_term$log2FoldChange <= -1*cap_limit, 25, 21)),
             alpha = 0.8) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6")) +
  scale_y_continuous(limits=c(-4,4)) +
  theme_alex() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
  xlab("Mean of counts") +
  ylab("log2-fold change") +
  ggtitle("MA plot between 0h and 4.5h (termseq)")

# 0h to 6h
ma6_term <- as.data.frame(results(dds_term, contrast = c("time","0h","6h"), cooksCutoff = F))

ma6_term <- ma6_term %>%
  dplyr::mutate(log2FoldChangeCapped= ifelse(abs(log2FoldChange) > cap_limit, 
                                             sign(log2FoldChange) * cap_limit, log2FoldChange))

ggplot(ma6_term) +
  geom_point(aes(x = baseMean, y = log2FoldChangeCapped),
             size = 0.8,
             color = ifelse(ma6_term$padj < 0.05, colours$alex_R2, "grey60"),
             fill = ifelse(ma6_term$padj < 0.05, colours$alex_R2, "grey60"),
             shape = ifelse(ma6_term$log2FoldChange >= cap_limit, 24, 
                            ifelse(ma6_term$log2FoldChange <= -1*cap_limit, 25, 21)),
             alpha = 0.8) +
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000,1000000),labels = c("1e0","1e1","1e2","1e3","1e4","1e5","1e6")) +
  scale_y_continuous(limits=c(-4,4)) +
  theme_alex() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.8) +
  xlab("Mean of counts") +
  ylab("log2-fold change") +
  ggtitle("MA plot between 0h and 6h (termseq)")

# heatmap
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste0(vsd$time,"-",vsd$rep)
colnames(sampleDistMatrix) <- NULL
colorMap <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colorMap)

vsd_term <- vst(dds_term, blind=FALSE)
sampleDists <- dist(t(assay(vsd_term)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste0(vsd_term$time,"-",vsd_term$rep)
colnames(sampleDistMatrix) <- NULL
colorMap <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colorMap)

vsd_no6 <- vst(dds_no6, blind=FALSE)
sampleDists <- dist(t(assay(vsd_no6)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste0(vsd_no6$time,"-",vsd_no6$rep)
colnames(sampleDistMatrix) <- NULL
colorMap <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colorMap)

# PCA
pca <- plotPCA(vsd, intgroup="time", returnData = T)
percentVar <- round(100 * attr(pca, "percentVar"))

ggplot(pca) +
  geom_point(aes(PC1,PC2, colour = time),
             size = 3) +
  theme_alex() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c(colours$darkgold,colours$azure,
                                 colours$purple,colours$green), name=NULL) 


pca_term <- plotPCA(vsd_term, intgroup="time", returnData = T)
percentVar <- round(100 * attr(pca_term, "percentVar"))

ggplot(pca_term) +
  geom_point(aes(PC1,PC2, colour = time),
             size = 3) +
  theme_alex() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c(colours$darkgold,colours$azure,
                                 colours$purple,colours$green), name=NULL) +
  ggtitle("PCA for termseq")

pca_no6 <- plotPCA(vsd_no6, intgroup="time", returnData = T)
percentVar <- round(100 * attr(pca_no6, "percentVar"))

ggplot(pca_no6) +
  geom_point(aes(PC1,PC2, colour = time),
             size = 3) +
  theme_alex() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c(colours$darkgold,colours$azure,
                                 colours$purple,colours$green), name=NULL) +
  ggtitle("PCA without 6h")



