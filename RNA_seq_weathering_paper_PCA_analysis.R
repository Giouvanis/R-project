install.packages("BiocManager")
library("BiocParallel")
library("DESeq2")
library("ggplot2")
library("magrittr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("rgl")
library("edgeR")


setwd("~/Desktop/weathering paper_22/RNA_seq")



##file import and test (2 factorial:symbiotic and mineral)
raw_counts <- as.data.frame(read.csv("raw_counts.csv",header = TRUE, row.names=1 ))

data.class(raw_counts)

samples <- read.csv("sample_list.csv", row.names=1)
samples <- samples[c("Symbiotic","Mineral")]
all(rownames(samples) %in% colnames(raw_counts))

raw_counts<- raw_counts[, rownames(samples)]
all(rownames(samples) == colnames(raw_counts))


set.seed(12345)

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = samples,
                              design= ~ Mineral+ Symbiotic + Mineral:Symbiotic)


###Warning message after run the above:
###In DESeqDataSet(se, design = design, ignoreRank) :
###some variables in design formula are characters, converting to factors

###This is a warning, not an error. In R, a warning means you can continue (the function will run) but it's telling you something that has happened.
###Here it is saying that some variable(s) in the design formula were characters, and that they have been converted to factors.

dds$Mineral <- relevel(dds$Mineral, ref="N","Y")
dds$symb <- relevel(dds$Symbiotic, ref="NS","S")

data.class(dds)

# first filtering, removing rows with read counts lower that and equal to 10

dds2 <- dds[ rowSums(counts(dds)) >= 500, ]

data.class(dds2)



#PCA_script_thesis

transformedreadcounts <- rlog(dds2, blind=TRUE)

data.class(transformedreadcounts)


groups <- c("NS_HAP_minus1",	"NS_HAP_minus2",	"NS_HAP_minus3",	"NS_HAP_minus4",	"S_HAP_minus1",	"S_HAP_minus2",	"S_HAP_minus3",	"S_HAP_minus4",	"NS_HAP_plus1",	"NS_HAP_plus2",	"NS_HAP_plus3",	"NS_HAP_plus4",	"S_HAP_plus1",	"S_HAP_plus2",	"S_HAP_plus3",	"S_HAP_plus4")

PCA <- plotPCA(transformedreadcounts, intgroup =c("Symbiotic","Mineral"), returnData=TRUE)

write.csv(PCA, "PCA.csv")

percentVar <- round(100 * attr(PCA, "percentVar"))
ggplot(PCA, aes(PC1, PC2, color=Mineral, shape=Symbiotic)) +
  geom_point(size=3) +
  geom_text(aes(label=rownames(samples)),hjust=0.5, vjust=-1)+
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() 



#PCA_test1
install.packages("tidyverse")
library(tidyverse)


dds <- DESeqDataSetFromMatrix(countData = raw_counts,
  colData = samples, design = ~1)
dds_norm <- vst(dds2)
plotPCA(dds_norm,intgroup = "Symbiotic")
plotPCA(dds_norm,intgroup = "Mineral")
plotPCA(dds_norm,intgroup = c("Symbiotic", "Mineral"))

pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("Symbiotic", "Mineral"),
    returnData = TRUE)



PCA <- plotPCA(transformedreadcounts, intgroup =c("Symbiotic","Mineral"), returnData=TRUE)
annotated_pca_plot <- ggplot(
  pca_results,
  aes(x = PC1, y = PC2, color = Mineral,shape = Symbiotic)) + geom_point(size=5)

annotated_pca_plot



###Uniform Manifold Approximation and Projection (UMAP)

#this is something i tried but i don't like the output, not informative in my case

install.packages("umap")
library(umap)

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = samples, design = ~1)
dds_norm <- vst(dds2)

normalized_counts <- assay(dds_norm)

umap_results <- umap::umap(normalized_counts)

umap_results[["layout"]]

umap_plot_df <- data.frame(umap_results$layout)

tibble::rownames_to_column(umap_plot_df)

table1 <- data.frame(umap_plot_df)


ggplot(umap_plot_df, aes(x = X1, y = X2)) +
  geom_point()


#Interpretation of UMAP plot and results
#Note that the coordinates of UMAP output for any given cell can change dramatically depending on parameters, and even run to run with the same parameters (Also why setting the seed is important). This means that you should not rely too heavily on the exact values of UMAP’s output.
#One particular limitation of UMAP is that while observed clusters have some meaning, the distance between clusters usually does not (nor does cluster density). The fact that two clusters are near each other should NOT be interpreted to mean that they are more related to each other than to more distant clusters. (There is some disagreement about whether UMAP distances have more meaning, but it is also probably safer to assume they don’t.)
#Playing with the parameters so you can fine-tune them is a good way to give you more information about a particular analysis as well as the data itself. Feel free to try playing with the parameters on your own in the code chunks above!
#In summary, a good rule of thumb to remember is: if the results of an analysis can be completely changed by changing its parameters, you should be very cautious when it comes to the conclusions you draw from it as well as having good rationale for the parameters you choose (adapted from Childhood Cancer Data Lab (2020) training materials).
### in my case i cannot see any useful info.



#this is the lfda which i am interested in...i think it will give a clearer image of my pca results...
###but code needs fixing

install.packages('lfda')

devtools::install_github('terrytangyuan/lfda')

library(lfda)

setwd("~/Desktop/weathering paper_22/RNA_seq")
raw_counts <- read.csv("raw_counts.csv",header = TRUE, row.names=1 )
raw_counts <- as.data.frame(read.csv("raw_counts.csv",header = TRUE, row.names=1 ))


data.class(raw_counts)
data.frame(raw_counts)
str(raw_counts)
summary(raw_counts)


data.matrix(raw_counts)

samples <- read.csv("sample_list.csv", row.names=1)
samples <- samples[c("Symbiotic","Mineral")]
all(rownames(samples) %in% colnames(raw_counts))

raw_counts<- raw_counts[, rownames(samples)]
all(rownames(samples) == colnames(raw_counts))


set.seed(12345)

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = samples,
                              design= ~ Mineral+ Symbiotic + Mineral:Symbiotic)

dds$Mineral <- relevel(dds$Mineral, ref="N","Y")
dds$symb <- relevel(dds$Symbiotic, ref="NS","S")

data.class(dds)


dds2<-dds [rowSums(counts(dds))>=50,]

k <- dds2[, -5]
y <- dds2[, 5]
r <- 3

model <- lfda(k,y , r=4, metric = "plain")
predict(model,dds2[, -5])

transformedData <- predict(model, raw_counts2[, -5])
plot(x = model, labels = raw_counts[, 5])

percentVar <- lfda(dds2[-5], dds2[, 5], r = 3, metric="plain")
autoplot(percentVar, data = dds2, frame = TRUE, frame.colour = 'Symbiotic')


####scree plot

install.packages("BiocManager")
library("BiocParallel")
library("DESeq2")
library("ggplot2")
library("magrittr")
library("matrixStats")

BiocManager::install("glmGamPoi")


setwd("~/Desktop/weathering paper_22/RNA_seq")
raw_counts <- read.csv("raw_counts.csv",header = TRUE, row.names=1 )
raw_counts <- as.data.frame(read.csv("raw_counts.csv",header = TRUE, row.names=1 ))


data.class(raw_counts)
data.frame(raw_counts)
str(raw_counts)
summary(raw_counts)


data.matrix(raw_counts)

samples <- read.csv("sample_list.csv", row.names=1)
samples <- samples[c("Symbiotic","Mineral")]
all(rownames(samples) %in% colnames(raw_counts))

raw_counts<- raw_counts[, rownames(samples)]
all(rownames(samples) == colnames(raw_counts))


set.seed(12345)

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = samples,
                              design= ~ Mineral + Symbiotic + Mineral:Symbiotic)

dds$Mineral <- relevel(dds$Mineral, ref="N","Y")
dds$symb <- relevel(dds$Symbiotic, ref="NS","S")

data.class(dds)


dds2<-dds [rowSums(counts(dds))>=50,]

transformedreadcounts <- rlog(dds2, blind=TRUE)

data.class(transformedreadcounts)


groups <- c("NS_HAP_minus1",	"NS_HAP_minus2",	"NS_HAP_minus3",	"NS_HAP_minus4",	"S_HAP_minus1",	"S_HAP_minus2",	"S_HAP_minus3",	"S_HAP_minus4",	"NS_HAP_plus1",	"NS_HAP_plus2",	"NS_HAP_plus3",	"NS_HAP_plus4",	"S_HAP_plus1",	"S_HAP_plus2",	"S_HAP_plus3",	"S_HAP_plus4")

PCA <- plotPCA(transformedreadcounts, intgroup =c("Symbiotic","Mineral"), returnData=TRUE)

percentVar <- round(100 * attr(PCA, "percentVar"))
ggplot(PCA, aes(PC1, PC2, color=Mineral, shape=Symbiotic)) +
  geom_point(size=3) +
  geom_text(aes(label=rownames(samples)),hjust=0.5, vjust=-1)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() 

## calculate the variance for each gene
variance <- rowVars(assay(transformedreadcounts))

## select the ntop genes by variance
topgenes <- order(variance, decreasing=TRUE)[seq_len(min(500, length(variance)))]

## perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(transformedreadcounts)[topgenes,]))

## the contribution to the total variance for each component
percentVar <- ((pca$sdev^2) / (sum(pca$sdev^2)))*100

##plot the "percentVar" of all components 
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:16)

colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity") 

barplot(percentVar, cex.names=1, xlab=paste("Principal component (PC), 1-", length(pca$sdev)), ylab="Proportion of variation (%)", main="", ylim=c(0,50))
 

####make a 3 components plot_NEEDS FIX

as.data.frame(topgenes)

PCA <- plot(topgenes, intgroup =c("Symbiotic","Mineral"), returnData=TRUE)
percentVar <- round(100 * attr(PCA, "percentVar"))


ggplot(PCA, aes("PC1", "PC2", "PC3", color=Mineral, shape=Symbiotic)) +
  geom_point(size=3) +
  geom_text(aes(label=rownames(groups)),hjust=0.5, vjust=-1)+
  xlab("PC1",percentVar[PC1],"% variance") +
  ylab("PC2",percentVar[PC2],"% variance") + 
  zlab("PC3",percentVar[PC3],"% variance")+
  coord_fixed()


####### Regularized log transformation for PCA
# The regularized log transform can be obtained using the rlog() function. 
# Regularized log transform is to stabilize the variance of the data and to make its distribution roughly symmetric
# Note that an important argument for this function is blind (TRUE by default). 
# The default “blinds” the normalization to the design. 
# This is very important so as to not bias the analyses (e.g. class discovery)



library(RcppArmadillo)
library(colorspace)
library(lattice)
library(RODBC)
library(Matrix)
library(survival)
library(Rcpp)
library(genefilter)
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(plyr)
library(DESeq2)
library(RColorBrewer)
library(stringr)
library(biomaRt)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pathview)
library(gage)
library(gageData)
library(Biobase)
library(S4Vectors)
library(stats4)
library(BiocGenerics)
library(parallel)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(SummarizedExperiment)


rld=rlog(dds2,blind=TRUE)


plotPCA(rld, intgroup = c("Symbiotic","Mineral"),ntop = 50000)
data <- plotPCA(rld, intgroup = c( "Mineral","Symbiotic"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
black.bold.16.text <- element_text(face = "bold", color = "black", size = 16)
ggplot(data=data, aes_string(x="PC1", y="PC2", color="Mineral",shape="Symbiotic")) + 
  geom_point(size=5) + 
  geom_text(aes(label=rownames(samples)),hjust=0.5,vjust=-1)+
  theme_bw() + 
  xlim(-25, 25) + 
  ylim(-20, 15) +
  theme(text = black.bold.16.text, 
        axis.text = black.bold.16.text,
        axis.line.x = element_line(color="black", size=0.25),
        axis.line.y = element_line(color="black", size=0.25),
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank(),
        legend.position=c(1,0),
        legend.justification=c(1.02,-0.02),) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) 


## Print 3D PCA plot

library(plotly)

plotPCA3D<-function (rld, intgroup = "samples", ntop = 500, returnData = FALSE)
  
  
  {rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(rld)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(rld)))) 
    {stop("the argument 'intgroup' should specify columns of colData(dds)")}
  intgroup.df <- as.data.frame(colData(rld)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) 
    {factor(apply(intgroup.df, 1, paste, collapse = " : "))}
  else {colData(rld)[[intgroup]]}
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  PC3 = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name = colnames(rld))
  if (returnData) {attr(d, "percentVar") <- percentVar[1:3]
    return(d)}
  message("Generating plotly plot")
  p <- plotly::plot_ly(data = dds_norm,
                       x = ~PC1,
                       y = ~PC2,
                       z = ~PC3,
                       color = group,
                       mode = "markers",
                       type = "scatter3d")
  return(p)}


PCA3D<-plotPCA3D(rld, intgroup = c("Symbiotic","Mineral"), ntop = 500, returnData = TRUE)
write.csv(PCA3D, "PCA.csv")


install.packages("scatterplot3d")
library(scatterplot3d)
par(mar=c(4,4,4,4), cex=1.0, cex.main=0.8, cex.axis=0.8)



scatterplot3d(PCA3D[,1:3], angle=70, color="black",pch = 17,
              xlim = NULL, ylim= NULL, zlim =NULL,
              xlab=paste("PC1:", percentVar[1],"% variance"), 
              ylab=paste("PC2:", percentVar[2],"% variance"), 
              zlab=paste("PC3:", percentVar[3],"% variance"),
              grid=TRUE, box=TRUE)

library(plotly)


fig <- plot_ly(PCA3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~group, colors = c("orange","chocolate3","limegreen","chartreuse4"))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC1:43% variance'),
                                   yaxis = list(title = 'PC2:14% variance'),
                                   zaxis = list(title = 'PC3:13% variance')))

fig






