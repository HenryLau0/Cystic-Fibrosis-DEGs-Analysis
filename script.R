
# Differential expression analysis with limma
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("annotate")
BiocManager::install("hgu133a.db")
BiocManager::install("GEOquery")

library(affy)
library(limma)
library(annotate)
library(hgu133a.db)
library(GEOquery)

#Read the raw dataset
CF <- ReadAffy()
CF
summary(CF)
summary(exprs(CF))

#load data from NCBI GEO
# load series and platform data from GEO
my_id <- "GSE15568"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse

pData(gse)[1:2,] ## print the sample information
fData(gse)[1,] ## print the gene annotation
exprs(gse)[1,] ## print the expression data

gset <- getGEO("GSE15568", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Data pre-processing
# check whether the dataset is normalized or not
pData(gset)$data_processing[1]

# log2 transformation
exprs(gset) <- log2(exprs(gset))

# check again the summary
summary(exprs(gset))

# perform RMA normalization
eset <- rma(CF)
write.exprs(eset, file="RMAnormalized.txt") 

#plot boxplot
colors = c(rep("lightpink",16), rep("lightblue",13))
par(mfrow=c(2,2))
boxraw <- boxplot(CF, col=colors, outline = FALSE, cex.axis=".8", las="3", 
                  main="Boxplot RAW")
boxnorm <- boxplot(log2(exprs(eset)), col=colors, outline = FALSE, cex.axis=".8", las="3",
                   main="Boxplot Normalized")
dev.off()

#plot histogram
colors = c(rep("lightpink",16), rep("lightblue",13))
par(mfrow=c(2,2))
histraw <- hist(CF, cex.axis=".8", las="3",  main="Histogram RAW")
histnorm <- plotDensity(exprs(gset), cex.axis=".8", las="3", xlab="log intensity", 
                        main="Histogram Normalised")
dev.off()

#mean intensity values of log2 transformed data (x-axis) and expression change (y-axis)
MAplot(CF)

#analysis of variance
#select CF and non-CF with CF disease
x <- c(rep(1,16), rep(2,13))
design <- model.matrix(~ -1+factor(x))
colnames(design) <- c("CF", "NonCF")
design

ID <- featureNames(gset)
Symbol <- getSYMBOL(ID, "hgu133a.db")
fData(gset) <- data.frame(ID=ID, Symbol=Symbol)

# fits a linear model for each gene based on the given series of arrays
fit <- lmFit(gset, design)

# dim shows dimension in rows and columns
dim(fit)

# colnames shows the column names
colnames(fit)

# display the first 10 rows of fit
rownames(fit)[1:10]

# display the last 10 rows of fit
tail(rownames(fit))

names(fit)

summary(fit)

#By topTable, 10 genes are selected with the smallest p-value adjusted for the 
#false discovery rate
contrast.matrix <- makeContrasts(CF-NonCF, levels=design)
contrast.matrix

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, 0.01)		
# computes empirical Bayes statistics for differential expression. 
limma_complete <- topTable(fit2, adjust="fdr", sort.by="B", number=62000)
# generates the list of differentially expressed genes sorted by B-values.
limma_complete[1:10,]

#How many DE genes are there?
length(which(limma_complete$adj.P.Val < 0.05))

#write a csv file
write.csv(limma_complete, file="limma_complete.csv")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# Q-Q plot
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)
