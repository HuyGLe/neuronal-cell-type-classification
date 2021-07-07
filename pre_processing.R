library("edgeR")
library("ggplot2")
library("matrixStats")
library(stringr)
library(DescTools)
memory.limit(60000)

setwd("C:/Users/Thomas Le/Desktop/Senior Design")

# Loads in data
intron_matrix <- read.csv("human_MTG_2018-06-14_intron-matrix.csv")
exon_matrix <- read.csv("human_MTG_2018-06-14_exon-matrix.csv")
full_matrix <- intron_matrix + exon_matrix


full_matrix_log2cpm <- log2(cpm(full_matrix)+1)

# Replace first row of gene id's
full_matrix[,1] <- intron_matrix[,1]
full_matrix_log2cpm[,1] <- intron_matrix[,1] 

# ------------------------------- External CSVs ------------------------------- #

# Get cell labels from external csv
# NOTE: This code was for only the highest cell subtype
# labels <- read.csv("final_two_sample_columns.csv")
# cell.types <- c("Exc", "Inh", "Oli", "Ast", "OPC", "End", "Mic")

# Create cluster columns
samples <- read.csv("human_MTG_2018-06-14_samples-columns.csv", stringsAsFactors = FALSE)
labels <- samples[c("sample_name", "cluster")]
labels <- labels[labels$cluster != "no class",]
labels <- cbind(labels, str_split_fixed(labels$cluster, " ", n=Inf))
colnames(labels) <- c("sample_name", "cluster", "higher", "layer", "intermediate", "granular")

cell.types <- unique(labels$cluster)

# Gets list of column indices for a particular cell type to subset
for(type in cell.types) {
  t(assign(paste0(type, ".celltypes"), labels[which(labels$clu == type),1]))
}

# Read in gene.rows file
gene.rows <- read.csv("human_MTG_2018-06-14_genes-rows.csv")

# Replace row names of full matrix with genes
rownames(full_matrix_log2cpm) <- gene.rows[["gene"]]


# ------------------------------- Main Matrix 0's ------------------------------- #

# Remove genes from whole matrix w/ 0 expression
keep <- rowSums(full_matrix_log2cpm[, 2:15929]) > 0
kept.genes <- full_matrix_log2cpm[keep,]


# Subset cell subtypes
for (n in cell.types) {
  temp.cell.types <- get(paste0(n, ".celltypes"))
  temp.subset <- kept.genes[, temp.cell.types]
  assign(paste0(n, ".mm.subset"), temp.subset)
}

# Create empty median common gene matrix
cg.median.mm.df <- data.frame(matrix(ncol = 0, nrow = length(kept.genes[,1])))
rownames(cg.median.mm.df) <- rownames(kept.genes)

# Create empty count common gene matrix
cg.count.mm.df <- data.frame(matrix(ncol = 0, nrow = length(kept.genes[,1])))
rownames(cg.count.mm.df) <- rownames(kept.genes)

# Fill median common gene matrix with row medians from every cell subtype
for (n in cell.types) {
  temp.subset.common <- get(paste0(n, ".mm.subset"))
  temp.subset.median.common <- rowMedians(temp.subset.common)
  cg.median.mm.df <- cbind(cg.median.mm.df, temp.subset.median.common)
  temp.subset.count.common <- rowSums(temp.subset.common > 0)
  cg.count.mm.df <- cbind(cg.count.mm.df, temp.subset.count.common)
}

colnames(cg.median.mm.df) <- cell.types

# Calculate variance of the medians and remove genes w/ 0 variance
cg.median.mm.var <- apply(cg.median.mm.df, 1, var)
keep.var <- which(cg.median.mm.var != 0)
cg.median.mm.var <- cg.median.mm.var[cg.median.mm.var != 0]

# Subset list of genes after removing variance criteria
cg.median.mm.df.var <- cg.median.mm.df[keep.var, ]


# This line of code takes a long time
full_matrix_log2cpm_subset <- rbind(subset(full_matrix_log2cpm, select = -X), samples[["cluster"]])

keep.var.genes <- names(keep.var)
rownames(full_matrix_log2cpm_subset)[50282] <- "Classification"
final.matrix <- full_matrix_log2cpm_subset[c(keep.var.genes, "Classification"),]

write.csv(final.matrix, file = "everything.csv", row.names = TRUE, col.names = TRUE)


set.seed(916)
trainingv <- vector()
testv <- vector()
validationv <- vector()
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".celltypes"))
  cut1 <- as.integer(length(temp.subset) * 0.6)
  cut2 <- as.integer(length(temp.subset) * 0.8 )
  temp.subset <- temp.subset[sample(1:length(temp.subset))]
  trainingv <- append(trainingv, temp.subset[0:cut1])
  testv <- append(testv, temp.subset[cut1:cut2])
  validationv <- append(validationv, temp.subset[cut2:length(temp.subset)])
}
  

training <- t(final.matrix[,trainingv])
test <- t(final.matrix[,testv])
validation <- t(final.matrix[,validationv])


write.csv(training, file = "training.csv", row.names = TRUE, col.names = TRUE)
write.csv(test, file = "test.csv", row.names = TRUE, col.names = TRUE)
write.csv(validation, file = "validation.csv", row.names = TRUE, col.names = TRUE)


### ENTROPY ###

# Create empty entropy common gene matrix
cg.median.entropy.mm.df <- data.frame(matrix(ncol = 0, nrow = length(keep.var)))
rownames(cg.median.entropy.mm.df) <- names(keep.var)

# Calculating the entropy for each gene
for (row in 1:nrow(cg.median.mm.df.var)) {
  entropy <- Entropy(as.numeric(cg.median.mm.df.var[row, ]))
  cg.median.entropy.mm.df[row,1] <- entropy
}

colnames(cg.median.entropy.mm.df) <- "Entropy"

# Plot histogram of entropy
hist(cg.median.entropy.mm.df[["Entropy"]], breaks = 100, xlab = "Entropy", ylab = "Frequency", main = "Entropy Distribution w/o NA")

### Coefficient of Variation ###

# Calculate coefficient of variation from median common gene matrix
cg.median.mm.mean <- apply(cg.median.mm.df.var, 1, mean)
cg.median.mm.cov <- cg.median.mm.var / cg.median.mm.mean

# Take care of nan values
cg.median.mm.cov[is.nan(cg.median.mm.cov)] <- 0


# Determine cov threshold
list_histo <- hist(cg.median.mm.cov[cg.median.mm.cov > 0], breaks=500, xlim=c(0,5), ylim=c(0,500), main = "Histogram of Coefficient of Variation", xlab = "Coefficient of Variation", freq=TRUE)
abline(v=median(cg.median.mm.cov), col="blue")
abline(v=thresh, col="red")
legend(3, 400, legend = c("Threshold", "Median"), col=c("red", "blue"), lty=1:1)

thresh <- list_histo$breaks[which(list_histo$counts == max(list_histo$counts))]

cov.kept.genes <- names(cg.median.mm.cov[cg.median.mm.cov > thresh])



# This line of code takes a long time
full_matrix_log2cpm_subset <- rbind(subset(full_matrix_log2cpm, select = -X), samples[["cluster"]])

rownames(full_matrix_log2cpm_subset)[50282] <- "Classification"
cov.final.matrix <- full_matrix_log2cpm_subset[c(cov.kept.genes, "Classification"),]



# Split into training, validation, test in 60, 20, 20

set.seed(916)
trainingv <- vector()
testv <- vector()
validationv <- vector()
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".celltypes"))
  cut1 <- as.integer(length(temp.subset) * 0.6)
  cut2 <- as.integer(length(temp.subset) * 0.8 )
  temp.subset <- temp.subset[sample(1:length(temp.subset))]
  trainingv <- append(trainingv, temp.subset[0:cut1])
  testv <- append(testv, temp.subset[cut1:cut2])
  validationv <- append(validationv, temp.subset[cut2:length(temp.subset)])
}





training <- t(cov.final.matrix[,trainingv])
test <- t(cov.final.matrix[,testv])
validation <- t(cov.final.matrix[,validationv])


write.csv(training, file = "cov_training.csv", row.names = TRUE, col.names = TRUE)
write.csv(test, file = "cov_test.csv", row.names = TRUE, col.names = TRUE)
write.csv(validation, file = "cov_validation.csv", row.names = TRUE, col.names = TRUE)



### Binary Score ###

final.m <- matrix(0, ncol = dim(cg.median.mm.df.var)[2], nrow = dim(cg.median.mm.df.var)[1])


for (val in 1:length(cell.types)) {
  temp.target <- replicate(length(cell.types), cg.median.mm.df.var[, val])
  temp.m <- cg.median.mm.df.var / temp.target
  temp.m[is.na(temp.m)] <- 1
  temp.m[temp.m > 1] <- 1
  sum.v <- rowSums(temp.m)
  final.m[,val] <- sum.v 
}


seventy5.matrix <- matrix(length(cell.types), ncol = dim(cg.median.mm.df.var)[2], nrow = dim(cg.median.mm.df.var)[1])
binary.score <- (seventy5.matrix  - final.m) / (dim(cg.median.mm.df.var)[2] - 1)
binary.score.max <- apply(binary.score, 1, max)

colnames(binary.score) <- colnames(cg.median.mm.df.var)
rownames(binary.score) <- rownames(cg.median.mm.df.var)



number.of.ones <- binary.score == 1

final.ones <- colSums(number.of.ones)


median.binary.score <- apply(binary.score,2,function(x){median(x[x>0])})


mean.binary.score <- apply(binary.score,2,function(x){mean(x[x>0])})



# Graph binary score distribution for each cluster
for (n in 1:length(cell.types)) {
  png(filename=paste0("C:/Users/Thomas Le/Desktop/Senior Design/Binary Score Histograms/", cell.types[n], " binary score hist.png"))
  hist(binary.score[,n], breaks=150, xlab="Binary Score Range", ylab="Frequency", main=paste0("Cluster ", cell.types[n]), ylim=c(0,150))
  abline(v=median.binary.score[n], col="blue")
  abline(v=mean.binary.score[n], col="red")
  dev.off()
}




## Comparing Brian's Binary Score's with Ours ##

# Rename columns of binary score
temp.names <- colnames(binary.score)
temp.names <- gsub("-", " ", temp.names)

# Make copy of binary score table for comparison
binary.score.temp <- binary.score
colnames(binary.score.temp) <- temp.names


# Get filtered sup3 table and rename cluster column
sup.table <- read.csv("sup_table3_filtered.csv")
sup.table.names <- sup.table$Cluster
sup.table.names <- gsub("_", " ", sup.table.names)
sup.table$Cluster <- sup.table.names

# Rename columns of binary score
temp.gene.names <- sup.table$Gene
temp.gene.names <- gsub("_", "-", temp.gene.names)
temp.gene.names <- gsub("\\*", "", temp.gene.names)

sup.table$Gene <- temp.gene.names

sup.table$Our_Binary_Score = NA

for (n in 1:dim(sup.table)[1]) {
  bs.col <- sup.table$Cluster[n]
  bs.row <- sup.table$Gene[n]
  our.binary.score <- binary.score.temp[bs.row, bs.col]
  sup.table$Our_Binary_Score[n] <- our.binary.score
}

sup.table$Diff <- sup.table$Binary_Score - sup.table$Our_Binary_Score

View(sort(abs(sup.table$Diff)))



# OPC: PCDH15, PDGFRA
# Oligo: ENPP2, CERCAM




# Calculate precursor matrices for binary score calculation
# seventy5.matrix <- matrix(length(cell.types), ncol = length(cell.types), nrow = dim(cg.median.mm.df.var)[1])
# gene.row.sums <- replicate(length(cell.types), apply(cg.median.mm.df.var, 1, sum))



# binary.score = (seventy5.matrix - (gene.row.sums/cg.median.mm.df.var)) / (length(cell.types) - 1)

# Calculate max value of binary score across each gene
# binary.score.max <- apply(binary_score, 1, max)

# Graph the binary score
# hist(binary.score.max, breaks = 100, xlab = "Binary Score", ylab = "Frequency", main = "Binary Score Distribution")


# BEVERLY'S Top 25% method
# top.25.list <- c()
# for (col in 1:length(cell.types)) {
#   temp.v <- binary.score[,col]
#   temp.v <- temp.v[0:length(temp.v)]
#   top.25.list <- c(top.25.list, temp.v)
# }
# top.25.list <- unique(top.25.list)
# 
# temp.v <- data.frame(0, nrow = 13945, ncol = 2)
# temp.v[,1] <- binary.score[,1]
# temp.v[,2] <- rownames(binary.score)
# temp.v[,1] <- sort(-temp.v[,1])
# names <- rownames(temp.v)
# top.25.list <- c(top.25.list, )




### Minimum set of Marker Genes ###

# Getting minimum set of marker gene
marker.genes.ref <- read.csv("Supplemental Table 2.csv", stringsAsFactors = FALSE)

rownames(marker.genes.ref) <- sort(cell.types)

marker.genes.ref <- marker.genes.ref[, c("X1", "X2", "X3", "X4", "X5", "X6")]
marker.genes.list <- c(marker.genes.ref[,1], marker.genes.ref[,2], marker.genes.ref[,3], marker.genes.ref[,4], marker.genes.ref[,5], marker.genes.ref[,6])



# Get unique genes from csv
marker.genes.list <- unique(marker.genes.list)
marker.genes.list <- marker.genes.list[1:158]
marker.genes.list <- marker.genes.list[marker.genes.list != ""]
marker.genes.list[marker.genes.list == "IFNG_AS1"] <- "IFNG-AS1"

# Check that all genes from the csv are in our pre-processed data
s.marker.genes.list <- rownames(cg.median.mm.df.var)
length(setdiff(marker.genes.list, s.marker.genes.list))



### Graphing marker genes from csv ###

marker.genes.subset <- full_matrix_log2cpm[which(rownames(full_matrix_log2cpm) %in% marker.genes.list), ]

# Creating 75 cell types subset matrices from the reference (mm.ref.subset)
for (n in cell.types) {
  temp.cell.types <- get(paste0(n, ".celltypes"))
  temp.subset <- marker.genes.subset[, temp.cell.types]
  rownames(temp.subset) <- rownames(marker.genes.subset)
  assign(paste0(n, ".mm.ref.subset"), temp.subset)
}

# Create empty median common gene data frame (75 columns)
cg.median.ref.df <- data.frame(matrix(ncol = 0, nrow = length(marker.genes.subset[,1])))
rownames(cg.median.ref.df) <- rownames(marker.genes.subset)

# Fill median common gene matrix with row medians from every cell subtype
for (n in cell.types) {
  temp.subset.common <- get(paste0(n, ".mm.ref.subset"))
  temp.subset.median.common <- rowMedians(temp.subset.common)
  cg.median.ref.df <- cbind(cg.median.ref.df, temp.subset.median.common)
}

colnames(cg.median.ref.df) <- cell.types


# Calculate variance of the medians and remove genes w/ 0 variance
cg.median.ref.var <- apply(cg.median.ref.df, 1, var)




# Create empty entropy common gene matrix
cg.median.entropy.ref.df <- data.frame(matrix(ncol = 0, nrow = length(marker.genes.subset[,1])))
rownames(cg.median.entropy.ref.df) <- rownames(marker.genes.subset)

# Calculating the entropy for each gene
for (row in 1:nrow(cg.median.ref.df)) {
  entropy <- Entropy(as.numeric(cg.median.ref.df[row, ]))
  cg.median.entropy.ref.df[row,1] <- entropy
}

colnames(cg.median.entropy.ref.df) <- "Entropy"

hist(cg.median.entropy.ref.df[["Entropy"]], breaks = 100, xlab = "Entropy", ylab = "Frequency", main = "Reference Entropy Distribution w/o NA")



# Calculate coefficient of variance from median common gene matrix
cg.median.ref.mean <- apply(cg.median.ref.df, 1, mean)
cg.median.ref.cov <- cg.median.ref.var / cg.median.ref.mean


hist(cg.median.ref.cov, breaks=500, xlim=c(0,5), ylim=c(0, 10), main = "Histogram of Coefficient of Variance", xlab = "Reference Coefficient of Variance", freq=TRUE)


# Binary Score

# Graph binary score distribution for each cluster with marker gene labels

cutoff = round(13945 * 0.15)

for (n in 1:length(cell.types)) {
  gene.list <- as.character(marker.genes.ref[cell.types[n],])
  gene.list <- gene.list[gene.list != ""]
  gene.list <- gene.list[gene.list != "NA"]
  gene.list[gene.list == "IFNG_AS1"] <- "IFNG-AS1"
  png(filename=paste0("C:/Users/Thomas Le/Desktop/Senior Design/Binary Score Histograms/", cell.types[n], " binary score hist.png"))
  hist(binary.score[,n], breaks=150, xlab="Binary Score Range", ylab="Frequency", main=paste0("Cluster ", cell.types[n]), ylim=c(0,150))
  for (gene in gene.list) {
    abline(v = binary.score[gene, n], col="blue")
    # print(binary.score[gene, n])
  }
  temp.binary.score = sort(binary.score[,n], decreasing = TRUE)
  abline(v=temp.binary.score[cutoff], col="pink")
  abline(v=median.binary.score[n], col="green")
  abline(v=mean.binary.score[n], col="red")
  dev.off()
}

clu.histogram = apply(binary.score, 2, var)



hist(binary.score.max.ref, breaks = 100, xlab = "Binary Score", ylab = "Frequency", main = "Reference Binary Score Distribution")

### End of graphing marker genes from csv ###


### Calculate Subset Sizes ###

subset.sizes = data.frame(matrix(ncol = 2, nrow = 0))

for (n in cell.types) {
  temp.subset.common <- get(paste0(n, ".mm.subset"))
  temp.size <- dim(temp.subset.common)[2]
  subset.sizes = rbind(subset.sizes, c(n, temp.size))
}


write.csv(subset.sizes, file = "cov_training.csv", row.names = TRUE, col.names = FALSE)




# ----------------------------- Individual Matrix 0's ----------------------------------- #


common.genes <- as.character(full_matrix$X)

# Subsets all cell types from full matrix
for (n in cell.types) {
  temp.cell.types <- get(paste0(n, ".celltypes"))
  temp.subset <- full_matrix_log2cpm[, temp.cell.types]
  rownames(temp.subset) <- intron_matrix[["X"]]
  # temp.subset.subtract <- temp.subset[!rowSums(temp.subset != 0),]
  # temp.subset <- temp.subset[!!rowSums(temp.subset != 0),]
  # genes.subtract <- rownames(temp.subset.subtract)
  # common.genes <- common.genes[! common.genes %in% genes.subtract]
  assign(paste0(n, ".subset"), temp.subset)
}

### Calculate Matrix Statistics ###

# Calculate means for every row
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  temp.mean <- rowMeans(temp.subset)
  assign(paste0(n, ".row.means"), temp.mean)
  # png(filename=paste0("C:/Users/Thomas Le/Desktop/Senior Design/", n, "_means_barplot.png"))
  # barplot(get(paste0(n, ".row.means")), xlab = "Genes", ylab = "Average expression (log2 CPM)", main = n)
  # dev.off()
}

# Calculate medians for every row
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  temp.median <- rowMedians(temp.subset)
  assign(paste0(n, ".row.medians"), temp.median)
  # png(filename=paste0("C:/Users/Thomas Le/Desktop/Senior Design/", n, "_medians_barplot.png"))
  # barplot(get(paste0(n, ".row.medians")), xlab = "Genes", ylab = "Median expression (log2 CPM)", main = n)
  # dev.off()
}

# Counts how many non-zeros are in each row
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  temp.count <- rowSums(temp.subset > 0)
  assign(paste0(n, ".row.counts"), temp.count)
  # png(filename=paste0("C:/Users/Thomas Le/Desktop/Senior Design/", n, "_counts_barplot.png"))
  # barplot(get(paste0(n, ".row.counts")), xlab = "Genes", ylab = "Number of cells expressed", main = n)
  # dev.off()
}

### Common Gene Pre-Processing ###

# Add in gene IDs as the first column
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  temp.subset <- cbind(full_matrix[,1], temp.subset)
  assign(paste0(n, ".subset"), temp.subset)
}

# Subset common genes from clusters
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  assign(paste0(n, ".subset.common"), subset(temp.subset, rownames(temp.subset) %in% common.genes))
}

# Create empty common gene matrices
cg.mean.df <- data.frame(matrix(ncol = 0, nrow = length(common.genes)))
rownames(cg.mean.df) <- common.genes

cg.median.df <- data.frame(matrix(ncol = 0, nrow = length(common.genes)))
rownames(cg.median.df) <- common.genes

cg.count.df <- data.frame(matrix(ncol = 0, nrow = length(common.genes)))
rownames(cg.count.df) <- common.genes

# Fill the common gene subset matrices
for (n in cell.types) {
  temp.subset.common <- get(paste0(n, ".subset.common"))
  temp.subset.mean.common <- rowMeans(temp.subset.common)
  cg.mean.df <- cbind(cg.mean.df, temp.subset.mean.common)
                   
  temp.subset.median.common <- rowMedians(temp.subset.common)
  cg.median.df <- cbind(cg.median.df, temp.subset.median.common)
  
  temp.subset.count.common <- rowSums(temp.subset.common > 0)
  cg.count.df <- cbind(cg.count.df, temp.subset.count.common)
}

# Rename the columns
colnames(cg.mean.df) <- cell.types
colnames(cg.median.df) <- cell.types
colnames(cg.count.df) <- cell.types

# Calculate row-wise variance
cg.mean.var <- apply(cg.mean.df, 1, var)
cg.median.var <- apply(cg.median.df, 1, var)
cg.count.var <- apply(cg.count.df, 1, var)

# Calculate coefficient of variance for cg.median
cg.median.mean <- apply(cg.median.df, 1, mean)
cg.median.cov <- cg.median.var / cg.median.mean

# Fill nan values w/ 0
cg.median.cov[is.nan(cg.median.cov)] <- 0


# Get list of genes to remove with variance < 0.1 in cg.median
cg.median.filtered <- cg.median.cov[cg.median.cov < 0.1]
genes.to.remove <- names(cg.median.filtered)

# Keeps track of number of genes before pre-processing
number.genes.original <- vector()
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  number.genes.original <- cbind(number.genes.original, dim(temp.subset)[1])
}

colnames(number.genes.original) <- cell.types

# Removes mid-var genes from all cell subtypes
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  temp.genes <- rownames(temp.subset)
  temp.genes <- temp.genes[!temp.genes %in% genes.to.remove]
  temp.subset <- temp.subset[temp.genes,]
  assign(paste0(n, ".subset"), temp.subset)
}

# Keeps track of number of genes after performing mid-var
number.genes.midvar <- vector()
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  number.genes.midvar <- cbind(number.genes.midvar, dim(temp.subset)[1])
}

colnames(number.genes.midvar) <- cell.types

# Cut out genes which are expressed in less than 50% of cells
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset"))
  temp.row.counts <- rowSums(temp.subset > 0)
  keep <- which(temp.row.counts >= dim(temp.subset)[2]*0.5)
  temp.subset.f <- temp.subset[keep,]
  assign(paste0(n, ".subset.f"), temp.subset.f)
}

# Keeps track of number of genes after mid-var and 50%
number.genes.fifty.cutoff <- vector()
for (n in cell.types) {
  temp.subset <- get(paste0(n, ".subset.f"))
  number.genes.fifty.cutoff <- cbind(number.genes.fifty.cutoff, dim(temp.subset)[1])
}

colnames(number.genes.fifty.cutoff) <- cell.types

# Set column names to NULL for bar plots
colnames(number.genes.original) <- NULL
colnames(number.genes.midvar) <- NULL
colnames(number.genes.fifty.cutoff) <- NULL

# ------------------------------ Miscellaneous Graphs and Other Code ------------------------------------ #

# Bar plot change in number of genes
barplot(number.genes.original, ylim = c(0, 50000), xlab = "Genes", ylab = "Count", main = "Number of Genes Before Pre-Processing")
barplot(number.genes.midvar, ylim = c(0, 50000), xlab = "Genes", ylab = "Count", main = "Number of Genes After Mid-Var")
barplot(number.genes.fifty.cutoff, ylim = c(0, 50000), xlab = "Genes", ylab = "Count", main = "Number of Genes After 50% Cutoff")

# Removes variables
for (n in cell.types) {
rm(list=paste0(n, ".row.means"))
rm(list=paste0(n, ".row.medians"))
rm(list=paste0(n, ".row.counts"))
rm(list=paste0(n, ".row.counts.f"))
rm(list=paste0(n, ".subset"))
rm(list=paste0(n, ".subset.f"))
rm(list=paste0(n, ".subset.common"))
}

# Example Histogram plots for normally distributed
hist(`Astro L1-2 FGFR3 GFAP.row.medians`[`Astro L1-2 FGFR3 GFAP.row.medians` != 0], breaks=150, xlim=c(0,13), xlab = "Gene Median cpm Range", main = "Astro L1-2 FGFR3 GFAP")
hist(`Astro L1-2 FGFR3 GFAP.row.means`[`Astro L1-2 FGFR3 GFAP.row.means` != 0], breaks=250, xlim=c(0,5), xlab = "Gene Mean cpm Range", main = "Astro L1-2 FGFR3 GFAP")


# Graph median gene values between all clusters
median.matrix.filtered <- cg.median.mm.df[!names(cg.median.mm.cov) %in% genes.to.remove.mm,]

for (n in 1:30) {
  png(filename=paste0("C:/Users/Thomas Le/Desktop/Senior Design/Histograms/", rownames(median.matrix.filtered)[n], "_gene_medians_hist.png"))
  hist(as.matrix(median.matrix.filtered[n,]), breaks=35, main = paste0("Histogram of gene ", rownames(median.matrix.filtered)[n]))
  dev.off()
}

# Histogram plot of coefficient of variance to choose threshold
hist(cg.median.mm.cov, breaks=500, xlim=c(0,5), ylim=c(0,500), main = "Histogram of Coefficient of Variance", xlab = "Coefficient of Variance")
