# Using EdgeR to perform diferential expression analysis

library (edgeR)

setwd("~/workspace/norovirus")

fi = read.delim("samples-counts.tsv", row.names="Symbol")
# Separating by spanish and brazilian samples
group = factor(c(1, 1, 1, 2, 2, 2, 2))

y = DGEList(counts=fi, group=group)
#y = calcNormFactors(y)
#design = model.matrix(~group)
#y = estimateDisp(y, design)

# Performing quasi-likelihood F-tests
#fit = glmQLFit(y, design)
#qlf = glmQLFTest(fit, coef=2)
#topTags(qlf)

# Filtering lowly expressed genes, this is because counts are so low, they will not be translated
#y$samples
#group lib.size norm.factors
#X3351     1    81805    1.5814930
#X3453     1   127389    1.0633715
#X3587     1   154109    0.5681635
#G53       2  1956747    0.7106078
#G71       2    48172    1.0264165
#H3        2   238020    1.2956658
#N99       2   141366    1.1074593

# Filter out lowly expressed genes 
# Wee keep genes that have one count per milion at least in two samples
keep = rowSums(cpm(y)>1) >=2
y = y[keep, , keep.lib.sizes=FALSE]

# Now the normalization factors are more normal
# $samples
# group lib.size norm.factors
# X3351     1    79667            1
# X3453     1   123760            1
# X3587     1   149692            1
# G53       2  1887245            1
# G71       2    47309            1
# H3        2   231080            1
# N99       2   137391            1

# Normalization, necessary for sample-specific effects
# We dont need to normalize by read length because the length is the same in all samples

# Perform TMM normalization, to normalize RNA composition by finding a set of scaling factors for the
# library sizes that minimize the log-fold changes between the samples for most genes
y = calcNormFactors(y)
y$samples

# Results:
# group lib.size norm.factors
# X3351     1    79667    1.7273815
# X3453     1   123760    1.0393038
# X3587     1   149692    0.4608158
# G53       2  1887245    0.7238030
# G71       2    47309    1.0262140
# H3        2   231080    1.4367111
# N99       2   137391    1.1326977

# A normalization facor below one indicates that a small number of high count genes are monopolizing the 
# sequencing, causing the counts of other genes to be lower than would be usual given the library size.
# As a result, the library size will be scaled down, analogous to scaling the counts upwards in that library.
# Conversely, a factor above one scales up the library size, analogous to downscaling the counts.

plotMDS(y)

# Estimate dispersions
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)

plotBCV(y)

# Perform statistical test
et = exactTest(y)
topTags(et)

#The total number of DE genes at 5% FDR is given by
summary(de <- decideTestsDGE(et))

# 1+2
# Down       0
# NotSig 15997
# Up         0

#Plot the log-fold-changes, highlighting the DE genes:
# The blue lines indicate 2-fold changes.
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")

# It is not significant



###########################################################
# Alejandros way
fi = read.delim("samples-counts.tsv", row.names="Symbol")
# Separating by spanish and brazilian samples
group = factor(c(1, 1, 1, 2, 2, 2, 2))

y = DGEList(counts=fi, group=group)
y = calcNormFactors(y)
design = model.matrix(~0+group)
y = estimateDisp(y, design)
plotMDS(y)

# Performing quasi-likelihood F-tests
fit = glmQLFit(y, design)
qlf = glmQLFTest(fit, coef=2)
topTags(qlf)
pru <- topTags(qlf,n=25000)
dat <- data.frame(pru)
up_reg <- subset(dat,logFC>1 & FDR<0.05)
down_reg <- subset(dat, logFC<(-1) & FDR<0.05)

# Testing another way
fi = read.delim("samples-counts.tsv", row.names="Symbol")
# Separating by spanish and brazilian samples
group = factor(c(1, 3, 1, 2, 2, 2, 1))

y = DGEList(counts=fi, group=group)
design = model.matrix(~0+group)

# Estimate dispersion of the dataset
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)

#Then we estimate gene-wise dispersion estimates, allowing a possible trend with averge count size:
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# determine differentially expressed genes.
fit <- glmFit(y, design)
# likelihood ratio tests for spain vs brazil  differences and show the top genes:
lrt <- glmLRT(fit)
topTags(lrt)
# counts-per-million in individual samples for the top genes
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
# total number of differentially expressed genes at 5% FDR
summary(de <- decideTestsDGE(lrt))
# group2
# Down    21872
# NotSig  38803
# Up          0
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")

# Not significant

