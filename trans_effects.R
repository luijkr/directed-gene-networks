options('stringsAsFactors' = FALSE)
library(BiocParallel)
library(GenomicRanges)
library(glmnet)

# set variables
data_path = '/../..'
out_path = '/../..'

nfolds <- 5 # K-fold cross-validation
type.measure <- 'mse' # loss function
alpha <- 1 # lasso

# read gene annotation
ensembl <- read.table(sprintf('%s/ensembl.txt', data_path), header = TRUE, sep = ',')
ensembl <- subset(ensembl,
    Genetype == 'protein_coding' &
    !duplicated(EnsemblGeneID) &
    ChromosomeName %in% 1:22
)

# load genotype annotation
load(sprintf('%s/SNPLocations.Rdata', data_path))

# load RNA-seq data
load(sprintf('%s/RNAseq.F2.Rdata', data_path))

# create GRanges objects
genes.locs <- GRanges(
    seqnames = paste0('chr', ensembl$ChromosomeName),
    ranges = IRanges(
    start = ensembl$GeneStart,
    end = ensembl$GeneEnd
    ),
    gene = ensembl$EnsemblGeneID
)
genes.locs <- genes.locs[genes.locs$gene %in% rownames(M)]

# read F-values
F.values <- read.table(sprintf('%s/fstats.txt', out_path), header = TRUE, sep = '\t')

# load genetic instruments
rdata.files <- list.files(out_path, pattern = 'ENSG*.gi.Rdata', full.names = TRUE)

gis <- bplapply(F.values$gene, FUN = function(gene) {
    file <- grep(gene, rdata.files, value = TRUE)
    load(file)
    out <- cbind(predicted.data$gi)
    colnames(out) <- gene
    rownames(out) <- predicted.data$run_id.y
    return(out)
}, BPPARAM = MulticoreParam(workers = 10, verbose = TRUE))

gis <- do.call(cbind, gis)
gis <- gis[match(cvrt$run_id.y, rownames(gis)), , drop = FALSE]

# transpose data
Y <- t(M)
rm(M); gc()

# define model matrix for covariates
mm <- model.matrix(~ ., data = cvrt)

# function to fit linear regression models for each index gene
run <- function(current.gene) {
    ## determine nearby genes to correct for
    nearby.genes <- subsetByOverlaps(ensembl, ensembl[ensembl$gene == current.gene], maxgap = 1e6)$gene
    nearby.genes <- setdiff(nearby.genes, current.gene)
    nearby.genes <- intersect(nearby.genes, colnames(gis)) # make sure they have a GI

    ## update model matrix by adding nearby index genes with a GI
    Z <- cbind(mm, gis[cvrt$run_id.y, nearby.genes, drop = FALSE])

    # code below is based on Sikorska et al. (2013)
    # https://doi.org/10.1186/1471-2105-14-166
    k <- ncol(Z) - 1
    n <- nrow(Y)

    U1 <- crossprod(Z, Y)
    U2 <- solve(crossprod(Z), U1)
    ytr <- Y - Z %*% U2

    X <- gis[cvrt$run_id.y, current.gene, drop = FALSE]
    U3 <- crossprod(Z, X)
    U4 <- solve(crossprod(Z), U3)
    Xtr <- X - Z %*% U4
    Xtr2 <- colSums(Xtr**2)

    b <- crossprod(ytr, Xtr)
    b <- b / matrix(Xtr2, nr = nrow(b), nc = ncol(b), byrow = TRUE)
    term1 <- colSums(ytr^2)
    term2 <- matrix(Xtr2, nc = ncol(b), nr = nrow(b), byrow = TRUE) * (b**2)
    term1 <- drop(term1)
    term2 <- drop(term2)
    sig <- (term1 - term2) / (n-k-2)
    err <- sqrt(sig * (1 / Xtr2))
    
    # create data frame with estimates and SE
    me <- data.frame(
        gi = current.gene,
        target = names(err),
        estimate = drop(b),
        se = drop(err)
    )
    
    # remove association with nearby genes (<10Mb)
    cis_genes <- subsetByOverlaps(ensembl, ensembl[ensembl$gene == current.gene], maxgap = 10e6)$gene
    me <- subset(me, !target %in% cis_genes)

    write.table(me, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t',
        file = sprintf('%s/%s.results.txt', current.gene))
}

# select genes to test
min.F <- 10 # minimum F-statistic
genes <- F.values$gene[F.values$f.new >= min.F]

max.cores <- 20
BPPARAM <- MulticoreParam(workers = min(max.cores, length(genes)), verbose = TRUE)
null <- bplapply(genes, FUN = run, BPPARAM = BPPARAM)
