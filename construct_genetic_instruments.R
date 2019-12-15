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

snps.locs <- GRanges(
    seqnames = paste0('chr', SNPLocations$SNPChr),
    ranges = IRanges(
        start = SNPLocations$SNPPos,
        end = SNPLocations$SNPPos
    ),
    snp = SNPLocations$SNPName
)

# function to create genetic instruments per chromosome
run <- function(chr) {
    # genes to create GI for
    genes <- ensembl$EnsemblGeneID[ensembl$ChromosomeName == chr]
    genes.tmp <- genes.locs[genes.locs$gene %in% genes]
    
    # identify nearby genes (<1Mb)
    overlaps <- findOverlaps(genes.tmp, snps.locs, maxgap = 1e6)
    overlaps <- data.frame(
        gene = genes.tmp$gene[queryHits(overlaps)],
        snp = snps.locs$snp[subjectHits(overlaps)]
    )

    # load genotype data
    load(sprintf('%s/genotypes.chr.%s.Rdata', data_path, chr))

    # loop over genes to construct GI
    out <- lapply(genes, function(current.gene) {
        # store nearby genes in new object
        snps <- overlaps$snp[overlaps$gene == current.gene]

        # genotype data for this gene, training data
        geno <- t(genotypes[rownames(genotypes) %in% snps, cvrt.train$ids2, drop = FALSE])
        
        # gene expression, training data
        y <- t(M[current.gene, cvrt.train$run_id.y, drop = FALSE])

        # model matrix, training data
        mm <- model.matrix(~ . + 0, data = cvrt.train)
        x <- cbind(mm, geno)

        # penalty factors for lasso
        penalty.factors <- c(rep(0, ncol(mm)), rep(1, ncol(x)-ncol(mm)))

        # create folds
        set.seed(10)
        foldid <- sample(rep(1:5, length = nrow(cvrt.train)))
        
        # fit penalized regression model
        cv.fit <- cv.glmnet(x = x, y = y, foldid = foldid, type.measure = type.measure, penalty.factor = penalty.factors, parallel = TRUE)

        # identify variants with non-zero coefficients
        beta <- drop(as.matrix(coef(cv.fit, s = 'lambda.min')))
        beta <- beta[beta != 0]
        snps <- names(beta)

        # check if any variants have been selected
        if (length(snps) > 0) {

            # genotype data, test data
            geno <- t(genotypes[snps, cvrt.test$ids2, drop = FALSE])
            
            # gene expression, test data
            y <- t(M[current.gene, cvrt.test$run_id.y, drop = FALSE])

            # create GI, test data
            gi <- drop(geno[, snps, drop = FALSE] %*% cbind(beta[snps]))

            # calculate F-statistic
            # fit reduced model and full model
            fit.reduced <- lm(y ~ ., data = cvrt.test)
            fit.full <- lm(y ~ . + gi, data = cvrt.test)

            # calculate F-value
            F.value <- anova(fit.reduced, fit.full)[2, 'F']

            ## predict data for all samples using model
            geno <- t(genotypes[snps, cvrt$ids2, drop = FALSE])
            y <- t(M[current.gene, cvrt$run_id.y, drop = FALSE])
            gi <- drop(geno[, snps, drop = FALSE] %*% cbind(beta[snps]))
 
            # create data frame with relevant data
            predicted.data <- data.frame(
                run_id.y = cvrt$run_id.y,
                gene = current.gene,
                y = drop(y),
                gi = gi
            )

            # save results
            save(F.values, predicted.data, file = sprintf('%s/%s.gi.Rdata', out_path, current.gene))
        }
    })
}

chrs <- 1:22
max.cores <- 30
nworkers <- min(max.cores, length(chrs))
BPPARAM <- MulticoreParam(workers = nworkers, verbose = TRUE)
empty <- bplapply(chrs, FUN = run, BPPARAM = BPPARAM)
