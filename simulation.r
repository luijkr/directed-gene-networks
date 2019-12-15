library(BiocParallel)
library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)

generate_correlated_variable <- function(x, rho) {
  # Generate new variable y with a fixed correlation rho
  # to the existing variable x.
  
  N <- length(x)
  
  theta <- acos(rho) # corresponding angle
  x2    <- rnorm(N) # new random data
  X     <- cbind(x, x2)
  Xctr  <- scale(X, center = TRUE, scale = FALSE) # centered columns (mean 0)
  
  Id   <- diag(N) # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop = FALSE])) # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q) # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2] # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o) # bind to matrix
  Y    <- Xc2 %*% diag(1 / sqrt(colSums(Xc2^2))) # scale columns to length 1
  
  y <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1] # final new vector
  y <- drop(scale(y)) # scale final vector

  return(y)
}

simulate_data <- function(gi_rho, cis_rho, trans_rho, n, n_gis) {
  # generate GIs
  mu <- numeric(n_gis)
  Sigma <- matrix(gi_rho, nr = n_gis, nc = n_gis)
  diag(Sigma) <- 1
  GI <- mvrnorm(n, mu, Sigma, empirical = TRUE)

  # generate index gene
  gi <- GI[, 1]
  index_gene <- generate_correlated_variable(x = gi, rho = cis_rho)
  
  # generate target gene
  target_gene <- generate_correlated_variable(x = index_gene, rho = trans_rho)
  
  return(list(
    GIs = GI,
    index_gene = index_gene,
    target_gene = target_gene
  ))
}


run_simulation <- function(gi_rho, cis_rho, trans_rho, n = 3071, n_gis = 2) {
	# Run the simulation under different conditions.

  # generate data
  data <- simulate_data(gi_rho, cis_rho, trans_rho, n, n_gis)
  gi <- data$GIs[, 1]
  gi_neighbor <- data$GIs[, -1]
  index_gene <- data$index_gene
  target_gene <- data$target_gene

  
  # Scenario 1A: index gene causally affects target gene
  # Uncorrected for LD/pleiotropy
  f1a <- coef(summary(lm(target_gene ~ gi)))
  if (all('gi' %in% rownames(f1a))) {
    f1a <- f1a['gi', 't value']
  } else {
    f1a <- NA
  }
  
  
  # Scenario 1B: index gene causally affects
  # Corrected for LD/pleiotropy
  f1b <- coef(summary(lm(target_gene ~ gi + gi_neighbor)))
  if (all(c('gi', 'gi_neighbor') %in% rownames(f1b))) {
    f1b <- f1b['gi', 't value']
  } else {
    f1b <- NA
  }
  
  
  # Scenario 2A: neighbor GI causally affects target gene
  # Uncorrected for LD/pleiotropy
  f2a <- coef(summary(lm(target_gene ~ gi_neighbor)))
  if (all('gi_neighbor' %in% rownames(f2a))) {
    f2a <- f2a['gi_neighbor', 't value']
  } else {
    f2a <- NA
  }


  # Scenario 2B: neighbor GI causally affects
  # Corrected for LD/pleiotropy
  f2b <- coef(summary(lm(target_gene ~ gi + gi_neighbor)))
  if (all(c('gi', 'gi_neighbor') %in% rownames(f2b))) {
    f2b <- f2b['gi_neighbor', 't value']
  } else {
    f2b <- NA
  }
  
	# extract and return t-values
	return(c(f1a, f1b, f2a, f2b))
}

# load correlations between different variables in data
load(sprintf('%s/simulation_rho_gis.Rdata', out_path))
load(sprintf('%s/simulation_rho_gi_index.Rdata', out_path))
load(sprintf('%s/simulation_rho_index_target.Rdata', out_path))

# calculate summary statistics
range_rho_gis <- summary(rho_gis)[c('Min.', '1st Qu.', 'Median', '3rd Qu.', 'Max.')]
range_rho_gi_index <- summary(rho_gi_index)[c('Min.', '1st Qu.', 'Median', '3rd Qu.', 'Max.')]
range_rho_index_target <- summary(rho_index_target)[c('Min.', '1st Qu.', 'Median', '3rd Qu.', 'Max.')]

# create data frame with all different settings
settings <- expand.grid(
  gi_rho = range_rho_gis,
  gi_index_rho = range_rho_gi_index,
  index_target_rho = range_rho_index_target
)

# function to extract settings and run simulation for those settings
run <- function(i) {
  gi_rho <- settings$gi_rho[i]
  cis_rho <- settings$gi_index_rho[i]
  trans_rho <- settings$index_target_rho[i]
  return(run_simulation(gi_rho = gi_rho, cis_rho = cis_rho, trans_rho = trans_rho))  
}

# set variables for simulation
max.cores <- 40
BPPARAM <- MulticoreParam(workers = max.cores, verbose = TRUE)
n_iter <- 500 # number of times a particular simulation is run

# create list for output
tvals <- list()
length(tvals) <- n_iter

# loop over all settings n_iter times
for (ni in 1:n_iter) {
  print(ni)
  out <- bplapply(1:nrow(settings), FUN = run, BPPARAM = BPPARAM)
  tvals[[ni]] <- do.call(rbind, out)
}

# combine results
out <- do.call(rbind, tvals)
out <- as.data.frame(out)
names(out)[1:4] <- paste('tval', 1:4, sep = '')

# melt data frame for ggplot
out <- melt(out, id.vars = c('gi_rho', 'gi_index_rho', 'index_target_rho'))
names(out)[names(out) == 'value'] <- 'tvalue'
names(out)[names(out) == 'variable'] <- 'setting'
out$pvalue <- 2 * pnorm(-abs(out$tvalue))

# calculate power and mean test statistics
out_summarized <- out %>%
  group_by(setting, gi_rho, gi_index_rho, index_target_rho) %>%
  summarize(
    tvalue = mean(tvalue),
    power = mean(pvalue < .05)
  )

# plot results
p <- ggplot(out_summarized, mapping = aes(
  x = index_target_rho, y = power,
  color = as.character(setting),
  group = as.character(setting)
)) +
  geom_smooth() +
  facet_grid(as.character(gi_rho) ~ as.character(gi_index_rho)) +
  theme_bw() +
  theme(legend.position = 'n')

pdf(sprintf('%s/simulation_power.pdf', out_path), width = 10, height = 8)
plot(p)
dev.off()
