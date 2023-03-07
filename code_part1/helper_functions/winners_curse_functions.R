# functions related to correcting for Winner's Curse

# function that corrects for Winner's Curse
# adapted and modified from https://github.com/cpalmer718/gwas-winners-curse
correct_winners_curse <- function (data, discovery.threshold = 5e-08) 
{
  #stopifnot(length(trait.mean) == 1)
  #stopifnot(is.numeric(trait.mean) || is.na(trait.mean))
  stopifnot(length(discovery.threshold) == 1)
  stopifnot(is.numeric(discovery.threshold) || is.na(discovery.threshold))
  #data <- read.table(input.file, header = header, sep = sep, stringsAsFactors = FALSE)
  stopifnot(nrow(data) > 0)
  #stopifnot(ncol(data) %in% c(8, 10))
  if (!(all(c("trait.mean","discovery.threshold") %in% colnames(data)))) {
    #data$trait.mean <- trait.mean
    data$discovery.threshold <- discovery.threshold
  }
  cols <- c("discovery.beta", "discovery.se",  "discovery.n", "discovery.freq",
            "discovery.threshold")
  stopifnot(all(colnames(data) %in% cols))
  
  data$debiased.beta.mle <- sapply(seq_len(nrow(data)), function(i) {
    debias.beta(data[i, "discovery.beta"], data[i, "discovery.se"], 
                data[i, "discovery.freq"], 
                data[i, "discovery.threshold"])
  })
  data$l95.mle <- calculate.ci(data[, "debiased.beta.mle"], 
                               data[, "discovery.se"], 0.025)
  data$u95.mle <- calculate.ci(data[, "debiased.beta.mle"], 
                               data[, "discovery.se"], 0.975)
  data$debiased.beta.mse <- compute.beta.mse(data[, "discovery.beta"], 
                                             data[, "discovery.se"], data[, "debiased.beta.mle"])
  data$l95.mse <- compute.ci.mse(data[, "discovery.beta"], 
                                 data[, "discovery.se"], data[, "l95.mle"], 0.025)
  data$u95.mse <- compute.ci.mse(data[, "discovery.beta"], 
                                 data[, "discovery.se"], data[, "u95.mle"], 0.975)
  for (colname in c("debiased.beta.mle", "l95.mle", "u95.mle", 
                    "debiased.beta.mse", "l95.mse", "u95.mse")) {
    data[, colname] <- signif(data[, colname], 5)
  }
  # write.table(data, output.file, row.names = FALSE, col.names = TRUE, 
  #             quote = FALSE, sep = ifelse(grepl("\\.csv$", output.file), 
  #                                         ",", "\t"))
  return(data)
}
# function pulled from https://github.com/cpalmer718/gwas-winners-curse 
debias.beta <- function(beta.biased,
                        stderr.biased,
                        freq,
                        p.threshold.init) {
  beta.nosign <- abs(beta.biased)
  p.actual <- 2.0 * (1.0 - pnorm(beta.nosign / stderr.biased))
  p.threshold <- 0
  ## deal with the fact that paper-reported thresholds for inclusion in replication
  ## are often greatly exaggerated, and if included variants do not in fact pass the threshold,
  ## the calculation needs adjustment to not explode.
  if (p.actual > p.threshold.init) {
    p.threshold <- p.actual
    warning("pair ", beta.biased, " (", stderr.biased, ") corresponds to p = ",
            p.actual, " > ", p.threshold.init, ", adjusting accordingly",
            sep = ""
    )
  } else {
    p.threshold <- p.threshold.init
  }
  
  x.lo <- -beta.nosign / 4
  x.hi <- beta.nosign
  
  f.solver <- function(x) {
    debiasing.func(x, p.threshold, beta.nosign, stderr.biased)
  }
  beta.mle <- uniroot(f.solver, c(x.lo, x.hi))$root * sign(beta.biased)
  beta.mle
}
##
debiasing.func <- function(beta.debiased.init, p.thresh,
                           beta.biased, stderr.biased) {
  ## for linear regression, the standard error is invariant with respect to
  ## regression coefficient. so skip the entire concept of standard error
  ## rescaling until other GLM linker support is attempted.
  beta.debiased <- ifelse(abs(beta.debiased.init) < 1e-16, sign(beta.debiased.init) * 1e-16, beta.debiased.init)
  adjusted.stderr <- stderr.biased
  q <- beta.debiased / adjusted.stderr - qnorm(1.0 - p.thresh / 2.0)
  r <- -beta.debiased / adjusted.stderr - qnorm(1.0 - p.thresh / 2.0)
  beta.debiased + adjusted.stderr * (dnorm(q) - dnorm(r)) / (pnorm(q) + pnorm(r)) - beta.biased
}
###
calculate.ci <- function(beta.debiased, stderr.debiased, p.thresh) {
  beta.debiased + qnorm(p.thresh) * stderr.debiased
}
###
compute.beta.mse <- function(biased.beta,
                             biased.se,
                             debiased.beta) {
  k <- biased.se^2 / (biased.se^2 + (biased.beta - debiased.beta)^2)
  k * biased.beta + (1 - k) * debiased.beta
}
###
compute.ci.mse <- function(biased.beta,
                           biased.se,
                           existing.ci,
                           p.thresh) {
  biased.ci <- calculate.ci(biased.beta, biased.se, p.thresh)
  k <- biased.se^2 / (biased.se^2 + (biased.ci - existing.ci)^2)
  k * biased.ci + (1 - k) * existing.ci
}