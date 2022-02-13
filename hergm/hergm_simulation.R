#installation also requires Rtools for the new version of R
#library(devtools)
#install_github("cran/hergm")
library(hergm)
set.seed(0)
d <- matrix(rbinom(n = 400, size = 1, prob = .1), 20, 20)#initial network
d <- as.network(d, directed = TRUE, loops = FALSE)
N <- 2

#networks with one neighbourhood
indicator <- rep.int(1, 20)
eta <- c(-1, -2, 1, 0)
edgelists <- simulate.hergm(d ~ edges_ij + transitiveties_ijk,max_number = 1, indicator = indicator, eta = eta, sample_size = N)#networks stored in a list form,edgelists$edgelist[[1]]

#networks with two neighbourhoods
indicator <- c(rep.int(1, 10), rep.int(2, 10))
eta <- c(-1, -1, -2, 1, 1, 0)
edgelists <- simulate.hergm(d ~ edges_ij + transitiveties_ijk,max_number = 2, indicator = indicator, eta = eta, sample_size = N)

d <- edgelists$edgelist[[1]]
d <- as.network(d, directed = TRUE, matrix.type = "edgelist")
object <- hergm(d ~ edges_ij+ transitiveties_ijk, max_number = 2,
                number_runs = 3, verbose = -1)

#estimation of the generated network
estimate <- function(i, edgelists) {
  d <- edgelists$edgelist[[i]]
  d <- as.network(d, directed = TRUE, matrix.type = "edgelist")
  object <- hergm(d ~ edges_ij + transitiveties_ijk, max_number = 2,
                  number_runs = 3, verbose = -1)
  indicator <- vector(length = d$gal$n)
  for (i in 1:d$gal$n) indicator[i] <- which.max(object$p_i_k[i, ])
  number <- length(unique(indicator))
  return(number)
  }
library("parallel")
RNGkind("L'Ecuyer-CMRG")
number <- mclapply(1:N, estimate, edgelists, mc.cores = 1)
number <- as.numeric(unlist(number))
barplot(table(factor(number, levels = 1:2)), ylim = c(0, N), space = 2)


d <- edgelists$edgelist[[1]]
d <- as.network(d, directed = TRUE, matrix.type = "edgelist")
object <- hergm(d ~ edges_ij + transitiveties_ijk, max_number = 2,
                number_runs = 3, verbose = -1)
#object <- hergm(kapferer ~ edges_ij, parametric = TRUE, max_number = 2, sample_size = 1e+3)

########################Application###################################
data("kapferer", package = "hergm")
set.seed(0)
object <- hergm(kapferer ~ edges_i, parametric = TRUE, max_number = 2, sample_size = 1e+3)
tmcmc.diagnostics(object)
object$mcmc.diagnostics
object


data("bali", package = "hergm")
set.seed(0)
indicator <- rep.int(1, 17)
object.m1 <- hergm(bali ~ edges_ij + triangle_ijk, max_number = 1,
                   indicator = indicator, sample_size = 1e+4)#default ml doesn't work for smaller nodes
object.m1 <- hergm(bali ~ edges_ij + triangle_ijk, max_number = 1, method="bayes",
                   indicator = indicator, sample_size = 1e+5)
set.seed(0)
indicator <- c(rep.int(1, 9), rep.int(2, 5), rep.int(1, 3))
object.m2 <- hergm(bali ~ edges_ij + triangle_ijk, max_number = 2,
                   indicator = indicator, sample_size = 2e+5)
set.seed(0)
indicator <- c(rep.int(NA, 9), rep.int(2, 5), rep.int(NA, 3))
object.m3 <- hergm(bali ~ edges_ij + triangle_ijk, max_number = 2,
                   indicator = indicator, sample_size = 2e+5)
set.seed(0)
object.m4 <- hergm(bali ~ edges_ij + triangle_ijk, max_number = 5,
                   sample_size = 2e+5, number_runs = 3)#default is ml

data("bunt", package = "hergm")
#bunt
#gplot(bunt, mode = "kamadakawai", displaylabels = FALSE)
object.local <- hergm(bunt ~ edges_ij + mutual_ij + ttriple_ijk, max_number = 32, relabel = 3, sample_size = 2e+5)#allocate to 32 groups

indicator <- rep.int(1, 32)
object.global <- hergm(bunt ~ edges_ij + mutual_ij + ttriple_ijk, max_number = 1, indicator = indicator, sample_size = 2e+5)