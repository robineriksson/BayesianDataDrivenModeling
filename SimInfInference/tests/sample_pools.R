## Robin Eriksson 2018
library(SimInfInference)

Svec <- sample(10:20, 5, TRUE)
Ivec <- sample(1:5, 5, TRUE)
mat <- t(matrix(c(Svec, Ivec), nrow = 5))

## does seeding work?
set.seed(1)
test1 <- sample_nodes_Rcpp(mat, get_seed())
print(test1)

set.seed(1)
test2 <- sample_nodes_Rcpp(mat, get_seed())
print(test2)

set.seed(2)
test3 <- sample_nodes_Rcpp(mat, get_seed())
print(test3)
## same seed -> same result
stopifnot(all(test1 == test2))

## different seed -> different result
stopifnot(!all(test1 == test3))
