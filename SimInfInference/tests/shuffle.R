## Robin Eriksson 2018
library(SimInfInference)
a <- 1:10
set.seed(1)
b <- testShuffle(a)

## print(a)
## print(b)
## does the shuffle work?
stopifnot(!all(a == b))

## does the seed work?
set.seed(2)
b2 <- testShuffle(a)
##print(b2)
stopifnot(!all(b == b2))
