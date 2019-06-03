# Bayesian Data-driven epidemiological modeling
### See the paper:
Bayesian epidemiological modeling over high-resolution network data: opportunities for optimized control.


## The release code version

### About
This collection of computational scripts were used for the results in
the paper: Bayesian data-driven epidemiological modeling.

The code is built for the purpose of i) to infer the parameters of the
epidemiological model given network- (transport) and surveillance-
data, and ii) to use the posterior for a Bayesian framework in a public
health setting, exploring detection and mitigation strategies.

In the code structure, there is an R-package which acts as a wrapper
for main simulator, SimInf. The package is written as a class (R6)
which can be called from scripts and it connects the three main parts
of the inference: a) the simulator, b) the summary statistics, and c)
the Bayesian sampler (estimator).

### A summary of the quantitative result
We are able to with the code and data in this repository to first
approximate the posterior of the STEC endemic in the cattle population

If we visually define our problem, we have a (a) local $\text{SIS}_E$
spread of the infection withing herds. The herds are (b) connected via
a a national transport network with other herds. In the timeframe of
bacterial sampling, the herds are (c) sampled for the bacteria were 1
indicates that the herd was found infected while 0 is a negative
result.  ![Visualization of the modeling problem](/Figures/introSingle.png)

On the network, we have observational data, and using SLAM (an
approximate Bayesian computation) approach to approximate the
parameter posterior (green). We then later affirm the found posterior
by re-iterating the same procedure with synthetic data generated from
the parameter posterior mean (blue & orange)
![Approximatedposteriordistribution](/Figures/realMultiPosterior.png)

Using this posterior we can simulate different scenarios, and easily
carry a Bayesian credible interval with us. First we explore different
detection procedures, were we try and find the best possible grouping
when deciding on sentinel node groups. We plot the probability of
detection mean and 95% CI for each node selection group against time.
![Bayesian modeling of detection experiments](/Figures/detection.png)

Latter we study some intervention procedures, and see how the
population prevalence decay with time. Here we also plot the mean and
95% CI for the different techniques considered.  ![Bayesian modeling
of intervention experiments](/Figures/intervention.png)


### Structure
The catalog structure follows,

* BPD
  - Scripts
    - DATA
    - INFERENCE
      - 1600system
      - realsystem
    - PREVDEC
  - README.md
  - SimInfInference
    - R
    - src
  - Figures

#### Scripts
All the scripts used for the computations, using the class from SimInfInference.

* **DATA**, The collection of data used in the paper.

* **INFERENCE**, Divided into three subfolders, the 1600 synthetic data
  (1600system) and the full dataset (realsystem).

* **PREVDEC**, The experiments that use the posterior generated from
  INFERENCE in the frame of public health.

#### SimInfInference Class inspired code structure presented as a
R-package. Acts as a wrapper around SimInf.

4 main compartments are used:
* Simulator: what simulator to use, we
  wrap around SimInf, but can be easily extended to other alternatives.
* SummaryStatistics: how to summarize the data.
* Proposal: how to propose new parameters.
* Estimator: what parameter estimator to use.

The code is written in two folders, R and src. Where R hold the R code
and src the RCPP. RCPP was used for some task for easy speed-up.

## Run the code (Replicate the result)
First install the SimInfInference package by in folder running `make rcpp` followed by `make install`

The inference for each dataset is perfomed by using the assigned script.
If one wish to extract the multi-set matrix plots, or tables, one also need
to save each file. The suggested location is /Scripts/DATA/posterior/.
In said directory, the results presented in the paper is stored and available to explore.

To further replicate our results, we suggest using a computational cluster, and
creating multiple (parallel) "Markov chains" and using these together for better
performance.

## References

* The Baysian inference methods

[1] Sisson SA, Fan Y, Beaumont M (2018) Handbook of Approximate
Bayesian Computation.  (Chapman and Hall/CRC).

[2] Wood SN (2010) Statistical inference for noisy nonlinear
ecological dynamic systems. Nature 466(7310):1102–1104.

[3] Haario H, Saksman E, Tamminen J, , et al. (2001) An adaptive
Metropolis algorithm. Bernoulli 7(2):223–242.

* The simulator, epidemiological models

[4] Widgren S, Bauer P, Eriksson R, Engblom S (2018) SimInf: An R
package for data-driven stochastic disease spread
simulations. Accepted for publication in J. Stat. Softw.

* The work is influensed and succeeds the series of previously published works

[5] Bauer P, Engblom S, Widgren S (2016) Fast event-based
epidemiological simulations on national scales. Int. J. High
Perf. Comput. Appl. 30(4):438–453.

[6] Widgren S, et al. (2016) Data-driven network modelling of disease
transmission using com- plete population movement data: spread of VTEC
O157 in Swedish cattle. Veterinary Res.  47(1):81.

[7] Widgren S, Engblom S, Emanuelson U, Lindberg A (2018)
Spatio-temporal modelling of vero- toxigenic E. coli O157 in cattle in
Sweden: Exploring options for control. Veterinary Res.  49(78).
