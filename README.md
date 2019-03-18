# Bayesian Data-driven epidemiological modeling
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

### Structure
The catalog structure follows,

* BPD
  - Scripts
    - DATA
    - INFERENCE
      - 1600system
      - realsystem
      - plotNtable
    - PREVDEC
  - README.md
  - SimInfInference
    - R
    - src

#### Scripts
All the scripts used for the computations, using the class from SimInfInference.

* **DATA**, The collection of data used in the paper.

* **INFERENCE**, Divided into three subfolders, the 1600 synthetic data
  (1600system), the full dataset (realsystem), and scripts for our plotting and
  table construction.

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
