# Instrumental variable identification of dynamic variance decompositions

Matlab code for inference on variance decompositions and the degree of invertibility/recoverability in a general Structural Vector Moving Average (SVMA) model identified by external instruments (IVs, also known as proxies)

**Reference:**
Plagborg-Møller, Mikkel and Christian K. Wolf (2021), "Instrumental Variable Identification of Dynamic Variance Decompositions", https://scholar.princeton.edu/mikkelpm/decomp_iv (paper + online appendix)

Tested in: Matlab R2020a on Windows 10 PC (64-bit)

Other versions: [Python code](https://github.com/jbduarte/svma_iv) written by João B. Duarte

## Contents

**[functions](functions):** Matlab routines
- [SVMAIV_estim.m](functions/SVMAIV_estim.m): main function for SVMA-IV inference
- [SVARIV_estim.m](functions/SVARIV_estim.m): SVAR-IV inference (assumes invertibility)

**[applications](applications):** empirical applications
- [gk/run_gk.m](applications/gk/run_gk.m): application based on Gertler & Karadi (2015)
- [kaenzig/run_kaenzig.m](applications/kaenzig/run_kaenzig.m): application based on Känzig (2021)

**[illustration](illustration):** numerical illustration
- [run_sw.m](illustration/run_sw.m): SVMA-IV and SVAR-IV analysis of Smets & Wouters (2007) model

**[simulations](simulations):** simulation study
- [run_sim.m](simulations/run_sim.m): run Monte Carlo experiments for various DGPs
- [print_results.m](simulations/print_results.m): display results for all DGPs

## Example

``` Matlab
addpath('functions');
% Given (T x n_y) data matrix Y with endogenous variables
% and (T x 1) data vector Z with external instrument
[bounds, id_recov, inv_test] = ...
  SVMAIV_estim(Y, Z, ...
               'ic', 'aic', ...   % Select lag length using AIC
               'signif', 0.1, ... % 10% significance level
               'n_boot', 500, ... % 500 bootstrap iterations
               'horiz', 1:24);    % Compute horizons 1-24 of FVR/FVD
```
Output:
- `bounds` structure: partial identification bounds
  - `bounds.estim`: point estimates of bounds
  - `bounds.ci`: confidence intervals for identified set
- `id_recov` structure: results under additional assumption of recoverability
  - `id_recov.estim`: point estimates of parameters
  - `id_recov.ci`: confidence intervals for parameters
- `inv_test` structure: pre-test for invertibility, implemented as a Granger casuality test, either in all equations jointly (subfield `all`) or in each *y* equation separately (subfield `eqns`)
  - `inv_test.wald_stat`: Wald statistics
  - `inv_test.df`: degrees of freedom
  - `inv_test.pval`: p-values

Parameter names:
- `alpha`: scale parameter
- `R2_inv`: degree of invertibility
- `R2_recov`: degree of recoverability
- `FVR`: forecast variance ratio
- `FVD`: forecast variance decomposition

See our [empirical applications](applications) for concrete examples of how to use the code. Additional optional arguments to [`SMVAIV_estim.m`](functions/SVMAIV_estim.m) are listed at the top of the function.

## Replication instructions

The figures and tables in our paper and supplement are produced as follows:
- Figures 1 and B.7 and Table 1: [applications/gk/run_gk.m](applications/gk/run_gk.m)
- Table 2: First run [simulations/run_sim.m](simulations/run_sim.m) nine separate times, with the variable `model.dgp` at the top set to each of the options `0` through `8`. Then run [simulations/print_results.m](simulations/print_results.m)
- Figures B.1-B.6 and Table B.1: [illustration/run_sw.m](illustration/run_sw.m)
- Figures B.8-B.11: [applications/kaenzig/run_kaenzig.m](applications/kaenzig/run_kaenzig.m)

## Acknowledgements
We are grateful to [Diego Känzig](https://github.com/dkaenzig) for allowing us to reproduce some of the data files used in his [2021 AER paper](https://www.aeaweb.org/articles?id=10.1257/aer.20190964).

Wolf acknowledges support from the Alfred P. Sloan Foundation and the Macro Financial Modeling Project. Plagborg-Møller acknowledges that this material is based upon work supported by the National Science Foundation under Grant \#1851665.
