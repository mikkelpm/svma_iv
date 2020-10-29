# Instrumental variable identification of dynamic variance decompositions

Matlab code for inference on variance decompositions and the degree of invertibility/recoverability in a general Structural Vector Moving Average (SVMA) model identified by external instruments (IVs)

**Reference:**
Plagborg-Møller, Mikkel and Christian K. Wolf (2020), "Instrumental Variable Identification of Dynamic Variance Decompositions", https://scholar.princeton.edu/mikkelpm/decomp_iv

Tested in: Matlab R2020a on Windows 10 PC (64-bit)

## Contents

**[functions](functions):** Matlab routines
- [SVMAIV_estim.m](functions/SVMAIV_estim.m): main function for SVMA-IV inference
- [SVARIV_estim.m](functions/SVARIV_estim.m): SVAR-IV inference (assumes invertibility)

**[application](application):** empirical example
- [run_gk.m](application/run_gk.m): example based on Gertler & Karadi (2015)

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
- `inv_test` structure: pre-test for invertibility, implemented either as a Granger casuality test for all equations jointly (subfield `all`) or in each *y* equation separately (subfield `eqns`)
  - `inv_test.wald_stat`: Wald statistics
  - `inv_test.df`: degrees of freedom
  - `inv_test.pval`: p-values

Parameter names:
- `alpha`: scale parameter
- `R2_inv`: degree of invertibility
- `R2_recov`: degree of recoverability
- `FVR`: forecast variance ratio
- `FVD`: forecast variance decomposition

See the [empirical application](application/run_gk.m) for a concrete example. Additional optional arguments to [`SMVAIV_estim.m`](functions/SVMAIV_estim.m) are listed at the top of the function.
