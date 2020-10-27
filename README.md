# Instrumental variable identification of dynamic variance decompositions

Matlab code for inference on variance decompositions and the degree of invertibility/recoverability in a general Structural Vector Moving Average (SVMA) model identified by external instruments (IVs)

**Reference:**
Plagborg-Møller, Mikkel and Christian K. Wolf (2020), "Instrumental Variable Identification of Dynamic Variance Decompositions", https://scholar.princeton.edu/mikkelpm/decomp_iv

Tested in: Matlab R2020a on Windows 10 PC (64-bit)

## Contents

**[functions](functions):** Matlab routines
- [SVMAIV_estim.m](functions/SVMAIV_estim.m): main function for SVMA-IV inference
- [SVARIV_estim.m](functions/SVARIV_estim.m): SVAR-IV inference (assumes invertibility)

**[application](application):** Empirical example
- [run_gk.m](application/run_gk.m): example based on Gertler & Karadi (2015)

**[illustration](illustration):** Numerical illustration
- [run_sw.m](illustration/run_sw.m): SVMA-IV and SVAR-IV analysis of Smets & Wouters (2007) model

**[simulations](simulations):** Simulation study

