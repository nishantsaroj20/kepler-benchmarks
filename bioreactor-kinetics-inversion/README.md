# bioreactor-kinetics-inversion

Fed-batch kinetic parameter estimation under a metabolic phase shift, followed by mass-balance-constrained dynamic feed optimization.

## Difficulty

Synthetic telemetry is generated from a stiff five-state fed-batch model with Haldane–Andrews kinetics, product inhibition, and a discontinuous metabolic switch, then sparsified with irregular sampling and sensor noise (including mildly noisy feed). Fermentation scientists and bioprocess engineers encounter the same workflow when characterizing strains and designing feeding policies. Hardness comes from jointly identifying non-smooth kinetics from limited noisy assays and solving a constrained dynamic optimization whose feasible set is shaped by substrate and dissolved-oxygen limits—not from withheld textbook facts.

## Reference solution

Fit the seven unknown kinetic parameters by simulating the unstructured mass-balance ODEs against both historical runs, using a coarse grid over weakly identifiable coordinates (`t_switch`, `P_crit`) followed by multi-start nonlinear least squares while holding measured feed as a noisy zero-order-hold input and reading known `gamma` from reactor specs. Then optimize a ten-interval piecewise-constant feed with SLSQP under volume, substrate, and dissolved-oxygen penalties to maximize final product mass $V(t_f)P(t_f)$.

## Verification

Ground truth parameters and a reference feed were produced by independent Radau integration with a fixed seed in `authoring/provenance/`. The verifier re-integrates the submitted feed under the true kinetics and checks relative parameter tolerances (10% on `mu_max`/`P_crit`/`alpha`/`beta`, 20% on `K_S`/`K_I`, `t_switch` within 0.5 h), $S\le 15$ g/L (0.5 g/L slack), $C_L\ge 0.15 C^\star$, $V\le V_{\max}$, and minimum product mass 249.33 g. An independent reimplementation recorded in `authoring/evidence/` recovers oracle reward 1 and nop reward 0 with matching constraint margins.
