# fea-contact-shift

Expert time estimate: 14 hours.

## Difficulty

This task mirrors a paid contact-mechanics workflow: recover end-of-history contact diagnostics from mesh-and-load packs for plane-strain frictional contact. The packs are synthetic, seeded, and continuum-faithful rather than raw lab files. Each pack supplies an interface gap from one or more parabolic patches, elastic constants, a Coulomb coefficient, and a discrete $ F_n,F_t $ history.

Frontier failures are scientific. Neighboring patches interact through the elastic log kernel, so treating each patch with an isolated Hertz formula misplaces peak pressure and contact length once clearances and spacings shift. Frictionless or tied-contact reductions match $ F_n $ but miss tangential reaction and mobilization under partial slip. Non-proportional and shear-reversal histories require an incremental stick-slip state; a final-step Coulomb projection without history drifts on reverse shear. Agents that overfit the four practice packs fail held-out packs whose patch layout, moduli, $ \mu $, and load-path family are deliberately shifted rather than reseeded.

Pointwise stick maps under pure force control are not unique for Coulomb friction; the graded observables are the unique integral and contact-pressure diagnostics named in the instruction.

## Reference solution

The reference builds the piecewise-constant plane-strain influence matrix from the integrated log kernel with the pack’s `rref`, then marches each `loads.csv` row with an active-set complementarity solve. Normal unknowns satisfy Signorini conditions with a rigid approach chosen to meet $ F_n $. Tangential unknowns use an accumulated slip state so stick nodes enforce incremental tangential compatibility while slip nodes sit on $ |q|=\mu p $ with consistent sign. After the final step it reports peak and RMS contact pressure, total reactions, contact length, pressure centroid, and friction mobilization.

## Verification

Held-out grading materializes eight packs inside the verifier, seeded from the bytes of `/app/solver/analyze.py`, using shifted geometry, stiffness, friction, and load-path families. Each run executes the submitted program with an 8 s wall-clock cap (about three times the reference’s worst practice pack on CPU) and `OMP`/`OPENBLAS`/`MKL` thread pins equal to 1. Outputs must match an independent truth solver sealed under `tests/`: same continuum conventions, separately coded active-set loop, interface upsampled by two, no import of the reference analyzer.

Tolerances follow native-versus-refined mesh spread on correct solves: relative 0.12 on `p_peak`, `p_rms`, and `contact_length`; 0.02 on `fn_total`; 0.08 on `ft_total`; absolute 0.05 on `mobilization`; absolute $ 5\times10^{-4} $ m on `x_p_centroid`. Zero, constant, and mean-like JSON baselines fail every held-out pack. Offline local validation without Harbor is recorded under `authoring/evidence/`.
