# fea-contact-shift

Expert time estimate: 14 hours.

## Difficulty

This task mirrors a paid contact-mechanics workflow: recover end-of-history contact diagnostics from synthetic, continuum-faithful mesh-and-approach packs for plane-strain frictional contact. Each pack supplies a multi-patch initial gap, elastic constants, a Coulomb coefficient, and a prescribed rigid-punch history $ (\delta_n,\delta_t) $. Total forces are unknowns.

Difficulty is scientific. Because loading is approach-controlled, $ F_n $ and $ F_t $ cannot be copied from the load file; they come from the contact solution. Tangential metrics ($ q $ peak, stick/slip lengths, accumulated slip $ \gamma $) depend on the incremental Coulomb history, so zeroing $ \delta_t $ during the solve or ignoring friction fails held-out shear and reversal packs even when normal pressures look plausible. Neighboring patches interact through the elastic log kernel, so isolated Hertz-per-patch shortcuts misplace peaks and contact length once clearances and spacings shift. No-contact packs require the open-gap branch.

## Reference solution

The reference builds the piecewise-constant plane-strain influence matrix from the integrated log kernel with the pack’s `rref`, then marches each `loads.csv` row with an active-set complementarity solve under prescribed $ \delta_n $ and $ \delta_t $. Normal unknowns satisfy Signorini conditions. Tangential unknowns use an accumulated slip state $ \gamma $ so stick nodes enforce $ u_t+\gamma=\delta_t $ while slip nodes sit on $ |q|=\mu p $ with consistent sign and update $ \gamma $. End-of-history diagnostics include pressure and force totals plus frictional-field metrics ($ q $ peak, stick fraction, slip length, $ |q| $-centroid, cumulative slip).

## Verification

Held-out grading materializes eight packs inside the verifier, seeded from the bytes of `/app/solver/analyze.py`, using shifted geometry, stiffness, friction, and approach-path families (including no-contact). Each run executes the submitted program as an unprivileged `runner` process with an 8 s wall-clock cap, minimal environment, `setsid` process-group cleanup, output opened with `O_NOFOLLOW`, and `/tests` sealed (`go-rwx`) so the child cannot read the truth module. Practice packs under `/app/data/practice/` are also executed.

Outputs must match an independent truth solver sealed under `tests/`: same continuum conventions, separately authored active-set loop, influence matrix from 4-point Gauss quadrature of the log kernel (not the reference antiderivative). Measured native-oracle versus Gauss-truth relative errors on practice and held-out draws stayed within about 0.02 on forces, 0.05 on peak pressure, 0.10 on contact length / $ q $ peak, 0.08 absolute on stick fraction, and 0.30 on cumulative slip; verifier tolerances sit above that spread. Zero/constant JSON baselines and frictionless ($ \delta_t=0 $) probes fail frictional metrics on contacting packs.
