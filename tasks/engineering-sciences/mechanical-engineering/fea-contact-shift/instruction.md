Write a reusable contact-diagnostics program at `/app/solver/analyze.py` that, given a mesh-and-load pack for plane-strain frictional contact against an elastic half-space, recovers a fixed set of contact-state metrics at the end of the pack’s load history.

The program must be invoked as `python3 /app/solver/analyze.py --pack <pack_dir> --out <out_json>`, where `<pack_dir>` and `<out_json>` are absolute paths. It must exit with status 0 and write a single JSON object to `<out_json>`. Anything on stderr is ignored. The program must work on every pack directory under `/app/data/practice/` and on further held-out packs that use the same file schema but deliberately shifted geometry, friction, stiffness, and load paths. Those held-out packs are not present in the image you see while developing.

Each pack directory contains exactly these files:

- `meta.json` with finite numbers `E` (Pa), `nu` (dimensionless), `mu` (Coulomb friction coefficient), and `rref` (m, reference length in the kernel). The boolean `plane_strain` is true and `kernel` is the string `plane_strain_log`.
- `interface.csv` with header `x,gap0`. Column `x` is the interface sample coordinate in metres. Column `gap0` is the initial gap in metres (positive when open). Samples are ordered by increasing `x`.
- `loads.csv` with header `step,Fn,Ft`. Rows are the discrete load history in order. `Fn` is the prescribed total normal force per unit depth in N/m and `Ft` is the prescribed total tangential force per unit depth in N/m. Apply the rows sequentially exactly as listed.

Use the following continuum conventions. The bodies are reduced to a rigid punch against a homogeneous isotropic linear-elastic half-space in plane strain with Young modulus `E` and Poisson ratio `nu`. Identical-material Dundurs coupling is neglected, so the normal and tangential surface operators share the same influence kernel and interact only through the Coulomb conditions. Piecewise-constant normal pressure $ p $ and tangential traction $ q $ live on the intervals implied by midpoints between consecutive `x` samples (end intervals mirrored). With $ A = 2(1-\nu^2)/(\pi E) $ and reference length `rref`, the elastic surface fields obey the integrated log kernel

$ u_n(x_i) = A \sum_j K_{ij} p_j $, $ u_t(x_i) = A \sum_j K_{ij} q_j $,

where $ K_{ij} = \int_{L_j}^{R_j} \ln(|x_i-s|/\mathrm{rref})\,ds $ over element $ j $ with ends $ L_j,R_j $. The gap at sample $ i $ after a rigid normal approach $ c_n $ is $ g_i = \mathrm{gap0}_i - c_n - u_n(x_i) $. Contact requires $ p_i \ge 0 $, $ g_i \ge 0 $, and $ p_i g_i = 0 $. Coulomb friction requires $ |q_i| \le \mu p_i $. On stick, the incremental tangential slip vanishes relative to the accumulated slip state carried from previous load steps; on slip, $ q_i $ sits on the friction bound with sign fixed by the incremental slip direction. The rigid values $ c_n $ and $ c_t $ are unknowns at each step chosen so that $ \sum_i p_i w_i = F_n $ and $ \sum_i q_i w_i = F_t $, where $ w_i = R_i-L_i $.

After the final load row, write JSON with exactly these keys and no others required for correctness beyond them:

- `schema_version`: integer `1`
- `p_peak`: maximum $ p_i $ (Pa)
- `p_rms`: contact-length-weighted root-mean-square of $ p $ over samples with $ p_i > 10^{-8}\max_j p_j $ (Pa)
- `fn_total`: $ \sum_i p_i w_i $ (N/m)
- `ft_total`: $ \sum_i q_i w_i $ (N/m)
- `contact_length`: $ \sum_i w_i $ over those contact samples (m)
- `x_p_centroid`: pressure-and-width-weighted mean of $ x $ over contact samples (m); use `0.0` if there is no contact
- `mobilization`: $ |\mathrm{ft\_total}| / (\mu\,|\mathrm{fn\_total}|) $ when $ \mathrm{fn\_total} \ne 0 $ and $ \mu > 0 $, otherwise `0.0`

All numeric values must be finite. Do not include NaN or Infinity.

The same program binary is re-run later on eight held-out packs drawn from the same modeling family with shifted patch layout, radii, clearances, $ E $, $ \nu $, $ \mu $, and load-path shape (normal-then-shear, proportional, and shear-reversal histories). Practice packs are not a reseed of those cases. Each program invocation receives at most 8 seconds of wall-clock time, 1 CPU thread, 2048 MB address space, 64 processes, and empty stdin, with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `MKL_NUM_THREADS` pinned to 1. Only `/app/solver/analyze.py`, the standard library, NumPy 1.26.4, SciPy 1.14.1, and the pack directory are available in that process. Working directory starts empty under `/tmp`. Network is unavailable during those runs.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
