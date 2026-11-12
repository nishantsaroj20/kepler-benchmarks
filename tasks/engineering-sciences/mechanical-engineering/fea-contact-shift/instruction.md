Write a reusable contact-diagnostics program at `/app/solver/analyze.py` that, given a mesh-and-approach pack for plane-strain frictional contact against an elastic half-space, recovers a fixed set of contact-state metrics at the end of the pack’s prescribed approach history.

The program must be invoked as `python3 /app/solver/analyze.py --pack <pack_dir> --out <out_json>`, where `<pack_dir>` and `<out_json>` are absolute paths. It must exit with status 0 and write a single JSON object to `<out_json>`. Anything on stderr is ignored. The program must work on every pack directory under `/app/data/practice/` and on further held-out packs that use the same file schema but deliberately shifted geometry, friction, stiffness, and approach paths. Those held-out packs are not present in the image you see while developing.

Each pack directory contains exactly these files:

- `meta.json` with finite numbers `E` (Pa), `nu` (dimensionless), `mu` (Coulomb friction coefficient), and `rref` (m, reference length in the kernel). The boolean `plane_strain` is true, `kernel` is the string `plane_strain_log`, and `control` is the string `approach`.
- `interface.csv` with header `x,gap0`. Column `x` is the interface sample coordinate in metres. Column `gap0` is the initial gap in metres (positive when open). Samples are ordered by increasing `x`.
- `loads.csv` with header `step,delta_n,delta_t`. Rows are the discrete rigid-punch history in order. `delta_n` is the prescribed normal approach in metres and `delta_t` is the prescribed tangential shift in metres. Apply the rows sequentially exactly as listed. Total forces are not prescribed; they are part of the solution.

Use the following continuum conventions. The bodies are reduced to a rigid punch against a homogeneous isotropic linear-elastic half-space in plane strain with Young modulus `E` and Poisson ratio `nu`. Identical-material Dundurs coupling is neglected, so the normal and tangential surface operators share the same influence kernel and interact only through the Coulomb conditions. Piecewise-constant normal pressure $ p $ and tangential traction $ q $ live on the intervals implied by midpoints between consecutive `x` samples (end intervals mirrored). With $ A = 2(1-\nu^2)/(\pi E) $ and reference length `rref`, the elastic surface fields obey the integrated log kernel

$ u_n(x_i) = A \sum_j K_{ij} p_j $, $ u_t(x_i) = A \sum_j K_{ij} q_j $,

where $ K_{ij} = \int_{L_j}^{R_j} \ln(|x_i-s|/\mathrm{rref})\,ds $ over element $ j $ with ends $ L_j,R_j $. After a prescribed approach $ \delta_n $, the gap at sample $ i $ is $ g_i = \mathrm{gap0}_i - \delta_n - u_n(x_i) $. Contact requires $ p_i \ge 0 $, $ g_i \ge 0 $, and $ p_i g_i = 0 $. Coulomb friction requires $ |q_i| \le \mu p_i $. Carry an accumulated tangential slip state $ \gamma_i $ across steps. On stick, $ u_t(x_i) + \gamma_i = \delta_t $. On slip, $ q_i $ sits on the friction bound with sign fixed by the incremental slip direction and $ \gamma_i $ is updated so that the slip decomposition remains consistent. If every sample stays open after a step, set all tractions and $ \gamma_i $ to zero for that step.

After the final history row, write JSON with exactly these keys and no others:

- `schema_version`: integer `1` (a JSON number that is an integer, not a boolean)
- `p_peak`: maximum $ p_i $ (Pa); `0.0` if there is no contact
- `fn_total`: $ \sum_i p_i w_i $ (N/m)
- `ft_total`: $ \sum_i q_i w_i $ (N/m)
- `contact_length`: $ \sum_i w_i $ over samples with $ p_i > 10^{-8}\max_j p_j $ (m); `0.0` if there is no contact
- `q_peak`: maximum of $ |q_i| $ (Pa); `0.0` if there is no contact
- `stick_fraction`: contact-length fraction of contacting samples that are not on the friction bound $ |q_i| \ge 0.99\,\mu p_i $; `0.0` if there is no contact
- `slip_length`: contact length on the friction bound as above (m); `0.0` if there is no contact
- `x_q_centroid`: $ |q| $-and-width-weighted mean of $ x $ over contact samples (m); `0.0` if there is no contact or zero tangential traction
- `cum_slip_l1`: $ \sum_i w_i |\gamma_i| $ after the final step (m²)

Here $ w_i = R_i - L_i $. All numeric values must be finite JSON numbers. Do not include NaN or Infinity. Do not include additional keys.

The same program is re-run later on eight held-out packs drawn from the same modeling family with shifted patch layout, radii, clearances, $ E $, $ \nu $, $ \mu $, and approach-path shape (normal-then-shear, proportional, shear-reversal, and no-contact histories). Practice packs are not a reseed of those cases. Each program invocation receives at most 8 seconds of wall-clock time, 1 CPU thread, 2048 MB address space, 64 processes, and empty stdin, with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `MKL_NUM_THREADS` pinned to 1. The process runs unprivileged. Only `/app/solver/analyze.py`, the standard library, NumPy 1.26.4, SciPy 1.14.1, and the pack directory are available in that process. Working directory starts empty under `/tmp`. Network is unavailable during those runs.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
