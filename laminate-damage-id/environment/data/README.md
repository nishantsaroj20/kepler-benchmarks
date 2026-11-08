# Finite prestressed laminate model

All records are synthetic. Units are metres, seconds, newtons, kilograms and degrees. Geometry, material constants, spring stiffnesses, bounds, load values, forecast forces and noise scales are in `plate_specs.json`. Coordinates use the lower-left origin, positive x right, positive y up the plate, and positive transverse displacement upward. The symmetric eight-ply layup is listed bottom to top; 0 degrees is the positive x direction.

## Breathing split

One interface is open to separation: the midplane between ply 4 and ply 5, counted from the bottom. It is the same interface in all three zones. Each zone has centre `(xc,yc)`, elliptical-support semi-axes `a,b` with `a>b`, counterclockwise major-axis angle `theta_deg`, and one severity `d`. Set

$ u_j=\cos\theta_j(x-x_{c,j})+\sin\theta_j(y-y_{c,j}) $,
$ v_j=-\sin\theta_j(x-x_{c,j})+\cos\theta_j(y-y_{c,j}) $,
$ g_j=\left(\max\left(0,\,1-(u_j/a_j)^2-(v_j/b_j)^2\right)\right)^2 $.

`g_j` is identically zero outside that ellipse. Angles are periodic modulo 180 degrees. `a` and `b` are support semi-axes, not e-fold radii, diameters or standard deviations. There is no further clipping and no rule that the ellipses must be disjoint. Evaluate the profiles only on the plate. The separated fraction is

$ s=1-\prod_{j=1}^{3}(1-d_j g_j) $.

This product is not a cross-ply or angle-ply knockdown. `D_intact` is the bending stiffness of the intact laminate about the plate midplane. `D_split` is the bending stiffness of the lower four plies about their own centroid plus the bending stiffness of the upper four plies about their own centroid. `A_intact` is the intact extensional stiffness. On a preparation with `q11>0` the split is open:

$ D(x,y)=D_{\mathrm{intact}}+s(x,y)\,(D_{\mathrm{split}}-D_{\mathrm{intact}}) $,

$ A(x,y)=A_{\mathrm{intact}}\,(1-\texttt{open\_membrane\_loss}\,s(x,y)) $.

On a preparation with `q11<0` the bending change is `closed_bending_retain` times the open-branch bending change, and the membrane stays intact:

$ D(x,y)=D_{\mathrm{intact}}+\texttt{closed\_bending\_retain}\,s(x,y)\,(D_{\mathrm{split}}-D_{\mathrm{intact}}) $,

$ A(x,y)=A_{\mathrm{intact}} $.

`closed_bending_retain` and `open_membrane_loss` are the numbers in `plate_specs.json`. With the published value 1, the closed and open bending fields are the same function of `s`. Build both tensors from the given unidirectional constants with the engineering shear strain `(exx, eyy, 2exy)` and reciprocal Poisson ratio `nu21=nu12*E2/E1`. `G13` and `G23` do not enter this Kirchhoff plate. Density and Poisson ratios do not change. Retain the `D16` and `D26` couplings. A bending-only map, or a closed preparation treated as fully intact, does not reproduce the holds.

## Reduced displacement space

The supplied finite forward model is defined by the following energy and spaces, not by a required solution algorithm. Let

$ \phi_{mn}=\sin(m\pi x/L_x)\sin(n\pi y/L_y),\quad m,n=1,2,3 $.

Order coefficients with m outermost and n innermost. Use all nine functions for each of the incremental fields `w=h*sum q_mn phi_mn`, `u=sum r^u_mn phi_mn` and `v=sum r^v_mn phi_mn`. Thus `q` is dimensionless, and the eighteen in-plane coefficients have units of length. The known stress-free imperfection is `w0=h*sum eta_mn phi_mn`, with eta from the specification. The total transverse position is `w0+w`. Incremental u, v and w vanish on all edges. Springs oppose the incremental normal slope relative to the stress-free shape; they do not fix that slope.

Define

$ \kappa=(-w_{,xx},-w_{,yy},-2w_{,xy})^T $,

$ e=(u_{,x}+w_{0,x}w_{,x}+w_{,x}^2/2,\ v_{,y}+w_{0,y}w_{,y}+w_{,y}^2/2,\ u_{,y}+v_{,x}+w_{0,x}w_{,y}+w_{0,y}w_{,x}+w_{,x}w_{,y})^T $.

Prescribed uniform compressive initial-stress resultants are `Nx` from the specification and `Ny=0.35 Nx`, both positive in compression. They remain prescribed during the test; redistribution of this initial-stress field is outside the reduced model. `F` is the known total transverse force, applied as uniform pressure `F/(Lx Ly)`. Its positive direction is upward. The potential is

$ \Pi(q,r;F)=\frac12\int[\kappa^T D\kappa+e^T A e-N_x((w_{0,x}+w_{,x})^2-w_{0,x}^2)-N_y((w_{0,y}+w_{,y})^2-w_{0,y}^2)]\,dA +\frac12\sum_{\mathrm{edges}}\int k_r(\partial_n w)^2\,ds-\frac{F}{L_x L_y}\int w\,dA $.

The area integral in this finite model means the 30 by 26 tensor Gauss-Legendre rule. Thickness integration and the sine-function boundary integrals are exact. This nine-term shallow-plate model is the target, not an instruction to converge a three-dimensional continuum mesh. The x=0 and y=Ly springs have unknown positive multipliers shared by every record; the other two edges remain nominal.

For fixed `q`, the in-plane energy has a unique minimizer. Define `Pi_*(q;F)=min_r Pi(q,r;F)`. Equilibria satisfy `grad_q Pi_*=0` and are stable when `H=grad_q^2 Pi_*` is positive definite. In-plane inertia, rotary inertia and damping are omitted. The dimensionless-coordinate mass is `M_q=rho*h*Lx*Ly*h^2*I/4`. Measured infinitesimal frequencies and shapes are eigenpairs of `H z = omega^2 M_q z` about the prepared equilibrium, with `f=omega/(2 pi)` and transverse sensor samples of `sum z_mn phi_mn`. Relaxation of `r` is part of `H`, not a frozen-in-plane approximation. `D` in the potential is the open or closed field of the preparation being traced, not a single field shared by both preparations.

## Preparation, observations and forecast

At F=0 each preparation starts at the minimum of the potential in one of the two half-spaces q11<0 or q11>0. The chosen state has nonzero q11 and is continued quasistatically to the specified held force without crossing a stability loss. These two preparations define the negative and positive branches. Each record comes from one preparation, but its branch label was not recorded. Only alternating-displacement modal shapes are observed: no static deflection measurements are supplied. A mode-shape column's arbitrary sign does not identify the equilibrium's branch.

Each record directory contains `frequencies.csv` (`resonance_id,frequency_Hz`) and `mode_shapes.csv` (`sensor_id,resonance_1,...,resonance_4`). These are four distinct modes among the lowest seven of the stable relaxed tangent. Omission and order may differ across records. A resonance ID associates a frequency with its shape column within that record only. Sensor positions are in `sensors.csv`; sensor IDs associate shape rows with positions. Each clean seven-sensor shape is normalized before component noise and a nonzero signed column scale are applied. Frequency noise is relative. The supplied noise law describes the generator, not an acceptance residual threshold.

In the output, a branch sign is the sign of the prepared equilibrium's q11. For `snap_load_N`, start on the negative preparation at F=0, where the split is closed, and increase F. `snap_load_N` is the first positive F at which this same stable branch loses stability: the smallest relaxed tangent eigenvalue reaches zero. It is not the force where the two minima have equal energies, a jump force under a finite loading rate, a flat-state buckling load, or a globally re-minimized energy path. All supplied held records precede their preparation's stability loss.

`open_branch_forecast_Hz` is a different calculation. At each force in `forecast.forces_N`, continue the positive preparation, with the split open, and report its lowest four natural frequencies in ascending order. Neither forecast force is one of the held records. The forecast uses the same finite model as the records.
