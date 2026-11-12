# Fortran port of PlasticIceSheet.jl

A dependency-free Fortran library (`libplastic.a`) with the forward solver and the
closures of [PlasticIceSheet.jl](https://github.com/fesmc/PlasticIceSheet.jl), for
coupling into kryos as an ice backend. Not ported: the AD / implicit-differentiation
extension and the tau_b inversion. Added: mode `"surface"` (not in the Julia package),
which solves the plastic surface equation for z_s directly; it is the recommended mode
for kryos (see [Mode "surface"](#mode-surface)).

## API

`use plastic_ice`. All reals are `real64` (exported as `plastic_dp`); arrays are
`(nx,ny)`, `i = x`, `j = y`, as in Julia.

```fortran
type(plastic_par_t)  :: par   ! rho_i=917, rho_w=1025, g=9.80665, H_min=1, mode="full"|"flat"|"surface",
                              ! max_sweeps=200, n_outer=100, tol=1e-6, outer_tol=1e-3, relax=0.5,
                              ! warm_frac_max=0.1, margin="face"|"centre", sweep="serial"|"wavefront"
type(plastic_t)      :: ps    ! par, nx, ny, dx, dy, w=H^2, H, z_srf, n_outer_used, converged,
                              ! stat, msg, n_sweeps, outer_change, warm_used, n_new_cells

call plastic_init(ps, nx, ny, dx, dy [, par, stat])
call plastic_solve(ps, z_bed, z_sl, tau_b, mask [, warm, stat, margin_frac])
call plastic_residual(ps, z_bed, z_sl, tau_b, mask, r [, stat, margin_frac])  ! Godunov residual (w, or z_s in "surface")
call plastic_stress_residual(ps, z_bed, z_sl, tau_b, mask, rmax [, r, stat, margin_frac])
                              ! | |grad z_s| rho_i g H - tau_b | / tau_b, upwind gradient, any mode
call plastic_reset(ps)                                            ! drop the warm start
call plastic_end(ps)
H_flt = plastic_flotation_thickness(z_bed, z_sl, par)             ! elemental
V = plastic_volume(ps); A = plastic_area(ps, mask); m = plastic_surface_misfit(ps, z_srf_obs, mask)

! closures
type(plastic_glen_t)     :: rheo   ! A=1e-16 [Pa^-n yr^-1], n=3
type(plastic_sliding_t)  :: slide  ! law="linear"|"weertman"|"coulomb", beta, m, C, u_0, u_max (coulomb speed cap;
                                   !   0: none, NaN at the cap)
                                   ! fric_map(nx,ny) (optional): the law's coefficient per cell in place of beta
                                   ! (linear, weertman) or C (coulomb; needs u_max > 0)
                                   ! (steady fluxes: face = mean of its cells, margin face = its ice cell;
                                   ! plastic_velocity; not with ctrl%method = "routing")
type(plastic_smb_ctrl_t) :: ctrl   ! method="steady"|"routing";
                                   ! steady: res_tol=1e-6 [m/yr], n_steady=500, dt0=1, dt_max=1e8 [yr],
                                   !   lin="auto"|"lis"|"builtin", lin_tol=1e-8, lin_maxit=2000, verbose;
                                   ! routing: mode="surface"|"flat" (inner solve; "flat" = Julia),
                                   !   n_outer=1000, relax=0.05, tol=1.0, tau_hi=5e5, max_sweeps=400,
                                   !   sweep_tol=1e-9; both: tau_init=5e4 (cold start)
call plastic_velocity(ps, tau_b, mask, rheo, slide, ux, uy [, speed, slope_eps, u_b, u_def, stat])
call plastic_smb_implied(ps, ux, uy, mask, smb [, speed, stat])
call plastic_balance_flux(z_srf, mask, smb, dx, dy, Q [, stat])
call plastic_balance_velocity(z_srf, H, mask, smb, dx, dy, speed, ux, uy [, Q, slope_eps, H_min, stat])
call plastic_implied_basal_velocity(z_srf, H, tau_b, mask, smb, dx, dy, rheo, u_basal &
                                    [, u_balance, u_def, Q, slope_eps, H_min, stat])
call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide [, warm, ctrl, H_init, tau_d, stat, margin_frac])
call plastic_margin_frac_from_sdist(sdist, mask, dx, dy, frac)  ! margin_frac from a signed distance
call plastic_smb_flux_divergence(ps, z_bed, z_sl, mask, rheo, slide, smb [, margin_frac, stat, &
                                  q_land, q_marine, us_x, us_y, tau_cap])
                                   ! implied SMB of ps%H for method "steady" (same face fluxes); optional:
                                   ! the outflow through each ice cell's margin faces to land / to the ocean
                                   ! (the ice-free neighbour below sea level) [m ice/yr, >= 0, sum = sum(smb)],
                                   ! the surface velocity of the same closure at the cell centres [m/yr];
                                   ! tau_cap(nx,ny) [Pa]: the driving stress of the fluxes at most tau_cap
                                   ! where > 0 (a face: the mean of its cells, a margin face: its ice cell),
                                   ! e.g. tau_b (a plastic bed carries at most its yield stress)
call plastic_smb_inversion(ps, z_bed, z_sl, mask, rheo, slide, smb_prior, dHdt, smb, speed [, ctrl, &
                           margin_frac, tau_cap, q_land, q_marine, u_law, u_def, info, stat])
                                   ! implied mass balance by joint inversion on the geometry of ps: the speed U
                                   ! per ice cell (along the plastic flow direction) and smb = dH/dt + div(H U n)
                                   ! (upwind face speed; land outflow ablated, left out), minimising the misfit
                                   ! of ln U to the closure's speed (sigma_u, smoothness over ctrl%smooth) and of
                                   ! smb to smb_prior (sigma_smb_rel, _min), a penalty on smb < 0 where the prior
                                   ! is >= 0 off the margin (cold), U >= the deformation speed; Gauss-Newton
                                   ! (plastic_inv_ctrl_t, plastic_inv_info_t); m ice/yr, m/yr
l = plastic_has_lis()              ! built with lis=1
u = plastic_basal_velocity(slide, tau_b); u = plastic_deformational_velocity(rheo, tau, H)
```

Units: m, Pa, velocities in m/yr (per-year `A` and sliding coefficients), smb in m ice / yr.

Differences from Julia:

- Sea level is a field `z_sl(nx,ny)`: `H_flt = (rho_w/rho_i) max(z_sl - z_bed, 0)`.
- `tau_b` and `smb` are fields (no scalar broadcast).
- `plastic_smb_implied` takes the velocity vector; the flux magnitude is `hypot(ux,uy)`,
  or `speed` when given. Passing the `speed` of `plastic_velocity` reproduces Julia's
  `smb_from_velocity` (which reads `vel.speed`) and makes the round trip with
  `plastic_balance_flux` exact.
- "No upwind value" in the sweep is a finite sentinel (`1e300`) instead of `Inf`.
- Ice cells in a region with no ice-free cell (unreachable by the sweep) get `H = 0`
  and `stat = PLASTIC_NOT_CONVERGED` with a message (Julia returns `Inf`).

Status: nothing stops inside the library. `stat` (and `ps%stat`, `ps%msg`) is
`PLASTIC_OK` (converged), `PLASTIC_NOT_CONVERGED` (result usable, tolerance not met)
or `PLASTIC_BAD_INPUT` (wrong shape, unknown mode or law; nothing computed).

All work arrays are allocated in `plastic_init`; a solve allocates nothing, and there is
no module state, so several instances can coexist (OpenMP-safe).

### Warm start

`plastic_t` keeps `w = H^2` and the mask of the last solve. With `warm=.true.` in
`"full"` mode the outer fixed point starts from the stored `w` instead of the flat
solve:

- ice cells that were ice before keep the stored `w`; ice-free cells take the new
  margin value `H_flt^2`;
- cells new to the mask take the flat solution (one extra eikonal solve);
- if more than `par%warm_frac_max` (default 10 %) of the ice cells are new, the solve
  starts cold.

Each inner sweep still starts from the sentinel on the ice cells (the monotone scheme
needs it), so `warm` has no effect in `"flat"` mode. `plastic_solve_smb` uses the same
rule with the stored `H` (cold value: flat solve at `ctrl%tau_init`).

Warm and cold solves both stop on an outer step `max|dH| < outer_tol`; with the linear
rate rho ~ 0.7 of the damped map each is within ~rho/(1-rho) outer_tol of the fixed
point, so they agree to a few `outer_tol` (and converge to the same fixed point as
`outer_tol -> 0`, see the tests).

## Margin position

The margin lies between an ice cell and an ice-free neighbour. `margin_frac(nx,ny)`
(optional, on the ice-free cells next to the ice) is the distance from the ice-free cell
centre to the margin as a fraction of the spacing, in [0, 1): the neighbour enters the
Godunov update at (1 - margin_frac) dx (dy), with its margin value (z_b + H_flt, or
H_flt^2). The factor is floored at `MARGIN_DMIN` = 0.01, so a margin closer than 1 % of a
cell to an ice-cell centre counts as at 1 %. Values outside [0, 1) on ice-free cells give
`PLASTIC_BAD_INPUT`. Without `margin_frac`, `par%margin` sets one value: `"face"`
(default, 0.5: margin at the cell face) or `"centre"` (0: at the ice-free cell centre,
as PlasticIceSheet.jl). The rule applies to all modes and to `plastic_solve_smb`, and to
the residual diagnostics.

Signed-distance convention (kryos): `sdist(nx,ny)` [m] is the distance to the margin,
positive in the ice and negative outside, at the cell centres.
`plastic_margin_frac_from_sdist(sdist, mask, dx, dy, frac)` returns
frac = -sdist / h, h = sqrt(dx dy), clamped to [0, 1 - MARGIN_DMIN], on ice-free cells
with an ice neighbour, 0 elsewhere. The per-cell value is isotropic: exact for margins
along a grid axis; for an oblique margin the exact axis distance would need the pair of
sdist values across each face (not implemented).

With unequal spacings on the two sides of an axis, the upwind term is the Godunov
max_s (u - a_s)+/d_s over the sides, so the update is the minimum over the side
combinations (up to four kernel calls, only at margin cells). Every candidate is monotone
in the neighbour values, hence so is the minimum: the flat (w-eikonal) and surface
schemes stay monotone. With equal spacings the update is the previous one: `"centre"`
is bit-identical to PlasticIceSheet.jl (comparison rerun), and the margin handling costs
~5 % (300 x 300: flat 5.6/6.1 ms, surface 8.9/9.2 ms, centre/face).

Tests (`tests/test_surface.f90`):

- 1-D Nye strip, land margins at f = 0, 0.5, (0.3, 0.77), (0.95, 0.05) of the end cells,
  dx = 4 km to 250 m: surface and flat equal sqrt(2 k d) to < 2e-13 (relative) for every f,
  with d the distance to the margin: no half-cell offset.
- Circular cap, R = 30 km, centre off the grid; volume relative error of the solution
  (against the exact H on the same cells):

| dx | centre | face | sdist |
|---|---|---|---|
| 4 km | +10.1 % | -0.83 % | -0.54 % |
| 2 km | +4.2 % | -1.9 % | -0.27 % |
| 1 km | +2.3 % | -0.96 % | -0.12 % |
| 500 m | +1.1 % | -0.56 % | -0.07 % |

  first order for centre and sdist; face has the opposite sign (a staircase margin is on
  average not at the face). Surface == flat with margin_frac on the flat bed.

### Native-grid check against ICESHEET

`compare/native` (copies of the drivers in `/work/ba1442/robinson/plastic-diag`):
ICESHEET (PaleoMIST 1.0) 5 km North America output, PlasticIceSheet `"surface"` on the
same bed, tau_b and mask, inputs block-averaged to 10-30 km (ice where >= half the block).
Volume ratio PlasticIceSheet / ICESHEET over the largest ice body; sdist = signed distance
to the 5 km ICESHEET margin at the coarse cell centres.

| dx | Greenland 0 ka: centre | face | sdist | NA 20 ka: centre | face | sdist |
|---|---|---|---|---|---|---|
| 5 km | 1.017 | 1.009 | 1.009 | 1.007 | 1.005 | 1.005 |
| 10 km | 1.063 | 1.047 | 1.050 | 1.029 | 1.025 | 1.025 |
| 20 km | 1.089 | 1.055 | 1.072 | 1.060 | 1.051 | 1.053 |
| 30 km | 1.113 | 1.060 | 1.087 | 1.082 | 1.070 | 1.073 |

The margin position accounts for 0.4-5 points of the excess at 10-30 km (centre to sdist or face). The rest is not the metric
(the ratio against the 5 km thickness sampled at the coarse centres is within 0.002 of
the block-mean ratio) but the coarse inputs and the first-order scheme: solving at 5 km
with the bed, mask and tau_b block-averaged to 30 km already gives 1.039 (0 ka) and 1.029
(20 ka) (plastic-diag variants). "face" beats "sdist" at 20-30 km on Greenland because the
coarse mask (>= 50 % of the block) leaves partly iced cells outside it, and moving the
margin inward offsets part of the other errors; sdist is the more faithful margin.

## Mode "surface"

Not in PlasticIceSheet.jl. The ICESHEET formulation: the plastic surface equation

    |grad z_s| = tau_b / (rho_i g (z_s - z_b))

is solved for u = z_s itself by one fast-sweeping solve, with no outer iteration and
no warm start (`warm` has no effect). Ice-free cells hold z_s = z_b + H_flt(z_b, z_sl);
H = z_s - z_b on the ice. `par%tol` applies to z_s [m].

Why: the right-hand side depends on z_s only through the cell's own thickness, so the
Hamiltonian |p| - tau_b/(rho_i g (u - z_b)) is increasing in u (proper) and the Godunov
update is monotone: the comparison principle and fast sweeping apply directly. The
`"full"` mode solves the same continuous equation for w = H^2, |grad w + 2 sqrt(w) grad z_b|
= 2 tau_b/(rho_i g), whose Hamiltonian is increasing in w only where grad z_s . grad z_b > 0
(not on reverse bed slopes), and handles the bed term by a lagged, central-difference
fixed point; that map has no attracting fixed point on GRL-32KM or NH-32KM.

Discretisation (`godunov_surface` in `src/plastic_eikonal.f90`): per axis, the upwind
difference times the edge-mean thickness,

    P_a(H) = (H - alpha)+ (H + alpha+) / 2,   alpha = a - z_b,ij   (a: lower x-neighbour surface)
    (P_a/dx)^2 + (P_b/dy)^2 = (tau_b/(rho_i g))^2.

P_a increases with H and does not increase with a (dP_a/da = z_b - a <= 0 for a >= z_b, and
-H/2 below), so the scheme is proper and monotone; the update is the unique root, in closed
form (the w-eikonal with w_a = alpha^2) when both neighbour surfaces are above the bed,
else by Newton from the smaller one-sided root (monotone, P^2 convex). On a flat bed it
is identical to `"flat"` (exact 1-D Nye profile).

The literal cell-centred form, (u - a)+^2/dx^2 + ... = (tau_b/(rho_i g (u - z_b)))^2, is also
monotone but inaccurate at the margin: on a flat bed its first ice cell is
H = sqrt(tau_b dx/(rho_i g)) instead of the Nye value sqrt(2 tau_b dx/(rho_i g)), 29 % low at
every dx, with order ~0.75 convergence at a fixed distance (test_surface). Hence the
edge mean.

Tests (`tests/test_surface.f90`):

| test | result |
|---|---|
| kernel, 2e5 random cases (neighbours above/below the bed, dx /= dy) | solves its equation to 5e-11, H > 0, nondecreasing in a and b |
| 1-D Nye strip, dx = 4 km ... 250 m | surface and flat exact (< 1e-13); cell-centred variant 29 % low at the first cell, 2.3 % at 8 km for dx = 250 m |
| circular cap, flat bed, dx = 2 km ... 250 m | surface == flat (< 7e-11 m); mean error r < R/2: 6.7, 2.4, 1.4, 0.66 m (order ~1) |
| 1-D sloped bed (5e-3) vs ODE, dx = 4 km ... 500 m | max error surface 9.3, 4.7, 2.4, 1.2 m (order 0.98); full 6.1, 3.3, 1.8, 0.9 m; flat 313 m |
| bed bump (300 m) where full converges, dx = 2, 1, 0.5 km | max abs(surface - full) (r < 24 km) 11.3, 6.2, 3.4 m; max abs(flat - full) 117 m |

The residual of `plastic_stress_residual` is cell-centred (H of the cell): it is exactly
1 where the edge mean is (H + 0)/2, i.e. at the first ice cell and where a neighbour
surface lies below the cell's bed, and O(dx) elsewhere. Its scheme-consistent
counterpart is `plastic_residual` (the Godunov residual of z_s), < 1e-8 m.

## SMB-driven reconstruction, method "routing": why it does not converge

`plastic_solve_smb` (Julia `diva_reconstruct`) iterates H -> T(H): route the SMB down
z_s = z_b + H with the MFD weights (flux Q), invert the DIVA flux law per cell for tau_d,
solve the plastic equation with tau_d (inner mode `ctrl%mode`, `"surface"` by default,
`"flat"` as Julia), and under-relax, H <- H + r (T(H) - H). The `"surface"` inner solve
does not make it converge. `compare/smb_diag.f90`, 2000 iterations at r = 0.05 from the
cold start:

| case (forcing) | inner | final max abs(T(H) - H) [m] | rate (2nd half) | time per iteration |
|---|---|---|---|---|
| synth150 (closure-B SMB of the flat solution) | surface | 351 | 0.9996 | 21 ms |
| synth150 | flat | 624 | 0.9999 | 19 ms |
| grl32 (MAR v3.14, Coulomb) | surface | 305 | 0.999 | 4.2 ms |
| grl32 | flat | 757 | 0.9996 | 3.9 ms |
| nh32_0ka (implied SMB; no SMB product on NH-32KM) | surface | 1719 | 1.0001 | 8.2 ms |

With the implied SMB of the plastic solution (`forcing implied`) H0 is an exact fixed
point (T(H0) - H0 < 4e-11 m), yet the iteration from the cold start ends 370-900 m rms
away from it. Its linear stability, from power iteration on the finite-difference
Jacobian J of T at H0:

| case | inner | dominant eigenvalue of J | eigenvector | stable for r < 2/(1 - lambda) |
|---|---|---|---|---|
| grl32 | surface | -188 | interior | 0.011 |
| grl32 | flat | -51 | interior | 0.038 |
| synth150 | surface | -596 | interior | 0.0034 |
| synth150 | flat | -234 | interior | 0.0085 |

The bound is exact: on grl32/flat a 1 m perturbation of H0 decays at r = 0.02 (6.5e-5 m
after 400 iterations) and grows to ~100 m at r = 0.05.

Cause: the MFD weights split a cell's flux by the surface drops to its neighbours, which
are slope x dx (a few metres). A change dH of the surface moves the split by O(dH / (slope
dx)), Q by the same fraction, tau_d ~ Q^(1/3)..., and the swept thickness downstream in the
opposite direction: |lambda| ~ H / (slope dx), large and growing with resolution. The
map is also discontinuous: routing changes (sort order, drops through zero) and cells
with Q <= 0 (tau_d = 0, the surface goes flat) carry the largest moves (`smb_diag`
classes). Under-relaxation would have to scale like slope dx / H, case by case, and still
fails from a far start. This is structural, not a tuning problem.

Structural fix (implemented as `ctrl%method = "steady"`, below; first prototyped with
explicit pseudo-time steps in `compare/smb_sia_proto.f90`):
solve mass conservation directly, div q = smb on the mask with the same membrane-free
DIVA flux, q = -H (u_b(tau) + u_def(tau, H)) grad z_s/|grad z_s|, tau = rho_i g H |grad z_s|.
This is the steady SIA with sliding, a monotone degenerate-elliptic problem (obstacle
H >= 0) with a unique solution; the flux depends on the local slope continuously, not on
a routing partition. The prototype uses explicit pseudo-time steps (staggered Mahaffy
fluxes, dt = 0.2 dx^2 / max D, ice-free cells fixed at z_b + H_flt):

| case (forcing) | result |
|---|---|
| grl32 (MAR) | steady: max abs(dH/dt) 1e-9 m/yr after 20552 steps (46 kyr), 6.4 s; 21 mask cells end ice-free |
| grl32, nh32_0ka (implied) | interior steady; margin cells cycle on the H = 0 obstacle (max abs(dH/dt) 24, 229 m/yr): the MFD-implied SMB puts all terminal flux into the margin cells as ablation, while the SIA flux leaves through the margin faces, so this forcing is inconsistent with it |
| synth150 (dx = 2 km, implied) | dt = 2e-3 yr: 1e6 steps (~30 min) reach only ~2 kyr, rms 35 m from the plastic solution (explicit stepping too slow at km resolution); same margin cycling |

The explicit prototype was too slow at km resolution, and the MFD-implied SMB did not match
its margin treatment; the steady method below solves both (implicit Newton, and an
implied SMB built on the same face fluxes).

## SMB-driven reconstruction, method "steady": direct formulation

`ctrl%method = "steady"` (default) solves steady mass conservation directly; `"routing"`
keeps the Julia iteration (conditionally stable, see above).

### Equations

Unknown: H >= 0 on the ice cells (mask); surface s = z_b + H. Membrane-free DIVA flux
(the closures' flux law with the driving stress taken from the geometry):

    q = -H U(tau, H) grad s / |grad s|,  tau = rho_i g H |grad s|,  U = u_b(tau) + 2A/(n+2) tau^n H

(Coulomb: tau capped at 0.999 C). Steady state with the obstacle H >= 0 is the
complementarity problem

    H >= 0,  F(H) = div q - smb >= 0,  H F(H) = 0

(F > 0 at H = 0: ablation exceeds the inflow, the cell stays ice-free).

### Discretisation (finite volumes, staggered fluxes)

Interior face between ice cells P and E (x shown, y alike): H_f = (1 - w) H_up + w H_dn
(up/down by surface), w = min(1/2, H_up/(4 |s_P - s_E|)): the mean (Mahaffy) where the
surface drop across the face is small against the thickness, limited towards the upstream
value where it is not. The limit keeps the face flux non-increasing in the downstream
thickness (dq/dH_dn <= 0 for q ~ H^(p+1) g^p, p >= 1: Glen and all sliding laws), so F_P is
non-increasing in its neighbours' H (an M-function). With the plain mean it is not, and
the active-set Newton stalled on cells next to H = 0 (GRL-32KM at 2e-2 m/yr); with full
upwinding it converges but the Vialov error grows from 0.3 % to 2.3 % at 2.5 km.
s_x = (s_E - s_P)/dx, s_y = mean of the central y-differences of s at P and E,
|grad s| = sqrt(s_x^2 + s_y^2), tau_f = rho_i g H_f |grad s|,
q_x = -H_f U(tau_f, H_f) s_x/|grad s| (0 where grad s = 0). Ice-free cells carry
s = z_b + H_flt (the plastic margin value). Faces on the array edge carry no flux.

Margin face (ice cell P, ice-free neighbour M). The margin lies d_c = (1 - margin_frac_M) dx
from the centre of P (the margin position of the plastic solve). Derivation: on the
segment from P to the margin the flux is the membrane-free flux with the one-sided
slope over that segment and the segment-mean thickness,

    g = (s_P - s_M)/d_c,  H_m = (H_P + H_flt,M)/2,  |grad s| = sqrt(g^2 + s_y,P^2),
    tau_m = rho_i g H_m |grad s|,  q_out = max(H_m U(tau_m, H_m) g/|grad s|, 0)

(s_y,P the central y-difference at P), leaving P through the face towards M; nothing
enters from an ice-free cell. Why not the plastic profile H ~ sqrt(2 k d) at bounded tau:
its flux tends to zero at a land margin, so with the margin at or inside the face no ice
could leave and a positive SMB up to the margin (the Vialov case) would have no steady
state. The one-sided SIA flux grows without bound with H_P, so every SMB has a steady
state; it is the Mahaffy face flux with the margin as the neighbour, and with
margin_frac it uses the same sub-cell margin position as the eikonal.

Cell balance: F_P = (q_e - q_w)/dx + (q_n - q_s)/dy - smb_P.

Implied SMB (closure B for this method): smb = (q_e - q_w)/dx + (q_n - q_s)/dy of the
same face fluxes on a given geometry (`plastic_smb_flux_divergence`). The plastic
solution forced by its own implied SMB is then a root of F, and the steady solve returns
it (fixed point; exact up to the solver tolerance).

### Solver

Pseudo-transient continuation with a semismooth Newton step per pseudo-time step:

    G(H) = (H - H_k)/dt_k + F(H),  Phi(H) = min(H, c_P G(H)) = 0,  c_P = 1/(1/dt_k + dF_P/dH_P).

- Active set A = {cells with H <= c_P G}: the cell's own Newton step would empty it (c_P = dt
  where the diagonal is not positive). A fixed c (1 yr) zeroed thick, slowly thinning
  cells at large dt; c = dt did the same. Equation H = 0 on A. Inactive cells: the Newton row
  of G, (1/dt) d_P + sum_K dF_P/dH_K d_K = -G_P, with d_K = -H_K for K in A moved to the
  right-hand side. New iterate max(H + lam d, 0), lam = 1, 1/2, ... (up to 6 halvings)
  until the rms residual decreases.
- Jacobian dF/dH by finite differences with 9 colours: F_P depends on H in the 3 x 3
  block around P, so perturbing every third cell in x and y at once gives one column
  per cell per residual evaluation (10 residual evaluations per Newton step; works for all
  three sliding laws).
- dt: switched evolution relaxation on the rms residual, dt_{k+1} = dt_k ||F_k||/||F_{k+1}|| (at least 1.5 dt_k
  after a step that did not increase the residual, so that a slowly decaying residual
  cannot stall dt), within
  [ctrl%dt0, ctrl%dt_max]; a step that increases the rms residual ten-fold is rejected and
  dt halved. The pseudo-time term keeps the system well posed at divides (where dF/dH
  vanishes with |grad s|) and at the start.
- Converged when max over the cells of |min(H, F)| (F alone on ice, max(-F, 0) on
  ice-free cells) < ctrl%res_tol [m/yr].
- Start: the stored H (warm) or a surface-mode plastic solve at ctrl%tau_init (cold).

Linear solve (9-point, non-symmetric): CSR matrix over the inactive ice cells.
- `lis=1` (ifx; fesm-utils' Lis): BiCGStab with ILU(0) through `lis_solver_set_option`
  ("-i bicgstab -p ilu -tol ... -maxiter ..."), as yelmo's `solver_linear`.
- Built-in fallback: BiCGStab with ILU(0), no dependency.

### OpenMP and determinism

Residual, Jacobian colouring, CSR assembly and the mat-vec run in parallel loops over
rows (each row written by one thread: deterministic). Dot products sum fixed blocks of
4096 entries in parallel and the block sums sequentially, so they do not depend on the
thread count; the ILU(0) triangular solves are sequential. The built-in solver is
therefore bit-reproducible for any thread count (tested: 1 vs 4 threads). Lis uses OpenMP
reductions: with `lis=1` results differ between thread counts at the level of the linear
tolerance (observed 6e-10 to 1.5e-9 m on GRL-32KM). The fast sweeping has an optional
wavefront ordering (`par%sweep = "wavefront"`, cells of one anti-diagonal in parallel),
the same Gauss-Seidel order and bit-identical to the serial sweep (tested in all modes).

### Results

Tests (`tests/test_steady.f90`, built-in solver unless noted):

| test | result |
|---|---|
| plastic fixed point under its implied SMB, bump 40 x 40 at 2 km | H_p is a root (0 steps from H_p); cold start 16 steps, 0.018 s, max abs(H - H_p) 2e-12 m |
| same, cap 150 x 150 at 2 km | cold start 37 steps, 1.2 s (1 thread), max abs(H - H_p) 2e-11 m |
| Lis (lis=1, 8 threads) vs built-in | 5e-13 m |
| built-in, 1 vs 4 threads | bit-identical |
| Vialov profile (n = 3, a = 0.3 m/yr, L = 300 km), dx = 20, 10, 5, 2.5 km | max rel error (x <= 0.9 L) 2.6, 1.1, 0.61, 0.32 % (order 0.9-1.2); volume +0.87, +0.45, +0.23, +0.12 %; 22-23 steps |

Real data and thread scaling (`compare/bench_threads.sh`, Levante shared node, ifx,
openmp=1 lis=1; wall time of plastic_solve_smb from the cold start):

| case | steps | residual [m/yr] | 1 thread | 4 | 8 | 16 (built-in / Lis) |
|---|---|---|---|---|---|---|
| grl32, MAR v3.14, Coulomb | 28 | 4e-8 | 0.09 / 0.10 s | 0.055 / 0.071 | 0.049 / 0.068 | 0.054 / 0.077 |
| synth150 (its input SMB) | 31 | 2e-10 | 1.37 / 1.27 s | 0.92 / 0.66 | 0.85 / 0.53 | 0.93 / 0.55 |
| nh32_0ka, implied SMB of its plastic solution | 32 | 2e-10 | 0.50 / 0.52 s | 0.30 / 0.34 | 0.27 / 0.28 | 0.28 / 0.29 |

On nh32_0ka the solver returns the plastic solution (max abs(H - H_plastic) 6e-10 m). On
grl32 with MAR, 17 of the 1710 mask cells end ice-free (ablation exceeds the inflow;
mean H 1812 m). The scaling is modest (the residual evaluations of the coloured Jacobian
dominate, the built-in ILU is sequential); Lis helps on the larger grid.

Eikonal, surface mode, best of 5 (serial / wavefront), same job:

| grid | 1 thread | 4 | 8 | 16 |
|---|---|---|---|---|
| synth300 (300 x 300) | 12.4 / 14.1 ms | 12.4 / 12.2 | 12.4 / 12.4 | 12.4 / 16.4 |
| nh32_20ka (276 x 276) | 14.0 / 16.2 ms | 14.0 / 13.8 | 14.0 / 14.6 | 18.7 / 16.4 |

An earlier run of the same job gave 8.6 / 6.9 / 8.5 ms (wavefront, 4 / 8 / 16 threads on
synth300): the gain is at most 1.8x and not reproducible on a shared node. A solve takes
~10 ms and ~10 sweeps; each sweep has nx + ny barriers. Not worth it at these sizes:
`par%sweep = "serial"` stays the default.

Status of `ctrl%method = "routing"`: kept for parity with PlasticIceSheet.jl; stable
only for relax < 2/(1 - lambda), lambda ~ -H/(slope dx) (0.003-0.04 in the cases above),
see the diagnosis above. Use "steady".

## Build

```bash
configme config PlasticIceSheet            # once per machine/compiler: ./Makefile from config/Makefile
make static                                # -> libplastic/include/libplastic.a + .mod
make check                                 # build and run the tests
make clean                                 # between compilers (.mod files differ)
```

The Makefile is generated by [configme](https://github.com/fesmc/configme) from the template
`config/Makefile`: configme's compiler and machine fragments (FC, FFLAGS, DFLAGS_*, the machine's
-march) in place of its placeholder, then `config/common.mk` (the dependency wiring). The repository
keeps no machine files. `config/common.mk` appends the compiler's value-safe floating-point flag
(`-fp-model precise` for ifx, `-fno-fast-math` for gfortran) to the machine's -Ofast: the library
tests NaN.

Output: `libplastic/include/libplastic.a` and the `.mod` files in `libplastic/include`.
Link a host with

```make
PLASTICROOT  = PlasticIceSheet
INC_PLASTIC  = -I$(PLASTICROOT)/libplastic/include
LIB_PLASTIC  = -L$(PLASTICROOT)/libplastic/include -lplastic
```

Options: `debug=1`, `openmp=1` (OpenMP in the steady solver, the built-in linear
solver and the optional wavefront sweep), `pic=1`, `lis=1` (Lis for the steady solver's
linear solves; fesm-utils' Lis: `LISROOT` defaults to `fesm-utils/lis/lis-omp` with
openmp=1, `lis-serial` otherwise, `fesm-utils` the checkout or link configme makes here). With lis=1 a host links `-L$(LISROOT)/lib -llis` after `-lplastic`:

```make
make -C PlasticIceSheet static openmp=1 lis=1 LISROOT=$(LISROOT)
LIB_PLASTIC = -L$(PLASTICROOT)/libplastic/include -lplastic -L$(LISROOT)/lib -llis
```

## Tests

`make check` runs `tests/test_core.f90`, `tests/test_closures.f90`,
`tests/test_surface.f90` (see [Mode "surface"](#mode-surface)) and
`tests/test_steady.f90` (the steady SMB-driven solver): the Julia test
set (godunov kernel, Nye profile, flat vs full, full over bed relief and its residual,
constitutive laws, DIVA velocity, closure B round trip `balance_flux(smb_implied) = u H`
to 1e-16, balance flux and A <-> B consistency, implied basal velocity, the SMB-driven
fixed point and cold-start recovery), plus a marine margin with a 2-D sea level
(analytic `H = sqrt(H_flt^2 + c d)`), warm == cold start (full and SMB-driven), and the
status reporting. All pass with ifx and gfortran, release and debug builds.

## Comparison with PlasticIceSheet.jl

`compare/run.sh` (sbatch, Levante shared partition, one core) writes the inputs of each
case from Julia (`compare/julia_ref.jl`), runs PlasticIceSheet.jl (1.12.6, after a JIT
warm-up) and the Fortran driver (`compare/plastic_compare.f90`, ifx 2023.2, `-O3
-march=core-avx2 -fp-model precise`) on them, and tabulates the differences over the ice
mask (`compare/compare.jl`). Times are wall clock, best of 3 for flat/full, one run for
the SMB-driven solve. All Fortran runs of the comparison use `par%margin = "centre"`
(Julia's margin position), including the surface-mode results below.

- `synth150`, `synth300`: synthetic cap (dx = 2 km, circular mask, 300 m bedrock high plus
  50 m undulations, tau_b = 1e5 Pa +- 30 %); full with relax = 0.4, n_outer = 200;
  SMB-driven (Weertman) from the closure-B SMB of the flat solution, cold start.
- `grl32`: Greenland GRL-32KM (FesmData: BedMachine v5 bed and grounded mask, MAR v3.14
  1961-1990 SMB), uniform tau_b = 1e5 Pa; full with relax = 0.4, n_outer = 500;
  SMB-driven with a regularized-Coulomb bed, C = rho_i g mean(H_obs), as
  `examples/reconstruct_greenland.jl`. Sea level 0. The Fortran SMB-driven runs use
  `ctrl%mode = "flat"` (Julia's inner solve).

| case | field | max H [m] | max abs dH [m] | rms dH [m] | t Julia [s] | t Fortran [s] | Julia / Fortran |
|---|---|---|---|---|---|---|---|
| synth150 | H_flat | 1603.6 | 0 | 0 | 0.0011 | 0.0010 | 1.1 |
| synth150 | H_full | 1467.7 | 0 | 0 | 0.1084 | 0.0840 | 1.3 |
| synth150 | H_smb | 1261.7 | 0 | 0 | 19.04 | 19.08 | 1.0 |
| synth300 | H_flat | 2261.9 | 0 | 0 | 0.0063 | 0.0054 | 1.2 |
| synth300 | H_full | 2135.6 | 0 | 0 | 0.8002 | 0.4690 | 1.7 |
| grl32 | H_flat | 3180.1 | 0 | 0 | 0.0004 | 0.0002 | 1.6 |
| grl32 | H_full | 3332.6 | 0 | 0 | 0.3321 | 0.2252 | 1.5 |
| grl32 | H_smb | 2418.1 | 0 | 0 | 4.546 | 3.858 | 1.2 |

The Fortran results are bit-identical to Julia in every cell, including the iterations
that did not converge. The speed is comparable (Julia's loops are already compiled
code); the gain of the port is the dependency-free coupling, not raw speed.

Convergence (identical in both codes):

| case | full: n_outer (converged) | SMB-driven: n_outer (converged) |
|---|---|---|
| synth150 | 44 (yes) | 1000 (no, cap) |
| synth300 | 62 (yes) | not run |
| grl32 | 500 (no) | 1000 (no) |

On GRL-32KM the `"full"` fixed point does not converge for any relax tested
(0.4 to 0.02, 2000 iterations): the last step stays at 3-200 m and ~2/3 of the ice cells
change by > 100 m between consecutive iterates, a global oscillation. Use `"flat"` at
32 km (as `examples/invert_taub_greenland.jl` does). The SMB-driven reconstruction on
GRL-32KM does not converge either (5000 iterations, last step ~950 m).

### Warm start

Fortran, coupling-step-like change, warm vs cold `"full"` solve (outer_tol = 1e-3):

| case | change | n_outer warm / cold | time warm / cold [s] | speedup | max abs (H_warm - H_cold) [m] |
|---|---|---|---|---|---|
| synth150 | 1 m GIA-like bed change | 27 / 44 | 0.050 / 0.084 | 1.7 | 1.3e-3 |
| synth300 | 1 m GIA-like bed change | 37 / 62 | 0.266 / 0.478 | 1.8 | 2.6e-3 |
| test (80x80) | 1 m bed change, outer_tol 1e-4 | 27 / 41 | | 1.5 | 1.8e-4 |
| test (80x80) | unchanged inputs | 1 / 41 | | | 0 |

The gain is limited by the linear rate (rho ~ 0.7) of the damped fixed point: from an
initial error e0 it needs ~log(e0/outer_tol)/log(1/rho) iterations, so even a small
change costs tens of iterations. Where the fixed point does not converge (grl32) the warm
start gains nothing.

## Mode "surface" on real data

Same driver and job (`compare/plastic_compare.f90`), sea level 0, uniform tau_b = 1e5 Pa:
GRL-32KM as above; NH-32KM (276 x 276) with the FesmData GEBCO2026 bed and the mask
PaleoMIST H_ice > 0 at 20 ka (`nh32_20ka`, 27669 ice cells) and 0 ka (`nh32_0ka`, 3055).

| case | surface: sweeps, t [s] | flat t [s] | full t [s] (converged) | Godunov residual [m] |
|---|---|---|---|---|
| synth150 | 0.0027 | 0.0010 | 0.084 (yes, 44) | |
| synth300 | 0.0105 | 0.0055 | 0.48 (yes, 62) | |
| grl32 | 9, 0.0007 | 0.0002 | 0.23 (no, 500) | 0 |
| nh32_20ka | 12, 0.0116 | 0.0044 | 9.5 (no, 500) | 7.5e-9 |
| nh32_0ka | 9, 0.0028 | 0.0012 | 1.2 (no, 500) | 0 |

Surface converges in one sweep solve everywhere, at 2-3 times the cost of flat and
40-800 times less than full. Where full converges (synthetic caps) surface agrees with
it to rms 2.7 m (dx = 2 km, 150 x 150) and 1.4 m (300 x 300, same dx).

Bed relief, surface - flat over the ice:

| case | mean [m] | rms [m] | max [m] | where (max) | mean abs by bed slope < 0.005 / 0.005-0.02 / > 0.02 [m] | cells > 100 m |
|---|---|---|---|---|---|---|
| synth150 | | 71.6 | -178 | | | |
| grl32 | +95 | 291 | +1106 | x = 240, y = -1690 km (north-east interior; bed -75 m, H_flat 2668 m) | 220 / 223 / 304 | 1183 / 1710 |
| nh32_0ka | -73 | 278 | +1054 | x = 252, y = -1688 km (same place; bed -29 m) | 195 / 222 / 368 | 2059 / 3055 |
| nh32_20ka | -157 | 413 | +2662 | x = -3620, y = 1960 km (bed -4414 m, slope 0.054) | 224 / 454 / 687 | 18174 / 27669 |

Bed relief changes H by a few hundred metres over most of the ice at 32 km, and more on
steep beds. Flat drapes the thickness over the bed (z_s follows z_b); surface keeps the
surface continuous, so it thickens the ice over inland bed depressions (the interior
behind the coastal mountains) and thins it over highs. The 20 ka maximum is a
deep-ocean cell inside the PaleoMIST mask (H_ice > 0 includes ice over 4 km of water),
where flat puts the surface 2 km below sea level.

Stress residual | |grad z_s| rho_i g H - tau_b | / tau_b (cell-centred, upwind), median
over the ice: surface 0.04 (grl32), 0.02 (nh32_20ka), 0.14 (nh32_0ka), max 1.0 (see
above); flat 0.49 / 0.39 / 0.76, max up to 22; full (not converged) 0.30 / 0.21 / 0.34,
max up to 16.

Recommendation for kryos: `mode = "surface"` (one sweep solve, no fixed point, bed
relief included, converges on GRL-32KM and NH-32KM). `"flat"` remains the cheapest and
equals `"surface"` on flat beds; `"full"` only reproduces PlasticIceSheet.jl.
