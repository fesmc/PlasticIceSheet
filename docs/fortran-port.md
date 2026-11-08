# Fortran port of PlasticIceSheet.jl

A dependency-free Fortran library (`libplastic.a`) with the forward solver and the
closures of [PlasticIceSheet.jl](https://github.com/fesmc/PlasticIceSheet.jl), for
coupling into kryos as an ice backend. Not ported: the AD / implicit-differentiation
extension and the tau_b inversion.

## API

`use plastic_ice`. All reals are `real64` (exported as `plastic_dp`); arrays are
`(nx,ny)`, `i = x`, `j = y`, as in Julia.

```fortran
type(plastic_par_t)  :: par   ! rho_i=917, rho_w=1025, g=9.80665, H_min=1, mode="full"|"flat",
                              ! max_sweeps=200, n_outer=100, tol=1e-6, outer_tol=1e-3, relax=0.5,
                              ! warm_frac_max=0.1
type(plastic_t)      :: ps    ! par, nx, ny, dx, dy, w=H^2, H, z_srf, n_outer_used, converged,
                              ! stat, msg, n_sweeps, outer_change, warm_used, n_new_cells

call plastic_init(ps, nx, ny, dx, dy [, par, stat])
call plastic_solve(ps, z_bed, z_sl, tau_b, mask [, warm, stat])
call plastic_residual(ps, z_bed, z_sl, tau_b, mask, r [, stat])  ! Godunov residual of ps%w
call plastic_reset(ps)                                            ! drop the warm start
call plastic_end(ps)
H_flt = plastic_flotation_thickness(z_bed, z_sl, par)             ! elemental
V = plastic_volume(ps); A = plastic_area(ps, mask); m = plastic_surface_misfit(ps, z_srf_obs, mask)

! closures
type(plastic_glen_t)     :: rheo   ! A=1e-16 [Pa^-n yr^-1], n=3
type(plastic_sliding_t)  :: slide  ! law="linear"|"weertman"|"coulomb", beta, m, C, u_0
type(plastic_smb_ctrl_t) :: ctrl   ! n_outer=1000, relax=0.05, tol=1.0, tau_hi=5e5,
                                   ! max_sweeps=400, sweep_tol=1e-9, tau_init=5e4
call plastic_velocity(ps, tau_b, mask, rheo, slide, ux, uy [, speed, slope_eps, u_b, u_def, stat])
call plastic_smb_implied(ps, ux, uy, mask, smb [, speed, stat])
call plastic_balance_flux(z_srf, mask, smb, dx, dy, Q [, stat])
call plastic_balance_velocity(z_srf, H, mask, smb, dx, dy, speed, ux, uy [, Q, slope_eps, H_min, stat])
call plastic_implied_basal_velocity(z_srf, H, tau_b, mask, smb, dx, dy, rheo, u_basal &
                                    [, u_balance, u_def, Q, slope_eps, H_min, stat])
call plastic_solve_smb(ps, z_bed, z_sl, smb, mask, rheo, slide [, warm, ctrl, H_init, tau_d, stat])
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

## Build

```bash
make static config=dkrz_levante_ifx        # or dkrz_levante_gfortran, macbook_gfortran
make check  config=dkrz_levante_ifx        # build and run the tests
make clean                                 # between compilers (.mod files differ)
```

Output: `libplastic/include/libplastic.a` and the `.mod` files in `libplastic/include`.
Link a host with

```make
PLASTICROOT  = PlasticIceSheet
INC_PLASTIC  = -I$(PLASTICROOT)/libplastic/include
LIB_PLASTIC  = -L$(PLASTICROOT)/libplastic/include -lplastic
```

Options: `debug=1`, `openmp=1` (adds the OpenMP flag; the library has no OpenMP
regions), `pic=1`. The fragments `config/<machine>_<compiler>` set `FC`, `FFLAGS`
and `DFLAGS_*`; ifx uses `-O3 -march=core-avx2 -fp-model precise`.

## Tests

`make check` runs `tests/test_core.f90` and `tests/test_closures.f90`: the Julia test
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
the SMB-driven solve.

- `synth150`, `synth300`: synthetic cap (dx = 2 km, circular mask, 300 m bedrock high plus
  50 m undulations, tau_b = 1e5 Pa +- 30 %); full with relax = 0.4, n_outer = 200;
  SMB-driven (Weertman) from the closure-B SMB of the flat solution, cold start.
- `grl32`: Greenland GRL-32KM (FesmData: BedMachine v5 bed and grounded mask, MAR v3.14
  1961-1990 SMB), uniform tau_b = 1e5 Pa; full with relax = 0.4, n_outer = 500;
  SMB-driven with a regularized-Coulomb bed, C = rho_i g mean(H_obs), as
  `examples/reconstruct_greenland.jl`. Sea level 0.

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
