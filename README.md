# PlasticIceSheet

Perfectly-plastic, steady-state ice-sheet reconstruction in Fortran: a dependency-free
port of [PlasticIceSheet.jl](https://github.com/fesmc/PlasticIceSheet.jl). The plastic
surface equation `|grad z_s| = tau_b / (rho_i g H)` is solved for `w = H^2` as a static
eikonal problem by fast sweeping, on a prescribed grounded-ice mask (`"flat"` mode, or
`"full"` with the bed-gradient correction by a damped fixed point and an optional warm
start). Mode `"surface"` (not in the Julia package, recommended) solves the surface
equation for z_s directly, with the bed relief and no outer iteration. Closures: Glen
deformation, linear / Weertman / regularized-Coulomb sliding, DIVA velocity, MFD balance
flux and velocity, implied SMB, and the SMB-driven reconstruction (steady mass
conservation solved directly by Newton; optional Lis and OpenMP).

```bash
configme config PlasticIceSheet       # the Makefile for this machine (fesmc/configme)
make static                           # -> libplastic/include/libplastic.a + .mod
make check                            # tests
```

```fortran
use plastic_ice
type(plastic_t) :: ps
call plastic_init(ps, nx, ny, dx, dy)
call plastic_solve(ps, z_bed, z_sl, tau_b, mask)   ! ps%par%mode = "surface" -> ps%H, ps%z_srf
```

API, build, tests and the comparison with the Julia package: [docs/fortran-port.md](docs/fortran-port.md).

Clean-room port of the MIT-licensed Julia package; no code from ICESHEET (GPL).
