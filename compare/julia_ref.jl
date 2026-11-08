# Reference run of PlasticIceSheet.jl for the Fortran port comparison.
#
# 1. Writes the inputs of each case to out/in_<case>.nc (the Fortran driver reads the
#    same files): z_bed, z_sl, tau_b, mask, smb, and the solver/closure settings as
#    global attributes.
# 2. Runs solve(:flat), solve(:full) and diva_reconstruct, writes out/julia_<case>.nc
#    with H_flat, H_full, H_smb and the wall-clock times (after a JIT warm-up).
#
# Run (Levante): julia --project=<env with PlasticIceSheet + NCDatasets> compare/julia_ref.jl [outdir]

using PlasticIceSheet
using NCDatasets
using Printf
using Statistics

const OUT = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "out")
const GRL_DIR = get(ENV, "GRL_DIR", "/work/ba1442/fesmdata/Greenland/GRL-32KM")
mkpath(OUT)

const RHO_I, G = 917.0, 9.80665

# --- cases ------------------------------------------------------------------------

# Synthetic ice cap on an n x n grid: a circular mask, a smooth bedrock high plus
# undulations, a spatially varying tau_b, and the closure-B SMB of its flat solution.
function synthetic(n; dx = 2000.0)
    z_b = [300.0 * exp(-((i - 0.44n)^2 + (j - 0.56n)^2) / (2 * (0.15n)^2)) +
           50.0 * sin(i / (0.0875n)) * cos(j / (0.1125n)) - 100.0 for i in 1:n, j in 1:n]
    tau = [1.0e5 * (1 + 0.3 * sin(i / (0.1375n)) * cos(j / (0.1625n))) for i in 1:n, j in 1:n]
    c = (n + 1) / 2
    mask = [hypot(i - c, j - c) <= 0.375n for i in 1:n, j in 1:n]
    rheo = GlenRheology(A = 1.0e-16, n = 3.0)
    slide = WeertmanSliding(β = 1.0e4, m = 3.0)
    z_s, H = solve(z_b, tau, mask; dx, dy = dx, mode = :flat)
    vel = diva_velocity(z_s, H, tau, mask, dx, dx; rheology = rheo, sliding = slide)
    smb = smb_from_velocity(vel, H, z_s, mask, dx, dx)
    return (; z_b, z_sl = zeros(n, n), tau, mask, smb, dx, dy = dx,
            x = collect(0:n-1) .* dx, y = collect(0:n-1) .* dx,
            rheo, slide, law = "weertman", relax_full = 0.4, n_outer_full = 200,
            do_smb = n <= 200)
end

# Greenland GRL-32KM (FesmData): BedMachine v5 bed and grounded-ice mask, uniform
# tau_b = 1e5 Pa, MAR v3.14 1961-1990 SMB, regularized-Coulomb bed with the cap at
# the mean overburden (as examples/reconstruct_greenland.jl).
function greenland()
    topo = joinpath(GRL_DIR, "GRL-32KM_TOPO-BedMachine-v5.nc")
    clim = joinpath(GRL_DIR, "GRL-32KM_SMB-MARv3.14-1961-1990-annual-mean.nc")
    rd(ds, v) = Float64.(coalesce.(ds[v][:, :], NaN))
    xc, yc, z_b, H_obs, m = NCDataset(topo) do ds
        (Array{Float64}(ds["xc"][:]), Array{Float64}(ds["yc"][:]),
         rd(ds, "z_bed"), rd(ds, "H_ice"), Int.(ds["mask"][:, :]))
    end
    smb = NCDataset(clim) do ds
        rd(ds, "smb")
    end ./ RHO_I                                     # kg m-2 yr-1 -> m ice / yr
    smb[isnan.(smb)] .= 0.0
    z_b[isnan.(z_b)] .= 0.0
    H_obs[isnan.(H_obs)] .= 0.0
    mask = m .== 2                                   # grounded ice
    nx, ny = size(z_b)
    dx = abs(xc[2] - xc[1]) * 1000.0
    dy = abs(yc[2] - yc[1]) * 1000.0
    H_ref = mean(H_obs[mask])
    slide = RegularizedCoulomb(C = RHO_I * G * H_ref, u_0 = 100.0, m = 3.0)
    return (; z_b, z_sl = zeros(nx, ny), tau = fill(1.0e5, nx, ny), mask, smb, dx, dy,
            x = xc .* 1000.0, y = yc .* 1000.0,
            rheo = GlenRheology(A = 1.0e-16, n = 3.0), slide, law = "coulomb",
            relax_full = 0.4, n_outer_full = 500, do_smb = true)
end

slide_params(s::WeertmanSliding) = (s.β, s.m, 0.0, 0.0)
slide_params(s::RegularizedCoulomb) = (0.0, s.m, s.C, s.u_0)

function write_inputs(path, c)
    isfile(path) && rm(path)
    NCDataset(path, "c") do ds
        defDim(ds, "x", length(c.x)); defDim(ds, "y", length(c.y))
        defVar(ds, "x", c.x, ("x",)); defVar(ds, "y", c.y, ("y",))
        defVar(ds, "z_bed", c.z_b, ("x", "y"))
        defVar(ds, "z_sl", c.z_sl, ("x", "y"))
        defVar(ds, "tau_b", c.tau, ("x", "y"))
        defVar(ds, "mask", Int32.(c.mask), ("x", "y"))
        defVar(ds, "smb", c.smb, ("x", "y"))
        beta, m, C, u0 = slide_params(c.slide)
        for (k, v) in ("dx" => c.dx, "dy" => c.dy, "glen_A" => c.rheo.A, "glen_n" => c.rheo.n,
                       "slide_beta" => beta, "slide_m" => m, "slide_C" => C, "slide_u0" => u0,
                       "relax_full" => c.relax_full)
            ds.attrib[k] = Float64(v)
        end
        ds.attrib["slide_law"] = c.law
        ds.attrib["n_outer_full"] = Int32(c.n_outer_full)
        ds.attrib["do_smb"] = Int32(c.do_smb)
    end
end

# Best of nrep wall-clock times, and the last result.
function timed(f, nrep)
    best = Inf; r = nothing
    for _ in 1:nrep
        t0 = time_ns(); r = f(); best = min(best, (time_ns() - t0) / 1e9)
    end
    return r, best
end

function run_case(name, c)
    @printf("%s: %d x %d, dx = %.0f m, ice cells %d\n", name, size(c.z_b)..., c.dx, count(c.mask))
    write_inputs(joinpath(OUT, "in_$name.nc"), c)
    kw = (; dx = c.dx, dy = c.dy)
    (_, Hf), tf = timed(() -> solve(c.z_b, c.tau, c.mask; kw..., mode = :flat), 3)
    (_, Hh), th = timed(() -> solve(c.z_b, c.tau, c.mask; kw..., mode = :full,
                                     relax = c.relax_full, n_outer = c.n_outer_full), 3)
    @printf("  flat %.4f s, full %.4f s\n", tf, th)
    Hs, ts = fill(NaN, size(c.z_b)), NaN
    if c.do_smb
        (_, Hs), ts = timed(() -> diva_reconstruct(c.z_b, c.mask, c.smb, c.dx, c.dy;
                                                    rheology = c.rheo, sliding = c.slide), 1)
        @printf("  smb  %.4f s\n", ts)
    end
    path = joinpath(OUT, "julia_$name.nc")
    isfile(path) && rm(path)
    NCDataset(path, "c") do ds
        defDim(ds, "x", size(c.z_b, 1)); defDim(ds, "y", size(c.z_b, 2))
        defVar(ds, "H_flat", Hf, ("x", "y"))
        defVar(ds, "H_full", Hh, ("x", "y"))
        defVar(ds, "H_smb", Hs, ("x", "y"))
        ds.attrib["t_flat"] = tf; ds.attrib["t_full"] = th; ds.attrib["t_smb"] = ts
        ds.attrib["threads"] = Threads.nthreads()
    end
end

# JIT warm-up on a small case, so the timings exclude compilation.
let c = synthetic(20)
    solve(c.z_b, c.tau, c.mask; dx = c.dx, dy = c.dy, mode = :flat)
    solve(c.z_b, c.tau, c.mask; dx = c.dx, dy = c.dy, mode = :full, n_outer = 2)
    diva_reconstruct(c.z_b, c.mask, c.smb, c.dx, c.dy; rheology = c.rheo, sliding = c.slide, n_outer = 2)
    g = greenland()
    diva_reconstruct(g.z_b, g.mask, g.smb, g.dx, g.dy; rheology = g.rheo, sliding = g.slide, n_outer = 2)
end

run_case("synth150", synthetic(150))
run_case("synth300", synthetic(300))
run_case("grl32", greenland())
println("julia_ref: done -> $OUT")
