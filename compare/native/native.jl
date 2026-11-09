# Native-grid check: ICESHEET (PaleoMIST 1.0) 5 km North America output vs PlasticIceSheet on the
# same inputs, coarsened by f (block means; ice where >= half the block has ice).
# Adapted from /work/ba1442/robinson/plastic-diag/scripts/native.jl (read-only inputs there).
#   native.jl prep <t yr BP> <f>        writes <OUT>/t<t>_f<f>/{hdr.txt,bed,tau,mask,Href,sdist}.bin
#   native.jl cmp  <t> <f> <label>      compares H_<label>.bin with the block-mean ICESHEET thickness
#                                       (volume ratio) and with the 5 km thickness at the coarse cell
#                                       centres (Hcen: point values, as the plastic H; ratio_centre)
# sdist: on coarse ice-free cells next to coarse ice, minus the distance from the cell centre to
# the nearest 5 km ICESHEET ice cell, less half a fine cell (the fine margin at its face).
# Bed: modern_topo - deform(t), z_sl = 0, tau_b: taub_t of the prepared product at t.
using NCDatasets, Statistics, Printf
const R = "/work/ba1442/fesmdata-sources/gowan2021_paleomist/ice_reconstruction/ice_reconstruction_files/North_America"
const T = "/work/ba1442/robinson/fesmdata-taub-prepared/Gowan2021_PaleoMIST_taub/Gowan2021_PaleoMIST_taub_NorthAmerica.nc"
const O = get(ENV, "OUT", "/work/ba1442/robinson/PlasticIceSheet-fortran/native")
rd(f) = NCDataset(ds -> Float64.(coalesce.(ds["z"][:, :], NaN)), f)
wr(f, a) = open(io -> write(io, Float64.(a)), f, "w")
rb(f, nx, ny) = (a = Array{Float64}(undef, nx, ny); open(io -> read!(io, a), f); a)
block(a, f, op) = (s = size(a) .÷ f; [op(view(a, (i-1)*f+1:i*f, (j-1)*f+1:j*f)) for i in 1:s[1], j in 1:s[2]])
function component(mask, seed)
    c = falses(size(mask)); st = [seed]; c[seed] = true
    while !isempty(st)
        p = pop!(st)
        for d in (CartesianIndex(1,0), CartesianIndex(-1,0), CartesianIndex(0,1), CartesianIndex(0,-1))
            q = p + d
            checkbounds(Bool, mask, q) && mask[q] && !c[q] && (c[q] = true; push!(st, q))
        end
    end
    c
end
function taub_at(t)
    td = NCDataset(T); dom = coalesce.(td["taub_domain"][:, :], -9999); ids = td["domain_id"][:]; tt = td["time"][:]
    tb = Float64.(td["taub_t"][:, :]); close(td)
    kt = argmin(abs.(tt .- t)); idx = Dict(ids[k] => k for k in eachindex(ids))
    [d > 0 && haskey(idx, d) ? tb[idx[d], kt] : 5000.0 for d in dom]
end
function sdist_coarse(ice, m, f)            # ice: fine mask, m: coarse mask
    nx, ny = size(m); sd = zeros(nx, ny); hf = 5000.0
    for j in 1:ny, i in 1:nx
        m[i, j] && continue
        any(m[ii, jj] for (ii, jj) in ((i-1,j), (i+1,j), (i,j-1), (i,j+1)) if 1 <= ii <= nx && 1 <= jj <= ny) || continue
        xc, yc = (i - 0.5) * f, (j - 0.5) * f          # coarse centre in fine-cell units
        dmin = Inf
        for jj in max(1, (j-2)*f+1):min(size(ice, 2), (j+1)*f), ii in max(1, (i-2)*f+1):min(size(ice, 1), (i+1)*f)
            ice[ii, jj] && (dmin = min(dmin, hypot(ii - 0.5 - xc, jj - 0.5 - yc)))
        end
        sd[i, j] = isfinite(dmin) ? -max(dmin - 0.5, 0.0) * hf : -f * hf
    end
    sd
end
cmd = ARGS[1]; t = parse(Int, ARGS[2]); f = parse(Int, ARGS[3])
D = "$O/t$(t)_f$(f)"
if cmd == "prep"
    mkpath(D)
    H = rd("$R/thickness/$t.nc"); H[isnan.(H)] .= 0; H[H .< 0] .= 0
    b = rd("$R/modern_topo/North_America.nc") .- rd("$R/deform/$t.nc"); b[isnan.(b)] .= 0
    tau = taub_at(t)
    bc = block(b, f, mean); fr = block(Float64.(H .> 0), f, mean); Hc = block(H, f, mean); tc = block(tau, f, mean)
    m = fr .>= 0.5; dx = 5000.0 * f
    sd = sdist_coarse(H .> 0, m, f)
    # 5 km thickness at the coarse cell centres (mean of the 1 or 2 x 2 central fine cells)
    c0 = isodd(f) ? (f + 1) ÷ 2 : f ÷ 2; c1 = isodd(f) ? c0 : c0 + 1
    Hcen = block(H, f, b -> mean(view(b, c0:c1, c0:c1)))
    nx, ny = size(bc)
    open("$D/hdr.txt", "w") do io; println(io, nx, " ", ny, " ", dx); end
    wr("$D/bed.bin", bc); wr("$D/tau.bin", tc); wr("$D/mask.bin", Float64.(m)); wr("$D/Href.bin", Hc); wr("$D/sdist.bin", sd); wr("$D/Hcen.bin", Hcen)
    @printf("prep t=%d f=%d: %d x %d, dx %.0f m, ice cells %d (fine %d)\n", t, f, nx, ny, dx, count(m), count(H .> 0))
else
    lbl = ARGS[4]
    nx, ny, dx = split(readline("$D/hdr.txt")); nx = parse(Int, nx); ny = parse(Int, ny); dx = parse(Float64, dx)
    Hp = rb("$D/H_$lbl.bin", nx, ny); Hr = rb("$D/Href.bin", nx, ny); m = rb("$D/mask.bin", nx, ny) .> 0.5
    Hc = rb("$D/Hcen.bin", nx, ny)
    big = component(m .| (Hr .> 0), argmax(Hr))
    d = (Hp .- Hr)[big .& (Hr .> 1)]
    @printf("RESULT t=%6d dx=%3.0f km %-16s ratio %.3f  ratio_centre %.3f  bias %6.1f rms %6.1f m (largest body)\n",
            t, dx / 1e3, lbl, sum(Hp[big]) / sum(Hr[big]), sum(Hp[big]) / sum(Hc[big]), mean(d), sqrt(mean(d .^ 2)))
end
