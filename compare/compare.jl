# Compare the Fortran port with PlasticIceSheet.jl: max / rms |dH| over the ice mask
# for the flat, full and SMB-driven solutions, and the wall-clock times. Writes the
# table to out/compare.md.
#
# Run: julia --project=<env with NCDatasets> compare/compare.jl [outdir] case...

using NCDatasets
using Printf

const OUT = ARGS[1]
const CASES = ARGS[2:end]

io = IOBuffer()
println(io, "| case | field | max H [m] | max abs dH [m] | rms dH [m] | t Julia [s] | t Fortran [s] | Julia / Fortran |")
println(io, "|---|---|---|---|---|---|---|---|")
for c in CASES
    mask = NCDataset(joinpath(OUT, "in_$c.nc")) do ds
        Array(ds["mask"][:, :]) .== 1
    end
    rj = NCDataset(joinpath(OUT, "julia_$c.nc"))
    rf = NCDataset(joinpath(OUT, "fortran_$c.nc"))
    for (f, t) in (("H_flat", "t_flat"), ("H_full", "t_full"), ("H_smb", "t_smb"))
        tj = rj.attrib[t]
        isfinite(tj) || continue
        Hj = Array(rj[f][:, :]); Hf = Array(rf[f][:, :])
        d = (Hf .- Hj)[mask]
        tf = rf.attrib[t]
        @printf(io, "| %s | %s | %.1f | %.3g | %.3g | %.4f | %.4f | %.1f |\n", c, f, maximum(Hj[mask]),
                maximum(abs, d), sqrt(sum(abs2, d) / length(d)), tj, tf, tj / tf)
    end
    println(io, "")
    @printf(io, "<!-- %s: Fortran n_outer full %d (converged %d), smb %d (converged %d); warm full %d vs cold %d (%.4f vs %.4f s); warm smb %d vs cold %d (%.4f vs %.4f s) -->\n",
            c, rf.attrib["n_outer_full"], rf.attrib["converged_full"], rf.attrib["n_outer_smb"],
            rf.attrib["converged_smb"], rf.attrib["warm_full_n_outer"], rf.attrib["cold_full_n_outer"],
            rf.attrib["warm_full_t"], rf.attrib["cold_full_t"], rf.attrib["warm_smb_n_outer"],
            rf.attrib["cold_smb_n_outer"], rf.attrib["warm_smb_t"], rf.attrib["cold_smb_t"])
    close(rj); close(rf)
end
s = String(take!(io))
print(s)
write(joinpath(OUT, "compare.md"), s)
