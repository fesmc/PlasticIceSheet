#!/bin/bash
#SBATCH --job-name=plastic-compare
#SBATCH --account=ba1442
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=compare/out/run_%j.log
#
# Fortran port vs PlasticIceSheet.jl on Levante. From the repository root:
#   mkdir -p compare/out && sbatch compare/run.sh
# JULIA_PKG: a checkout of PlasticIceSheet.jl; the Julia environment (the package +
# NCDatasets) is created in compare/out/env on the first run.
set -e
JULIA=${JULIA:-$HOME/.juliaup/bin/julia}
JULIA_PKG=${JULIA_PKG:-/work/ba1442/robinson/PlasticIceSheet-fortran/julia-ref/PlasticIceSheet.jl}
OUT=compare/out
ENV_DIR=$OUT/env
CASES="synth150 synth300 grl32"

mkdir -p $OUT
if [ ! -f $ENV_DIR/Manifest.toml ]; then
    JULIA_PKG=$JULIA_PKG $JULIA --project=$ENV_DIR -e \
        'using Pkg; Pkg.develop(path=ENV["JULIA_PKG"]); Pkg.add("NCDatasets"); Pkg.precompile()'
fi

make compare config=${config:-dkrz_levante_ifx}
$JULIA --project=$ENV_DIR compare/julia_ref.jl $OUT
./libplastic/bin/plastic_compare.x $OUT $CASES
$JULIA --project=$ENV_DIR compare/compare.jl $OUT $CASES
