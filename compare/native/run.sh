#!/bin/bash
#SBATCH --partition=shared
#SBATCH --account=ba1442
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --job-name=plastic-native
#SBATCH --output=compare/out/native-%j.log
# From the repository root (after make compare): sbatch compare/native/run.sh
# Margin at the ice-free cell centre, at the face, and from the signed distance to the 5 km margin.
JULIA="${JULIA:-$HOME/.juliaup/bin/julia} --project=${JULIA_ENV:-compare/out/env}"
X=./libplastic/bin/plastic_native.x
export OUT=${OUT:-/work/ba1442/robinson/PlasticIceSheet-fortran/native}
for t in ${TIMES:-0 20000}; do for f in ${FS:-1 2 4 6}; do
  $JULIA compare/native/native.jl prep $t $f
  for mg in centre face sdist; do
    $X $OUT/t${t}_f${f} surface $mg 917.0 9.80665
    $JULIA compare/native/native.jl cmp $t $f surface_$mg
  done
done; done
