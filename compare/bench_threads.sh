#!/bin/bash
#SBATCH --job-name=plastic-threads
#SBATCH --account=ba1442
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=compare/out/threads_%j.log
# Thread scaling of the eikonal (serial vs wavefront) and of the steady SMB-driven solve.
# From the repository root, after compare/run.sh made compare/out/in_*.nc:
#   sbatch compare/bench_threads.sh
set -e
make clean > /dev/null
make compare openmp=1 lis=1 > /dev/null
export OMP_PROC_BIND=close OMP_PLACES=cores
for t in 1 4 8 16; do
    OMP_NUM_THREADS=$t ./libplastic/bin/bench_threads.x compare/out
done
