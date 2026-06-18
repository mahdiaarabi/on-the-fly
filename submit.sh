#!/bin/bash
#SBATCH -J iron
#SBATCH -o iron.o%j
#SBATCH -e iron.e%j

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A CHE21041
#SBATCH -p small
#SBATCH -t 48:00:00

export OMP_NUM_THREADS=56

rm foundling.txt
module purge
module use /work2/01114/jfonner/frontera/modulefiles
module load gaussian/16rC.01

./a.out &
sleep 10
for i in {1..30000}
  do
    g16 gin.gjf
    sleep 10
  done

