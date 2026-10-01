#!/bin/bash

for p1 in 0 0.2 0.4 0.6 0.8 1; do # prop DE gs
  for p2 in 0 0.4 0.8; do # corr
    for p3 in 100; do # gene set size
      for p4 in 0.5; do # prop DE mat gene expression
      sbatch citcdf_gsa/comp_sea_citcdf/new_all_ok/code/new_set_comp_DE.sh $p1 $p2 $p3 $p4
      done
    done
  done
done

