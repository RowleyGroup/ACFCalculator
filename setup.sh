#!/bin/bash

for j in `seq 1 5`
do
  for i in `seq 0 39`
   do
    ./ACFcalculator -i ${HOME}/permeation_${j}/up-nve/xvg_files/sim_${i}/pullx.xvg -t gromacs -a ${HOME}/acf_up${j}/acf_up_sim_${i}.dat -o ${HOME}/out_up${j}/out_up_sim_${i}.out -f 2 -s 2 --maxcorr 2000
  done 
  for k in `seq 0 39`
   do
    ./ACFcalculator -i ${HOME}/permeation_${j}/down-nve/xvg_files/sim_${k}/pullx.xvg -t gromacs -a ${HOME}/acf_down${j}/acf_down_sim_${k}.dat -o${HOME}/out_down${j}/out_down_sim_${k}.out -f 2 -s 2 --maxcorr 2000
  done
done

