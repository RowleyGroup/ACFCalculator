#!/bin/bash

for i in `seq 0 39`
  do
    ./ACFcalculator -i ${HOME}/permeation/up-nve/xvg_files/sim_${i}/pullx.xvg -t gromacs -a ${HOME}/acf_up/acf_up_sim_${i}.dat -o ${HOME}/out_up/out_up_sim_${i}.out -f 2 -s 2 --maxcorr 2000
  done 
for k in `seq 0 39`
  do
    ./ACFcalculator -i ${HOME}/permeation/down-nve/xvg_files/sim_${k}/pullx.xvg -t gromacs -a ${HOME}/acf_down/acf_down_sim_${k}.dat -o${HOME}/out_down/out_down_sim_${k}.out -f 2 -s 2 --maxcorr 2000
  done

