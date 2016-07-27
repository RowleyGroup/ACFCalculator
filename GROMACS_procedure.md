The following is a procedure for performing GROMACS simulations and ACFCalculator calculations on a system composed of a DPPC
membrane with a harmonically restrained solute. 

###GROMACS Simulations

All the necessary files for running a GROMACS simulation are in the `permeation/` directory which is set up as:

```
permeation/  
  down/  
    sim_0/  
      /* GROMACS input files (ex: .mdp, .gro, ...)  
    sim_1/ 
    .... 
    sim_39/  
    makemdp.py  
    gromacssubmit-list
    pull-template.mdp  
  down-nve/  
    sim_0/ 
    ... 
    sim_39/  
    makemdp.py  
    gromacssubmit-list  
    gromacssubmit-list-restart  
    pull-template.mdp  
  up/  
    ...  
  up-nve/  
    ...  
```    

Sample GROMACS input files are included in this Git. Note that these are not sufficent to set up the system, for that see a 
GROMACS guide. The following commands create individual `.mdp` files for each equilibration simulation with their unique 
restraint reference positions, and submits the jobs to the queue:

```
python makemdp.py
./gromacssubmit-list [job name e.g., down/up] -np [number of processors per simulation, 8] -c [number of windows, 40]
```

Run these for both the `down/` and `up/` directories. This procedure assumes use of the PBS queueing system. The 
gromacssubmit-list file will have to be modified if using a different queueing system.

After the equilibration simulations are complete (check with qstat), copy the directories:

```
cp -r up up-nve
cp -r down down-nve
```

The following commands regenerate the necessary `.mdp` files and submits the NVE jobs to the queue:

```
python makemdp.py
./gromacssubmit-list [job name e.g., down-nve/up-nve] -np [number of processors per simulation, 8] -c [number of windows, 40]
```

The NVE jobs might not complete (check using `wc -l pullx.xvg`, the files should contain 5000016 lines). If they're incomplete,
resubmit with:

```
./gromacssubmit-list-restart [job name e.g., down-nve/up-nve] -np [number of processors per simulation, 8] -c 
[number of windows, 40]
```

Once the simulations are complete, we need to copy the `.xvg` files into a separate directory to run ACFCalculator:

```
cd permeation/down-nve
mkdir xvg_files
find . -name "pullx.xvg" -exec cp --parents \{\} \xvg_files \;

cd xvg_files
rm -r down
rm -r xvg_files
```
The `find` command adds the `.xvg` files from the `down/` directory, and those it already copied to `xvg_files`, which is why 
we remove the redundant files. Repeat this process for `up-nve/`.

###ACFCalculator

```
sh setup.sh
```

The script `setup.sh` runs ACFCalculator and is located in `ACFCalculator/build/` with the ACFCalculator executable. It assumes
the input files can be found as `/${HOME}/permeation/down-nve/xvg_files/sim_${i}/pullx.xvg` for i:[0,39] (this will be true if 
following the above procedure). It will save the ACF calculations as `/${HOME}/acf_down/acf_down_sim${i}.dat` and the output 
from ACFCalculator as `/${HOME}/out_down/out_down_sim${i}.out`. Explanation of the options used in this file can be found in 
the README. ACFCalculator prints some data to screen, all of which is saved in the `.out` files.

```
python getDall.py
```

The `getDall.py` script compiles all of the VACF D(s) values from the `sim_${i}/` directories into a single `out_ds.out` file. 
The output file is a space-delineated file containing two columns: z-position (in Angstroms) and D(s) in cm^2/s. Any graphing
utility can be used to easily plot this data.
