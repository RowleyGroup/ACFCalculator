# ACFCalculator
Usage: ACFcalculator [options]
##Allowed Options
```
  -h [ --help ]                produce help message
  -i [ --input ] arg           file name of time series
  -t [ --type ] arg            type of time series (namd, charmm, gromacs, general)
  -c [ --cutoff ] arg (=0)     cutoff to integrate ACF
  -a [ --acf ] arg             file name to save autocorrelation functions in
  -o [ --output ] arg          file name to save output to
  -s [ --timestep ] arg (=1)   time between samples in time series file (fs)
  -m [ --maxcorr ] arg (=1000) maximum number of time steps to calculate correlation functions over
  -f [ --field ] arg (=1)      index of field to read from time series file
  -r [ --factor] arg (=1)      factor for converting time to fs when using --type general
```

##Requirements

Boost C++ Libraries (http://www.boost.org/)

##Compilation
```
mkdir build
cd build
cmake ..
make
```

##References
*Gaalswyk, K., Awoonor-Williams, E., Rowley, C. N. Generalized Langevin Methods for Calculating Transmembrane Diffusivity, J. Chem. Theory Comput. 2016, doi: [10.1021/acs.jctc.6b00747](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00747)

*Awoonor-Williams, E., Rowley, C.N. Molecular simulation of nonfacilitated membrane permeation, Biochim. Biophys. Acta - Biomembranes 2016, doi: [10.1016/j.bbamem.2015.12.014](https://www.sciencedirect.com/science/article/pii/S0005273615004125)
