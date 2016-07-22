# ACFCalculator
Usage: ACFcalculator [options]
Allowed options:

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
  

##Requirements:

Boost C++ Libraries (http://www.boost.org/)

##Compilation:

mkdir build

cd build

cmake ..

make

