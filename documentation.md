ACFCalculator
=============

[ACFCalculator](https://github.com/RowleyGroup/ACFCalculator) implements methods for calculating the position-dependent diffusion coefficient of a solute from an MD time series using methods based on the Generalized Langevin Equation. These methods estimate the diffusion coefficient of a solute from the time series of a solute restrained by a harmonic potential. These methods calculate these coefficients from the friction coefficient of the system, which are related to the Velocity Autocorrelation Function (VACF) and Position Autocorrelation Function (PACF). See [T. B. Woolf and B. Roux](http://pubs.acs.org/doi/abs/10.1021/ja00092
a048) and [G. Hummer](http://iopscience.iop.org/article/10.1088/1367-2630/7/1/034/meta;jsessionid=E4202C77BEAAF418645D578AA7EE8
FC5.c1.iopscience.cld.iop.org). The theoretical background and justification for these methods can be found in the following 
papers: for the VACF [T. B. Woolf](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC45285/), [M. F. Schumaker](http://www.sciencedir
ect.com/science/article/pii/S0006349500765222), [T. B. Woolf and B. Roux](http://pubs.acs.org/doi/abs/10.1021/ja00092a048), and 
for the PACF [G. Hummer](http://iopscience.iop.org/article/10.1088/1367-2630/7/1/034/meta;jsessionid=E4202C77BEAAF418645D578AA7E
E8FC5.c1.iopscience.cld.iop.org), [C. T. Lee](http://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00022).

Use of the code can be referenced/acknowledged as : Gaalswyk, K. , Awoonor-Williams, E.,  Rowley, C.N., 
[Generalized Langevin Methods for Calculating Transmembrane Diffusivity](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00747) *J. Chem. Theory. Comput.* **2016**  doi: 10.1021/acs.jctc.6b00747


Installation
------------

Untar `ACFCalculator.tar.gz`, (ex: `tar -xvzf ACFCalculator.tar.gz`) which creates the directory `ACFCalculator/` containing 
the necessary build and helpful run files. The [Boost C++ libraries](http://www.boost.org/) are required . To build:

```
		cd ACFCalculator
		mkdir build
		cd build
		cmake ..
		make
```

which will create the ACFCalculator executable.

Command Line Arguments and File Formats
---------------------------------------

### Arguments

```
	ACFCalculator -h [--help] -i [--input] -t [--type] -c [--cutoff] -a [--acf] 
	-o [--output] -s [--timestep] -m [--maxcorr] -f [--field] -r [--factor]
```

The required arguments are: `-i [--input]` and `-t [--type]` which will exit with printed error message if omitted, the rest 
are optional.

|Argument|Description|
|---|---|
|`help` | produces the help message listing all available command line arguments and their description.|
|`input` | specifies the file name of the time series ACFCalculator is going to read. If no time series file is provided, the program exits with printed error message.|
|`type` | specifies the type of time series input file. The options are: namd, charmm, gromacs, or general which handles generic delimiter-separated files. If no file type is given, the program exits with printed error message.|
|`cutoff` | specifies where to terminate integration when using the PACF method. If no cutoff is given, then the method integrates using all of the given time series.|
|`acf` | is the file name for saving the PACF and VACF of the whole time series. If no file name is given, these are not saved.
|`output` | is the file name for saving the results of ACFCalculator, which includes the position-dependent diffusion coefficient using the PACF and VACF methods, and the results of the linear least squares extrapolation method. If this is not given, these are not saved.|
|`timestep` | specifies the time between samples in the time series file in femtoseconds (fs). Default is 1.|
|`maxcorr` | specifies the maximum number of time steps to calculate the correlation functions over. A longer maxcorr (2000+) may be required to ensure complete sampling. Default is 1000.|
|`field` | specifies the field index to read from the time series file. Note that NAMD and the general type begin indexing at 0, while GROMACS begins indexing at 1. Default is 1.|
|`factor` | is used for converting time to fs when using the general type. A factor of 10 would convert from ps to fs. Default is 1.|

###File Formats

NAMD formatted files are expected to be space-delimited files with comments prefaced by a "\# ". The field argument specifies 
which column contains the time-series data. The first field contains the time in fs. In a two-column file, the second column 
(field=1) contains the position of the solute at that time. The commented preamble identifies the fields.

CHARMM formatted files are expected to contain a single field with comments prefaced by a "\# ". The field argument is 
irrelevant for CHARMM files.

GROMACS formatted files are expected to be space-delimited files with comments prefaced by a "\# '", and field-specific 
comments prefaced by a " @ ". The first field contains the time in fs. In a two-column file, the second column (field=2, note 
that with NAMD formatted-files this column is field=1) contains the position of the solute at that time. The commented preamble 
identifies the fields. Note that GROMACS specifies time in ps and position in nm. The code handles this automatically.

The general type expects files that are delimited using any combination of spaces, tabs, or " , ". It expects lines containing 
time-series data to begin with a numeric values; any line beginning with a non-numeric field is considered a comment, and 
ignored. The field argument is used, and indexing begins at 0 (in a two-column file, the second column (field=1) contains the 
position of the solute at that time). Make sure to use the factor argument to convert to fs if necessary.

The script `setup.sh` is an example of how to run the program for several sets of data using GROMACS formatted files. The 
script will have to be modified depending on how the input is saved, but it works in tandem with the 'getDall.py' script 
mentioned later. 

Output
------

ACFCalculator produces an [output].out file containing the position-dependent diffusion coefficient for the solute using both 
methods. The .out file is formatted as:

```
		#m [slope from linear least squares method]
		#b [Ds intercept from linear least squares method]
		#r2 [R^2 from linear least squares method]
		#ddDs_min [minima of the second derivative of the VACF]
		#avg [average position]
		#D (ACF) [diffusion coefficient using the PACF method]
		#Ds (VACF) [diffusion coefficient using the VACF method]
```		

The script `getDall.py` is included to compile all of the the diffusion coefficients from several simulations into a single 
file, which is particularly useful for doing position-dependent diffusivity where there can be several simulations. Some 
modification of this script may be necessary, especially if not used in conjunction with the `setup.sh` script. This script 
also prints to screen the average and standard deviation of the positional diffusion coefficients when compiling several sets 
of data.

Discussion
----------

Limitations with the individual methods have been discussed in [G. Hummer](http://iopscience.iop.org/article/10.1088/1367-2630/
7/1/034/meta;jsessionid=E4202C77BEAAF418645D578AA7EE8FC5.c1.iopscience.cld.iop.org), and [C. T. Lee](http://pubs.acs.org/doi/
abs/10.1021/acs.jcim.6b00022). Care should be taken when using the code for heterogeneous systems like a bilipid membrane. A 
un-equilibrated system may cause the code to be unable to identify minimas in the D(s) function, which will produce an error. 
Or, more likely, the code will run but will significantly over/under estimate the diffusion coefficient. This effect will be 
noted in the membrane interior. Check that the system has been fully equilibrated before using this code to ensure correct 
function.

Results are dependent on the choice of cutoff and maxcorr.


## References
* Gaalswyk, K., Awoonor-Williams, E., Rowley, C. N. Generalized Langevin Methods for Calculating Transmembrane Diffusivity, J. Chem. Theory Comput. 2016, doi: [10.1021/acs.jctc.6b00747](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00747)

* Awoonor-Williams, E., Rowley, C.N. Molecular simulation of nonfacilitated membrane permeation, Biochim. Biophys. Acta - Biomembranes 2016, doi: [10.1016/j.bbamem.2015.12.014](https://www.sciencedirect.com/science/article/pii/S0005273615004125)

