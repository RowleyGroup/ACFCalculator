#include <boost/program_options.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/tokenizer.hpp>

#ifdef FFTW
#include <complex.h>
#include <fftw3.h>
#endif

using boost::math::policies::policy;
using boost::math::tools::newton_raphson_iterate;
using boost::math::tools::halley_iterate;
using boost::math::tools::eps_tolerance; // Binary functor for specified number of bits.
using boost::math::tools::bracket_and_solve_root;
using boost::math::tools::toms748_solve;

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <regex>
#include <vector>

const double seg_min=0.2;
const double s_start=0.000001;
const double s_end=1.0;
const double s_increment=0.0001;

const double fit_width1=1000.0;
const double fit_width2=1000.0;

const double kB=1.38064852E-23;
const double default_temperature=298.15;

const int default_maxcorr=1000;
const double default_timestep=1.0;
const double R=8.314;

// thresholds for warnings about decay of ACF and VACF
const double acf_warning=0.02;
const double vacf_warning=0.02;

// compile: mkdir build  cd build  cmake ..  make

// syntax: ./ACFCalculator -i[infile] -t[type] -a[acf_outfile] -o[outfile] -f[field] -s[timestep] -r[ts_factor]

// the time series can be in the format saved by NAMD solvars module, GROMACS, CHARMM, or others (use general mode)

// calculate diffusion coefficient of a harmonically-restrained degree of freedom by 
// calculating the autocorrelation function a time series of a molecular dynamics 
// simulation then integrating it 

// Hummer, G. Position-dependent diffusion coefficients and free energies from Bayesian
// analysis of equilibrium and replica molecular dynamics simulations. New J. Phys. 2005,
// 7, 34. 

// internal units:
// distance = angstrom
// time = femptoseconds
// velocity = Angstrom / femptoseconds


namespace po = boost::program_options;

double *vacf;
int nCorr;
double var;
double varVel;
double timestep;

// calculate laplace transform of a series at time s
double laplace(double *series, double deltat, double s, int length)
{
  double F=0.0;

  for(int i=0;i<length;++i)
    {
      F+=exp(-s*i*deltat)*series[i]*deltat;
    }

  return(F);
}

// calculate denominator of D(s) function
double denom_func(double s)
{
  double laplaceVACF=laplace(vacf, timestep, s, nCorr);
  double f=laplaceVACF*(s*var+varVel/s)-var*varVel;
  return(f);
}

// calculate D(s) function
double Ds_func(double s)
{
  double laplaceVACF=laplace(vacf, timestep, s, nCorr);
  double f=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);
  return(f);
}

// Calculate correlation function by direct summation
//
// Allen, M.; Tildesley, D. Computer Simulation of Liquids; Oxford Science Publications, Clarendon Press: Oxford, 1989. 
// @y time series
// @nSamples number of data points in time series
// @nCorr length of correlation function to calculate
// @return pointer to double array containing correlation function

double *calcCorrelation(double *y, int nSamples, int nCorr)
{
  double *corr=new double[nCorr];
  int t;
  int ttoMax;

  for(int i=0;i<nCorr;++i)
    {
      corr[i]=0.0;
    }

  for(int i=0;i<nSamples;++i)
    {
      ttoMax=nSamples;

      if(i+nCorr<nSamples)
        ttoMax=i+nCorr;

      for(int j=i;j<ttoMax;++j)
        {
          t=j-i;
          corr[t]+=y[i]*y[j];
        }
    }

  for(int i=0;i<nCorr;++i)
    {
      corr[i]=corr[i]/(nSamples-i);
    }

  return(corr);
}

#ifdef FFTW
/*
//Work in progress
double *calcCorrelation_FFT(double *y, int nSamples, int nCorr)
{
  double *corr=new double[nCorr];
  fftw_complex *out;
  fftw_plan p;

  out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nSamples);
  
  p = fftw_plan_dft_r2c_1d(nSamples, y, out, FFTW_ESTIMATE);

  // calculate fft
  
  fftw_execute(p);
  
  // calculate complex conjugate

  fftw_plan p_inv;
  
  double *cc_in=new double[nSamples];
    
  for(int i=0;i<nSamples;++i)
    {
      cc_in[i]=out[i][0]*out[i][0]+out[i][1]*out[i][1];
      std::cout << i << " " << cc_in[i] << std::endl;
    }

  fftw_destroy_plan(p);
  
  //  fftw_complex *cc_out=(fftw_complex *) fftw_malloc(sizeof(fftw_complex) * nSamples);
  double *cc_out=new double[nSamples];
  //calculate inverse fft
  p_inv=fftw_plan_r2r_1d(nSamples, cc_in, cc_out, FFTW_REDFT00, FFTW_ESTIMATE);

  fftw_execute(p_inv);

  // copy correlation function into output and normalize
  for(int i=0;i<nCorr;++i)
    {
      corr[i]=cc_out[i]/nSamples;
      std::cout << i << " " << corr[i] << std::endl;
    }
  fftw_destroy_plan(p_inv);
  
  fftw_free(cc_out);
  
  delete[] cc_in;
  return(corr);
}
*/
#endif

// calculate variance in series
double variance(double *y, int nSamples)
{
  double v2=0.0;
  for(int i=0;i<nSamples;++i)
    v2+=y[i]*y[i];
  v2/=nSamples;
  return(v2);
}

// calculates average and subtracts it from series, returns average
double calc_subtract_average(double *y, int nSamples)
{
  double avg=0.0;

  for(int i=0;i<nSamples;++i)
    {
      avg+=y[i];
    }
  avg/=nSamples;

  for(int j=0; j<nSamples; ++j)
    {
      y[j]-=avg;
    }
  return(avg);
}

// calculate the velocity time series by numerically differentiating the position time series
double *velocity_series(double *y, int nSamples, double deltat)
{
  double *vel=new double[nSamples-2];
  double avg=0.0;

  // calculate velocity by finite difference approach. units: A/fs
  for(int i=1;i<nSamples-1;++i)
    {
      vel[i-1]=(y[i+1]-y[i-1])/(2.0*deltat);
    }

  avg=calc_subtract_average(vel, nSamples-2);
  return(vel);
}

// numerically integrates correlation series upto given cutoff
double integrateCorrCutoff(double *acf, int nCorr, double deltat,  double cutoff)
{
  double acfInt=0.0;

  for(int i=0;i<nCorr-1;++i)
    {
      if(cutoff!=0)
        if(acf[i]<cutoff*acf[0])
	  break;
      acfInt+=0.5*(acf[i]+acf[i+1])*deltat;
    }
  return(acfInt);
}

// read time series from filename fname
// file is formed like NAMD colvar traj file
// store the number of points in numSamples
// each line should store one time step
std::vector<double> readSeriesNAMD(char *fname, int &numSamples, int field)
{
  std::ifstream datafile(fname, std::ifstream::in);
  std::vector<double> series;

  std::string line;
  std::istringstream iss;
  int begin;

  if(field==1)
    begin=15;
  else if(field==2)
    begin=37;
  else
    begin=61;

  numSamples=0;

  while(getline(datafile,line))
    {
      if(line.at(0)!='#')
        {
          try
            {
              std::string str2=line.substr(begin,23);
              series.push_back(atof(str2.c_str()));
              ++numSamples;
            }
          catch (std::exception& e)
            {
              break;
            }
        }
    }

  return(series);
}

// read time series from CHARMM file fname
// store the number of points in numSamples
// each line should store one time step
std::vector<double> readSeriesCHARMM(char *fname, int &numSamples)
{
  std::ifstream datafile(fname, std::ifstream::in);
  std::vector<double> series;

  std::string line;
  std::istringstream iss;

  numSamples=0;

  while(getline(datafile,line))
    {
      if(line.at(0)!='#')
        {
          try
            {
              series.push_back(atof(line.c_str()));
              ++numSamples;
            }
          catch (std::exception& e)
            {
              break;
            }
        }
    }

  return(series);
}

// read time series from GROMACS file fname
// store the number of points in numSamples
// each line should store one time step
std::vector<double> readSeriesGROMACS(char *fname, int &numSamples, int field)
{
  std::ifstream datafile(fname, std::ifstream::in);
  std::vector<double> series;
  std::string line;
  std::string item;
  std::istringstream iss;
  
  numSamples=0;
  
  while(getline(datafile,line))
    {
      if(line.at(0)!='#'&& line.at(0)!='@')
	{
	  try
	    {
	      double dbl;
	      int nchar=0;
	      char *cline=(char *) line.c_str();
	      
	      for(int i=1;i<=field;++i)
		{
		  if(sscanf(cline, "%lf%n", &dbl, &nchar) == 0)
		    break;
		  cline+=nchar;
		}
	      // convert nm to Angstrom
	      series.push_back(dbl*10.0);
	      ++numSamples;
	    }
	  catch (std::exception& e)
	    {
	      break;
	    }
	}
    }

  return(series);
}


// read time series from filename fname
// works for NAMD, GROMACS, CHARMM, and other formated files given correct field
// store the number of points in numSamples
// each line should store one time step
// factor for converting time to fs
std::vector<double> readSeriesGeneral(char *fname, int &numSamples, int field, double ts_factor, double p_factor)
{
  std::ifstream datafile(fname, std::ifstream::in);
  std::vector<double> series;
  std::string line;

  std::regex e ("[-]?[0-9]*\\.?[0-9]+");
  boost::char_separator<char> sep{" ", "\t"};

  numSamples=0;
  while(getline(datafile,line))
    {
      // tokenizes a line by stripping spaces and tabs into traj vector
      boost::tokenizer<boost::char_separator<char>> tok(line, sep);
      std::vector<std::string> traj;

      for(const auto& t:tok)
  	{
	  traj.push_back(t);
	}

      // adds value in field to series if start of line is numeric    
      if(std::regex_match(traj[0], e))
        { 
          try
            {
              series.push_back(atof(traj[field].c_str())*p_factor);
              ++numSamples;
            }
           catch (std::exception& e)
            {
              break;
            }
        }
        
    }
  return(series);
}

// calculates intercept by linear interpolation for subset of D(s)
double leastsquares_subset(const std::vector<double>& x, const  std::vector<double>& y, int lower, int upper, double &m, double &b, double &r2)
{
  int n=upper-lower;
  double s_x  = accumulate(x.begin()+lower, x.begin()+upper, 0.0);
  double s_y  = accumulate(y.begin()+lower, y.begin()+upper, 0.0);

  double s_xx = inner_product(x.begin()+lower, x.begin()+upper, x.begin()+lower, 0.0);
  double s_xy = inner_product(x.begin()+lower, x.begin()+upper, y.begin()+lower, 0.0);
  double s_yy = inner_product(y.begin()+lower, y.begin()+upper, y.begin()+lower, 0.0);
 
  m=(n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
  b=(s_y-m*s_x)/n;

  // calculate residuals, use to set r2
  double resSum=0.0;
  
  for(int i=lower;i<upper;++i)
    {
      double res=y.at(i)-(m*x.at(i)+b);
      resSum+=res*res;
    }
  r2=1.0-resSum/s_yy;

  return(b); 
}

// searches backwards from second singularity until line tangent
// to the curve would give a non-negative intercept
int find_negative_bound(const std::vector<double>& x, const std::vector<double>& y)
{
  for(unsigned int i=x.size()-2;i>0;--i)
    {
      double slope=(y[i+1]-y[i])/(x[i+1]-x[i]);
      double intercept=y[i]-slope*x[i];

      if(intercept>0)
	{
	  return(i);
	}
    }
  // search  didn't work for some  reason, do not shrink area
  return(x.size());
}

// calculate difference in inital slope
std::vector<double> diff(const std::vector<double>& x, const std::vector<double>& y)
{
  std::vector<double> df;

  double init_slope1=(y.at(1)-y.at(0))/(x.at(1)-x.at(0));
  double init_slope2=(y.at(2)-y.at(1))/(x.at(2)-x.at(1));
  
  df.push_back( init_slope1-(init_slope2-init_slope1));
  
  for(unsigned int i=0;i<y.size()-1;++i)
    {
      df.push_back( (y.at(i+1)-y.at(i))/(x.at(i+1)-x.at(i)) );
    }

  return(df);
}

// finds the minimum
int findmin(const std::vector<double>& y)
{
  int min=0;
  double minval=fabs(y.at(0));
  
  for(unsigned int i=0;i<y.size();++i)
    {
      if(fabs(y.at(i))<minval)
	{
	  min=i;
	  minval=fabs(y.at(i));
	}
    }
  return(min);
}

// return 5 pt running average
std::vector<double> smooth(const std::vector<double>& y)
{
  std::vector<double> ySmooth;
  ySmooth.push_back( (y.at(0) + y.at(1) + y.at(2))/3.0);
  ySmooth.push_back( (y.at(0) + y.at(1) + y.at(2) + y.at(3))/4.0);
	      
  for(unsigned int i=2;i<y.size()-2;++i)
    {
      ySmooth.push_back( (y.at(i-2) + y.at(i-1) + y.at(i) + y.at(i+1) + y.at(i+2))/5.0);
    }
  ySmooth.push_back(y.at(y.size()-1));
  ySmooth.push_back(y.at(y.size()-2));
  return(ySmooth);
}

// calculates Ds from s to s_max
void calcDs(double s, double s_max, double s_delta, std::vector<double>& s_vals, std::vector<double>& ds_vals)
{
  while(s<=s_max)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double denom=(laplaceVACF*(s*var+varVel/s)-var*varVel);
      double Ds_val=-(laplaceVACF*var*varVel)/denom;

      // adds if positive
      if(Ds_val>0.0)
        {
          s_vals.push_back(s);
          ds_vals.push_back(Ds_val);
        }
      s=s+s_delta;
    }
}

int main(int argc, char *argv[])
{
  enum InputType {gromacs, namd, charmm, general};
  InputType type;
  
  bool write_acf;
  std::vector<double> series, seriesVel;
  double *velSeries;
  double *acf, *timeSeries;
  double acfInt;
  char *fname, *acf_fname, *output_fname;
  int field=1;
  int numSamples;
  double cutoff;
  int bits = std::numeric_limits<double>::digits;
  double s1, s2;
  double ts_factor, p_factor;
  double k=0.0;
  double temperature=default_temperature;
  
  double mass=1.0;
  
  po::options_description desc("Allowed options", 1024, 512);

  desc.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<std::string>()->required(), "file name of time series")
    ("type,t", po::value<std::string>(), "type of time series (namd, charmm, gromacs, general)")
    ("cutoff,c", po::value<double>()->default_value(0), "cutoff to integrate ACF")
    ("acf,a", po::value<std::string>()->required(), "file name to save autocorrelation functions in")
    ("output,o", po::value<std::string>()->required(), "file name to save output to")
    ("timestep,s", po::value<double>(&timestep)->default_value(default_timestep), "time between samples in time series file (fs)")
    ("maxcorr,m", po::value<int>(&nCorr)->default_value(default_maxcorr), "maximum number of time steps to calculate correlation functions over")
    ("field,f", po::value<int>(&field)->default_value(1), "index of field to read from time series file")
    ("ts_factor,r", po::value<double>(&ts_factor)->default_value(1), "factor for time conversion to fs when using general, default 1. ex: for ps to fs use 10")
    ("p_factor,p", po::value<double>(&p_factor)->default_value(1), "factor for position to Angstrom when using general, default 1. ex: for nm to Angstrom use 10")
    ("spring,k", po::value<double>(&k)->default_value(0), "spring constant of harmonic restraint (units:kcal A-2 /mol, optional)") 
    ("mass,w", po::value<double>(&mass)->default_value(0), "mass of the solute (units amu, optional)") 
    ("temperature,e", po::value<double>(&temperature)->default_value(0), "temperature (optional)")  
#ifdef FFTW
    ("fft", po::value<std::string>(), "use fft to calculate correlation functions")
#endif
    ;

  try
    {
      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);

      if (vm.count("help"))
        {
	  std::cout << "Usage: ACFcalculator [options]" << std::endl;
	  std::cout << desc;
          return 1;
        }
      po::notify(vm);
      
      if (vm.count("input"))
        {
          const std::string fname_str=vm["input"].as<std::string>();
          fname = new char [fname_str.length()+1];
          std::strcpy(fname, fname_str.c_str());
        }
      else
        {
          std::cout << "Time series file must be provided" << std::endl;
          return(1);
        }

     if(vm.count("type"))
        {
          const std::string type_str=vm["type"].as<std::string>();
	  if(type_str.compare("gromacs")==0)
	    {
	      type=gromacs;
	    }
	  else if(type_str.compare("charmm")==0)
	    {
	      type=charmm;
	    }
	  else if(type_str.compare("namd")==0)
	    {
	      type=namd;
	    }
          else if(type_str.compare("general")==0)                                                           
            { 
              type=general;
            }
	  else
	    {
	      std::cout<<"Type of file not recognized"<<std::endl;
	      return(1);
	    }
	}
     else
       {
	 std::cout << "#Type of file not provided, defaulting to general file reader" << std::endl;
	 type=general;
       }

      if (vm.count("acf"))
        {
          const std::string acf_fname_str=vm["acf"].as<std::string>();
          acf_fname = new char [acf_fname_str.length()+1];
          std::strcpy(acf_fname, acf_fname_str.c_str());
          write_acf=true;
        }
      else
        {
          std::cout << "#ACF will not be saved" << std::endl;
          write_acf=false;
        }

#ifdef FFTW
      bool use_fft=false;
      if (vm.count("fft"))
	{
	  std::cout << "#Correlation functions will be calculated using FFT" << std::endl;
	  use_fft=true;
	}      
#endif
      std::cout << "#Time step of " << timestep << " fs will be used." << std::endl;

      if(vm.count("cutoff"))
	{
          cutoff=vm["cutoff"].as<double>();
	  if(cutoff>0.0)
	    std::cout<<"#ACF cutoff of " << cutoff << " will be used." << std::endl;
	  else
	    {
	      std::cout<<"#No cutoff will be used to integrate the ACF." << std::endl;
	      cutoff=0;
	    }
	}
      else
	{
          std::cout<<"#No cutoff will be used in integraction of ACF." << std::endl;
	}
      
      if (vm.count("output"))
        {
          const std::string output_fname_str=vm["output"].as<std::string>();
          output_fname = new char [output_fname_str.length()+1];
          std::strcpy(output_fname, output_fname_str.c_str());
          std::cout << "#D(s) will be saved in " << output_fname << std::endl;
        }
      else
        {
          std::cout << "#D(s) will not be saved" << std::endl;
        }
    }
  catch ( const std::exception& e )
    {
      std::cerr << e.what() << std::endl;
      return 1;
    }


  if (type==gromacs)                                                                                
    {                                                                                                    
      series=readSeriesGROMACS(fname, numSamples, field);                                                
    }                                                                                                    
  else if (type==charmm)                                                                                 
    {
      series=readSeriesCHARMM(fname, numSamples);                                                        
    } 
  else
    { 
      std::cout<<"File type could not be determined from input, using General with time series factor of "<<ts_factor <<std::endl;
      series=readSeriesGeneral(fname, numSamples, field, ts_factor, p_factor);                                                
    }

  if(series.size()==0)
    {
      std::cerr << "Error: Time series could not be read from file." << std::endl;
      return(1);
    }
  
  timeSeries=&series[0];

  double avg=calc_subtract_average(timeSeries, numSamples);

  velSeries=velocity_series(timeSeries, numSamples, timestep);

  numSamples=numSamples-2;

  acf=calcCorrelation(timeSeries, numSamples, nCorr);
  //  calcCorrelation_FFT(timeSeries, numSamples, nCorr);
  vacf=calcCorrelation(velSeries, numSamples, nCorr);

  var=variance(timeSeries, numSamples);
  varVel=variance(velSeries, numSamples);
  std::cout << "#var " << var << std::endl;
  std::cout << "#varVel " << varVel << std::endl;

  // if a spring constant is defined, use it to calculate the variance
  if(k>0.0)
    {
      double factor;
      
      // convert mass to kg
      mass=1.66054e-27*mass;
      k=k*4184.0;
      
      var=R*temperature/k;
      
      factor=var/acf[0];
      for(int i=0;i<nCorr;++i)
        acf[i]*=factor;
      // units of A2/fs2
      
      varVel=1.0E-10*kB*temperature/mass;
      factor=varVel/vacf[0];
      
      for(int i=0;i<nCorr;++i)
        vacf[i]*=factor;
      std::cout << "#varNew " << var << std::endl;
      std::cout << "#varVelNew " << varVel << std::endl;
  
    }
  if(acf[nCorr-1]/var > acf_warning)
    {
      std::cout << "#WARNING: ACF has only decayed to " << acf[nCorr-1]/var*100 << "% of initial value at the end of the region of integration. " <<  std::endl;
    }
 
  if(vacf[nCorr-1]/var > vacf_warning)
    {
      std::cout << "#WARNING: VACF has only decayed to " << vacf[nCorr-1]/var*100 << "% of initial value at the end of the region of integration. " <<  std::endl;
    }

  // writes ACF and VACF to acf_outfile
  if(write_acf)
    {
      std::ofstream acf_file(acf_fname);

      for(int i=0;i<nCorr;++i)
        acf_file << i << " " << acf[i] << " " << vacf[i] << std::endl;

      acf_file.close();
    }

  acfInt=integrateCorrCutoff(acf, nCorr, timestep, cutoff);
  std::ofstream out_file(output_fname,  std::ofstream::out);
  
  // Step 1: find sigularities by root finding algorithms
  double s=s_increment;
  boost::uintmax_t max_iter=500;
  boost::math::tools::eps_tolerance<double> tol(30);
  
  // find minimum in denominator
  std::pair<double, double> denom_min = boost::math::tools::brent_find_minima(denom_func, s_start, s_end, bits);
  std::pair<double, double> singularity_one = boost::math::tools::toms748_solve(denom_func, s_start, denom_min.first, tol, max_iter);
  std::pair<double, double> singularity_two = boost::math::tools::toms748_solve(denom_func, denom_min.first, s_end, tol, max_iter);


  // using the singularities as the bounds for the root finding can cause numerical instabilities or
  // the minimization to search in the wrong direction
  // this section adjusts s1 and s2 so that they bracket the minimum on the negative side

  // set a tolerance 
  double denom_tol=denom_min.second/100;

  // contract range so that D(s) is positive range are between roots
  s1=singularity_one.first;

  while(denom_func(s1)>denom_tol)
    {
      s1+=s_increment;
    }

  s2=singularity_two.first;
  while(denom_func(s2)>denom_tol)
    {
      s2-=s_increment;
    }
  
  std::pair<double, double> Ds_min = boost::math::tools::brent_find_minima(Ds_func, s1, s2, bits);

  std::cout<<"#s1 "<<s1<<std::endl;
  std::cout<<"#s2 "<<s2<<std::endl;
  std::cout<<"#Ds_min "<<Ds_min.second<<std::endl;

  // Step 2: calculate s over range
  double s_range=s2-Ds_min.first;
  double s_max=s2;
  
  if(s_range<=0.0)
    {
      throw std::out_of_range("No minimum found within singularities");
    }

  double s_delta=s_range/fit_width1;
  
  std::vector<double> s_values;
  std::vector<double> intDs;
 
  calcDs(s, s_max, s_delta, s_values, intDs); 

  // Step 3: smoothening
  int upper=find_negative_bound(s_values, intDs);
  int lower=0;
  double m, b, r2;
  
  std::vector<double> df=diff(s_values, intDs);
  std::vector<double> df_smooth=smooth(df);
  
  std::vector<double> ddf=diff(s_values, df);
  std::vector<double> ddf_smooth=smooth(ddf);

  std::vector<double> s_ddf_smooth_positive;
  std::vector<double> ddf_smooth_positive;
  
  for(unsigned int i=0;i<ddf_smooth.size();++i)
    {
      if(ddf_smooth.at(i)>0.0)
	{
	  s_ddf_smooth_positive.push_back(s_values.at(i));
	  ddf_smooth_positive.push_back(ddf_smooth.at(i));
	}
    }
  
  int min_ddf=findmin(ddf_smooth_positive);  
  lower=min_ddf-1;
  double ddf_min=ddf_smooth_positive.at(min_ddf);  
 
  while(lower>0 && ddf_smooth_positive.at(lower)<ddf_smooth_positive.at(min_ddf)*2.0)
    {
      --lower;
    }

  upper=min_ddf+1;
  while(ddf_smooth_positive.at(upper) < ddf_smooth_positive.at(min_ddf)*3.0 && upper < ddf_smooth_positive.size()-1)
    {
      ++upper;
    }
  
  // ensure range is sufficiently wide
  while((upper-lower)<seg_min*fit_width1)
    {
      if(ddf_smooth_positive.at(upper) < ddf_smooth_positive.at(lower) && lower > 0)
	{
	  --lower;
	}
      else if(upper<ddf_smooth_positive.size())
	{
	  ++upper;
	}
      else
	{
	  break;
	}
    }
  
  s=s_ddf_smooth_positive.at(lower);
  s_max=s_ddf_smooth_positive.at(upper);
  s_delta=(s_max-s)/fit_width2;
  
  std::vector<double> s_final;
  std::vector<double> ds_final;

  calcDs(s, s_max, s_delta, s_final, ds_final);
 
  // Step 4: linear fit over range
  double Ds_intercept=0.0;
  
  leastsquares_subset(s_final, ds_final, 0, s_final.size(), m, Ds_intercept, r2);
  out_file << "#m= " << m << std::endl; 
  out_file << "#b= " << Ds_intercept << std::endl; 
  out_file << "#r2= " << r2 << std::endl; 
  out_file << "#ddDs_min = " << ddf_min << std::endl; 
  out_file << "#avg = " << avg << std::endl;
  out_file << "#D (ACF) = " << var*var/acfInt*0.1 << " cm2/s " << std::endl;
  out_file << "#Ds (VACF) = " << Ds_intercept*0.1 << " cm2/s " << std::endl;

  out_file.close();

  delete[] acf;
  delete[] vacf;
}
