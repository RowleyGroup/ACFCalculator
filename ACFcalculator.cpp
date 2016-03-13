#include <boost/program_options.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/vector.hpp>

using boost::math::policies::policy;
using boost::math::tools::newton_raphson_iterate;
using boost::math::tools::halley_iterate;
using boost::math::tools::eps_tolerance; // Binary functor for specified number of bits.
using boost::math::tools::bracket_and_solve_root;
using boost::math::tools::toms748_solve;

#include <stdexcept>

#include <iostream>
#include <fstream>

#include <vector>

#define DEFAULT_MAXCORR 1000
#define DEFAULT_TIMESTEP 1
#define SEG_MIN 0.2
#define S_INCREMENT 0.0001

#define LINEAR_BINS 10 // define range into bins

// compile: g++ correlation.cpp -o correlation 

// syntax: ./correlation time-series.traj 

// the time series is assumed to be in the format saved by the NAMD colvars module 

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

double laplace(double *series, double timestep, double s, int length)
{
  double F=0.0;

  for(int i=0;i<length;++i)
    {
      F+=exp(-s*i*timestep)*series[i]*timestep;
    }

  return(F);
}


double denom_func(double s)
{
  double laplaceVACF=laplace(vacf, timestep, s, nCorr);
  double f=laplaceVACF*(s*var+varVel/s)-var*varVel;
  return(f);
}

double Ds_func(double s)
{
  double laplaceVACF=laplace(vacf, timestep, s, nCorr);
  double f=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);
  std::cout << "s " << s << " " << laplaceVACF << " " << f << " " << var << " "<< varVel << std::endl;
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

double variance(double *y, int nSamples)
{
  double v2=0.0;
  for(int i=0;i<nSamples;++i)
    v2+=y[i]*y[i];
  v2/=nSamples;
  return(v2);
}

double calc_average(double *y, int nSamples)
{
  double avg=0.0;

  for(int i=0;i<nSamples;++i)
    {
      avg+=y[i];
    }
  avg/=nSamples;

  return(avg);
}

void subtract_average(double *y, int nSamples)
{
  double avg=calc_average(y, nSamples);

  for(int i=0;i<nSamples;++i)
    y[i]-=avg;
}

// calculate the velocity time series by numerically differentiating the position time series

double *velocity_series(double *y, int nSamples, double timestep)
{
  double *vel=new double[nSamples-2];
  double avg=0.0;
  // calculate velocity by finite difference approach. units: A/fs
  for(int i=1;i<nSamples-1;++i)
    {
      vel[i-1]=(y[i+1]-y[i-1])/(2.0*timestep);
    }

  for(int i=0;i<nSamples-2;++i)
    {
      avg+=vel[i];
    }
  avg/=(nSamples-2);

  for(int i=0;i<nSamples-2;++i)
    vel[i]-=avg;

  return(vel);
}

double integrateCorr(double *acf, int nCorr, double timestep)
{
  double I=0.0;

  for(int i=0;i<nCorr-1;++i)
    {
      I+=0.5*(acf[i]+acf[i+1])*timestep;
    }
  return(I);
}

double integrateCorrCutoff(double *acf, int nCorr, double timestep,  double cutoff)
{
  double I=0.0;

  for(int i=0;i<nCorr-1;++i)
    {
      if(acf[i]<cutoff*acf[0])
	break;
      I+=0.5*(acf[i]+acf[i+1])*timestep;
    }
  return(I);
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

// read time series from filename fname
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
// // store the number of points in numSamples
// // each line should store one time step

//  xxx change to split fields

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
      if(line.at(0)!='#'|| line.at(0)!='@')
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

// calculates intercept by linear interpolation

double leastsquares(const std::vector<double>& x, const std::vector<double>& y)
{
    unsigned int n=x.size();
    double s_x  = accumulate(x.begin(), x.begin(), 0.0);
    double s_y  = accumulate(y.begin(), y.end(), 0.0);
    double s_xx = inner_product(x.begin(), x.end(), x.begin(), 0.0);
    double s_xy = inner_product(x.begin(), x.end(), y.begin(), 0.0);
    double m=(n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    double b=(s_y-m*s_x)/n;

    return(b);
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

int main(int argc, char *argv[])
{
  enum InputType {gromacs, namd, charmm};
  InputType  type=namd;
  
  bool write_acf;
  std::vector<double> series, seriesVel;
  double *velSeries;
  double *acf, *timeSeries;
  double I;
  char *fname, *acf_fname, *output_fname;
  int field=1;
  int numSamples;
  double cutoff;
  int bits = std::numeric_limits<double>::digits;
  
  po::options_description desc("Allowed options", 1024, 512);

  desc.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<std::string>()->required(), "file name of time series")
    ("type,t", po::value<std::string>(), "type of time series (namd, charmm, gromacs)")
    ("cutoff,c", po::value<double>()->default_value(0), "cutoff to integrate ACF")
    ("acf,a", po::value<std::string>()->required(), "file name to save autocorrelation functions in")
    ("output,o", po::value<std::string>()->required(), "file name to save output to")
    ("timestep,s", po::value<double>(&timestep)->default_value(DEFAULT_TIMESTEP), "time between samples in time series file (fs)")
    ("maxcorr,m", po::value<int>(&nCorr)->default_value(DEFAULT_MAXCORR), "maximum number of time steps to calculate correlation functions over")
    ("field,f", po::value<int>(&field)->default_value(1), "index of field to read from time series file")
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
	  else
	    {
	      type=namd;
	    }
	}
     else
       {
	 std::cout << "#Type of file not provided, defaulting to NAMD" << std::endl;
	 type=namd;
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

      
      std::cout << "#Time step of " << timestep << " fs will be used." << std::endl;

      if(vm.count("cutoff"))
	{
          cutoff=vm["cutoff"].as<double>();
	  if(cutoff>0.0)
	    std::cout<<"#ACF cutoff of " << cutoff << " will be used." << std::endl;
	  else
	    std::cout<<"#No cutoff will be used to integrate the ACF." << std::endl;
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

  if(type==namd)
    { 
      series=readSeriesNAMD(fname, numSamples, field);
    }
  else if(type==gromacs)
    {
      series=readSeriesGROMACS(fname, numSamples, field);
    }
  else
    {
      series=readSeriesCHARMM(fname, numSamples);
    }

  if(series.size()==0)
    {
      std::cerr << "Error: Time series could not be read from file." << std::endl;
      return(1);
    }
  
  timeSeries=&series[0];

  double avg=calc_average(timeSeries, numSamples);
  
  subtract_average(timeSeries, numSamples);
  velSeries=velocity_series(timeSeries, numSamples, timestep);

  numSamples=numSamples-2;

  acf=calcCorrelation(timeSeries, numSamples, nCorr);
  vacf=calcCorrelation(velSeries, numSamples, nCorr);

  var=variance(timeSeries, numSamples);
  varVel=variance(velSeries, numSamples);

  if(write_acf)
    {
      std::ofstream acf_file(acf_fname);

      for(int i=0;i<nCorr;++i)
        acf_file << i << " " << acf[i] << " " << vacf[i] << std::endl;

      acf_file.close();
    }

  if(cutoff!=0)
    {    
      I=integrateCorrCutoff(acf, nCorr, timestep, cutoff);
    }
  else
    {
      I=integrateCorr(acf, nCorr, timestep);
    }
    
  std::ofstream out_file(output_fname,  std::ofstream::out);

  double s=S_INCREMENT;

  double s_min_num=1000;
  double s_min_loc=-1.0;
  
  out_file << "#" << std::setw(15) << std::left << "s"<< std::setw(15) << std::left <<"Ds" << std::setw(15) << std::left << "laplaceVACF" << std::setw(15) << std::left << "Ds numerator" << std::setw(15) << std::left << "Ds denominator" << std::endl;

  // find sigularities by root finding algorithms

  // find minimum in denominator

  boost::uintmax_t max_iter=500;
  boost::math::tools::eps_tolerance<double> tol(30);
  
  std::pair<double, double> denom_min = boost::math::tools::brent_find_minima(denom_func, S_INCREMENT, 1.0, bits);

  std::pair<double, double> singularity_one = boost::math::tools::toms748_solve(denom_func, S_INCREMENT, denom_min.first, tol, max_iter);
  std::pair<double, double> singularity_two = boost::math::tools::toms748_solve(denom_func, denom_min.first, 1.0, tol, max_iter);

  std::cout << "sigularity " << singularity_one.first <<  " " << singularity_two.first << std::endl;

  std::pair<double, double> Ds_min = boost::math::tools::brent_find_minima(Ds_func, singularity_one.first, singularity_two.first, bits);
  

  out_file << "#Ds_min " << Ds_min.first << std::endl;
  
  // Step2: calculate s over range

  double s_range=singularity_two.first-Ds_min.first;
  double s_max=singularity_two.first;

  double s_delta=s_range/1000.0;
  
  std::vector<double> s_values;
  std::vector<double> intDs;

  out_file << std::setw(15) << std::left << "#s"<< std::setw(15) << std::left <<"Ds" << std::setw(15) << std::left << "laplaceVACF" << std::setw(15) << std::left << "Ds numerator" << std::setw(15) << std::left << "Ds denominator" << std::endl;
  
  while(s<=s_max)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double Ds_val=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);
      out_file << std::setw(15) << std::left << s<< std::setw(15) << std::left <<  Ds_val << std::setw(15) << std::left << laplaceVACF << std::setw(15)<< std::left << -(laplaceVACF*var*varVel) << std::setw(15) << std::left << (laplaceVACF*(s*var+varVel/s)-var*varVel) << std::endl;
      
      if(Ds_val>0)
	{
	  s_values.push_back(s);
	  intDs.push_back(Ds_val);
	}
      s=s+s_delta;
    }
  
  int upper=find_negative_bound(s_values, intDs);
  int lower=0;
  
  // Step 3
  double m, b, r2;
  
  std::vector<double> df=diff(s_values, intDs);
  std::vector<double> df_smooth=smooth(df);
  
  std::vector<double> ddf=diff(s_values, df);
  std::vector<double> ddf_smooth=smooth(ddf);

  std::vector<double> s_ddf_smooth_positive;
  std::vector<double> ddf_smooth_positive;
  
  for(unsigned int i=0;i<ddf_smooth.size();++i)
    {
      std::cout << s_values.at(i) << " " << intDs.at(i) << " " << df.at(i) << " " << df_smooth.at(i) << " " << ddf.at(i) << " " << ddf_smooth.at(i) << std::endl;
      if(ddf_smooth.at(i)>0)
	{
	  s_ddf_smooth_positive.push_back(s_values.at(i));
	  ddf_smooth_positive.push_back(ddf_smooth.at(i));
	}
    }
  
  int min_ddf=findmin(ddf_smooth_positive);
  std::cout << "#min " << min_ddf << std::endl;
  lower=min_ddf-1;
  
  while(lower>0 && ddf_smooth_positive.at(lower)<ddf_smooth_positive.at(min_ddf)*2.0)
    {
      --lower;
    }
  std::cout << "#lower " << lower <<  std::endl;

  upper=min_ddf+1;
  while(ddf_smooth_positive.at(upper) < ddf_smooth_positive.at(min_ddf)*3.0 && upper < ddf_smooth_positive.size()-1)
    {
      ++upper;
    }
  
  std::cout << "#lower " << lower << " upper " <<  upper << std::endl;
  std::cout << "#lower_s " << s_ddf_smooth_positive.at(lower) << " upper_s " <<  s_ddf_smooth_positive.at(upper) << std::endl;

  // ensure range includes at least 50 points (xxx replace with fraction xxx)
  while((upper-lower)<50)
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
  std::cout << "#lower extend " << lower << " upper " <<  upper << std::endl;
  std::cout << "#lower_s extend " << s_ddf_smooth_positive.at(lower) << " upper_s " <<  s_ddf_smooth_positive.at(upper) << std::endl;
  
  s=s_ddf_smooth_positive.at(lower);
  s_max=s_ddf_smooth_positive.at(upper);
  s_delta=(s_max-s)/100.0;

  
  std::vector<double> s_final;
  std::vector<double> ds_final;
  
  while(s<=s_max)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double denom=(laplaceVACF*(s*var+varVel/s)-var*varVel);

      double Ds_val=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);
            
      s_final.push_back(s);
      ds_final.push_back(Ds_val);

      s=s+s_delta;
    }
  
  // linear fit over range
  double Ds_intercept=0.0;
  
  leastsquares_subset(s_final, ds_final, 0, s_final.size(), m, Ds_intercept, r2);

  out_file << "#avg = " << avg << std::endl;
  out_file << "#D (ACF) = " << var*var/I*0.1 << " cm2/s " << std::endl;
  out_file << "#Ds (VACF) = " << Ds_intercept*0.1 << " cm2/s " << std::endl;
  
  out_file.close();

  delete[] acf;
  delete[] vacf;
  //  series.earse();
}
