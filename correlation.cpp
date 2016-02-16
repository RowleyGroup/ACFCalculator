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

#define INTCUTOFF 0.05
#define DEFAULT_MAXCORR 1000
#define DEFAULT_TIMESTEP 1
#define SEG_MIN 0.2
#define R2_THRESHOLD 0.9999
#define NSEGBINS 10

#define LINEAR_BINS 10 // define range into bins

/* compile: g++ correlation.cpp -o correlation */

/* syntax: ./correlation time-series.traj 

/* the time series is assumed to be in the format saved by the NAMD colvars module */

/* calculate diffusion coefficient of a harmonically-restrained degree of freedom by */
/* calculating the autocorrelation function a time series of a molecular dynamics */
/* simulation then integrating it */

/* Hummer, G. Position-dependent diffusion coefficients and free energies from Bayesian
  analysis of equilibrium and replica molecular dynamics simulations. New J. Phys. 2005,
  7, 34. */

/* internal units:
distance = angstrom
time = femptoseconds
velocity = Angstrom / femptoseconds
*/

namespace po = boost::program_options;

// http://vilipetek.com/2013/10/07/polynomial-fitting-in-c-using-boost/

std::vector<double> polyfit( const std::vector<double>& oX,
                        const std::vector<double>& oY, int nDegree )
{
  using namespace boost::numeric::ublas;

  if ( oX.size() != oY.size() )
    throw std::invalid_argument( "X and Y vector sizes do not match" );

  // more intuative this way
  nDegree++;

  size_t nCount =  oX.size();
  matrix<double> oXMatrix( nCount, nDegree );
  matrix<double> oYMatrix( nCount, 1 );

  // copy y matrix
  for ( size_t i = 0; i < nCount; i++ )
    {
      oYMatrix(i, 0) = oY[i];
    }

  // create the X matrix
  for ( size_t nRow = 0; nRow < nCount; nRow++ )
    {
      double nVal = 1.0f;
      for ( int nCol = 0; nCol < nDegree; nCol++ )
        {
          oXMatrix(nRow, nCol) = nVal;
          nVal *= oX[nRow];
        }
    }

  // transpose X matrix
  matrix<double> oXtMatrix( trans(oXMatrix) );
  // multiply transposed X matrix with X matrix
  matrix<double> oXtXMatrix( prec_prod(oXtMatrix, oXMatrix) );
  // multiply transposed X matrix with Y matrix
  matrix<double> oXtYMatrix( prec_prod(oXtMatrix, oYMatrix) );

  // lu decomposition
  permutation_matrix<int> pert(oXtXMatrix.size1());
  const std::size_t singular = lu_factorize(oXtXMatrix, pert);
  // must be singular
  BOOST_ASSERT( singular == 0 );

  // backsubstitution
  lu_substitute(oXtXMatrix, pert, oXtYMatrix);

  // copy the result to coeff
  return std::vector<double>( oXtYMatrix.data().begin(), oXtYMatrix.data().end() );
}

double *vacf;
int nCorr;
double var;
double varVel;
double timestep;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double laplace(double *series, double timestep, double s, int length)
{
  double F=0.0;

  for(int i=0;i<length;++i)
    {
      F+=exp(-s*i*timestep)*series[i]*timestep;
    }
  return(F);
}


double denom(double s)
{
  double laplaceVACF=laplace(vacf, timestep, s, nCorr);
  double f=laplaceVACF*(s*var+varVel/s)-var*varVel;
  return(f);
}

double Ds(double s)
{
  double laplaceVACF=laplace(vacf, timestep, s, nCorr);
  double f=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);
  return(f);
}

double bisect(double lower, double upper)
{
  int it=0;
  double fa=denom(lower);
  double fb=denom(upper);
  double fmid=0.0;
  double mid;

  while(it<1000)
    {
      mid=(lower+upper)/2.0;
      fmid=denom(mid);

      if(fabs(fmid)<1E-10)
        {
          return(mid);
        }

      if(sgn(fmid)==sgn(fa))
        {
          lower=mid;
          fa=fmid;
        }
      else
        {
          upper=mid;
          fb=fmid;
        }
      ++it;
    }
  return(mid);
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
  int min;
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

void subtract_average(double *y, int nSamples)
{
  double avg=0.0;

  for(int i=0;i<nSamples;++i)
    {
      avg+=y[i];
    }
  avg/=nSamples;

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
  double threshhold;

  threshhold=acf[0]*acf[0]*INTCUTOFF;

  for(int i=0;i<nCorr-1;++i)
    {
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
  double *timeSeries;
  int i=0;
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

std::vector<double> readSeriesRaw(char *fname, int &numSamples)
{
  std::ifstream datafile(fname, std::ifstream::in);
  std::vector<double> series;
  double *timeSeries;
  int i=0;
  std::string line;
  std::istringstream iss;
  int begin;

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
//

std::vector<double> readSeriesGROMACS(char *fname, int &numSamples)
{
  std::ifstream datafile(fname, std::ifstream::in);
  std::vector<double> series;
  double *timeSeries;
  int i=0;
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
	      while(getline(iss,item,' '))
		{
		  series.push_back(atof(item.c_str()));
		}
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
    int n=x.size();
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

// iteratively finds linear segment of D(s)
// trims front or back percentage
// stops when r2 is > 0.99 or length is now less than SEG_MIN fraction

void find_linear_range(std::vector<double>& s, std::vector<double>& Ds, int &lower, int &upper)
{
  double m, b, r2=0.0;
  int min_length=(upper-lower)*SEG_MIN;
  
  do
    {
      leastsquares_subset(s, Ds, lower, upper, m, b, r2);
      
      int delete_length=(upper-lower)*0.05;
      
      double mM=0.0, bM=0.0, r2M=0.0;
      leastsquares_subset(s, Ds, lower+delete_length, upper-delete_length, mM, bM, r2M);

      double mF=0.0, bF=0.0, r2F=0;
      leastsquares_subset(s, Ds, lower, lower+delete_length, mF, bF, r2F);

      double mB=0.0, bB=0.0, r2B=0.0;
      leastsquares_subset(s, Ds, upper-delete_length, upper, mB, bB, r2B);

      // test which section to delete.  Delete the least linear segment
      if(fabs(r2B) < fabs(r2F))
	{
	  // the front is closer to the middle than the back. delete the back segment
	  upper-=delete_length;
	}
      else
	{
	  // the back is closer to the middle than the front. delete the front segment
	  lower+=delete_length;
	}
    }
  while(r2<R2_THRESHOLD && (upper-lower)>min_length);
}

double exp_fit_linearregression(const std::vector<double>& x, const std::vector<double>& y)
{
  std::vector<double> lny;
  int n=x.size();

  for(int i=0;i<x.size();++i)
    {
      lny.push_back(log(y[i]));
    }
  double intercept=leastsquares(x, lny);

  return(exp(intercept));
}

// searches backwards from second singularity until line tangent
// to the curve would give a non-negative intercept

int find_negative_bound(const std::vector<double>& x, const std::vector<double>& y)
{
  for(int i=x.size()-2;i>0;--i)
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

int main(int argc, char *argv[])
{
  enum InputType {gromacs, namd, raw};
  InputType  type;
  
  bool write_acf;
  std::vector<double> series, seriesVel;
  double *velSeries;
  double *acf, *timeSeries;
  double I;
  char *fname, *acf_fname, *output_fname;
  double var_m2;
  int field=1;
  int numSamples;
  double varAnalytical;
  double k=10.0*4.184*1000;
  double varVelAnalytical;
  double varVel;
  double cutoff;

  timestep=1.0;
  varAnalytical=8.314*298.15/k;
  varVelAnalytical=8.314*298.15/(18.01/1000.0)*1E-10;

  po::options_description desc("Allowed options");

  desc.add_options()
    ("help,h", "produce help message")
    ("input,i", po::value<std::string>()->required(), "file name of time series")
    ("type,t", po::value<std::string>(), "type of time series (namd, charmm, gromacs, amber, txt)")
    ("cutoff,c", po::value<double>()->default_value(INTCUTOFF), "cutoff to integrate ACF")
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
      po::notify(vm);

      if (vm.count("help"))
        {
          return 1;
        }

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
	  if(type_str.compare("gromacs"))
	    type=gromacs;
	  else if(type_str.compare("raw"))
	    type=raw;
	  else 
	    type=namd;
	}
     else
       {
	 std::cout<<"Type of file not provided"<<std::endl;
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

      if (vm.count("timestep"))
        {
          timestep=vm["timestep"].as<double>();
          std::cout << "#Time step of " << timestep << " fs will be used." << std::endl;
        }
      else
        {
          timestep=DEFAULT_TIMESTEP;
          std::cout << "#Default timestep of " << std::endl;
        }

      if(vm.count("cutoff"))
      {
          cutoff=vm["cutoff"].as<double>();
          std::cout<<"#Cutoff of "<<cutoff<<" will be used."<<std::endl;
      }
      else
      {
          cutoff=INTCUTOFF;
          std::cout<<"Default cutoff of "<<cutoff<<" will be used."<<std::endl;
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
  else if(type=gromacs)
    {
      series=readSeriesGROMACS(fname, numSamples);
    }
  else
    {
      series=readSeriesRaw(fname, numSamples);
    }
  timeSeries=&series[0];

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

  I=integrateCorr(acf, nCorr, timestep);

  std::ofstream out_file(output_fname,  std::ofstream::out);
  std::cout << "#output " << output_fname << std::endl;

  out_file << "#var " << var << std::endl;
  out_file << "#varVel " << varVel << std::endl;
  out_file << "#varVelAnalytical " << varVelAnalytical << std::endl;
  out_file << "#nCorr " << nCorr << std::endl;

  double s=0.0001;

  double s_min_num=1000;
  double s_min_loc=-1.0;
  out_file << "#" << std::setw(15) << std::left << "s"<< std::setw(15) << std::left <<"Ds" << std::setw(15) << std::left << "laplaceVACF" << std::setw(15) << std::left << "Ds numerator" << std::setw(15) << std::left << "Ds denominator" << std::endl;

  while(s<=1.0)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double denom=(laplaceVACF*(s*var+varVel/s)-var*varVel);
      double Ds=-(laplaceVACF*var*varVel)/denom;

      if(denom<s_min_num)
        {
          s_min_num=denom;
          s_min_loc=s;
        }

      out_file << std::setw(15) << std::left << s<< std::setw(15) << std::left <<  Ds << std::setw(15) << std::left << laplaceVACF << std::setw(15)<< std::left << -(laplaceVACF*var*varVel) << std::setw(15) << std::left << (laplaceVACF*(s*var+varVel/s)-var*varVel) << std::endl;
      s=s+0.0001;
    }

  // find singularity numerically
  s=0.0001;

  double root=1.0;
  double s2=0.0;
  double root_one=s_min_loc;
  double root_two=0.0;
  int both = 0;

  while(s<=0.1)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double denom=(laplaceVACF*(s*var+varVel/s)-var*varVel);
      if(denom<0)
        {
          root=denom;
          root_one=s;
	  break;
        }
      s=s+0.0001;
    }

  // find second singularity
  // occurs when denominator becomes positive again
  while(s<=1.0)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double denom=(laplaceVACF*(s*var+varVel/s)-var*varVel);

      if(denom>0)
        {
	  root_two=s;
	  break;
        }

      s=s+0.0001;
    }
  
  double Ds=1.0;
  double Ds_prev=1.0;
  
  out_file << "#1st singularity " << root_one << std::endl;
  out_file << "#2nd singularity "<< root_two << std::endl;

  // find minimum in Ds between first and second signularities
  s=root_one;
  while(s<=root_two)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double denom=(laplaceVACF*(s*var+varVel/s)-var*varVel);
      Ds_prev=Ds;
      Ds=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);
      if(Ds>Ds_prev)
        {
          break;
        }
      s=s+0.00001;
    }

  double Ds_min=s;

  out_file << "#Ds_min " << Ds_min << std::endl;
  
  // Step2: calculate s over range
  //  s = (root_two - Ds_min)/3.3;
  //  double s_max = (root_two - Ds_min)/2.4;
  //  double s_delta=1.0*(s_max - s)/100; // calculate D(s) at 100 point

  double s_range=root_two-Ds_min;
  double s_min=Ds_min;
  double s_max=root_two;
  
  //  double s_max=Ds_min+s_range*0.4;
  //  double s_max=(root_two - Ds_min)/2.4;
  s=s_min;
  double s_delta=s_range/1000;
  
  std::vector<double> s_values;
  std::vector<double> intDs;

  out_file << std::setw(15) << std::left << "#s"<< std::setw(15) << std::left <<"Ds" << std::setw(15) << std::left << "laplaceVACF" << std::setw(15) << std::left << "Ds numerator" << std::setw(15) << std::left << "Ds denominator" << std::endl;
  while(s<=s_max)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double Ds=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);
      double logDs=log(Ds);
      out_file << std::setw(15) << std::left << s<< std::setw(15) << std::left <<  Ds << std::setw(15) << std::left << laplaceVACF << std::setw(15)<< std::left << -(laplaceVACF*var*varVel) << std::setw(15) << std::left << (laplaceVACF*(s*var+varVel/s)-var*varVel) << std::endl;
      s_values.push_back(s);
      intDs.push_back(Ds);
      s=s+s_delta;
    }
  
  int upper=find_negative_bound(s_values, intDs);
  int lower=0;
  
  find_linear_range(s_values, intDs, lower, upper);
  
  // Step 3
  double m, b, r2;
  
  leastsquares_subset(s_values, intDs, lower, upper, m, b, r2);

  out_file << "#I = " << I << std::endl;
  out_file << "#var = " << var << std::endl;
  out_file << "#D = " << var*var/I << " A2/fs " << std::endl;
  out_file << "#D = " << var*var/I*0.1 << " cm2/s " << std::endl;

  out_file << "#Ds (sub)= " << b*0.1 << " cm2/s " << std::endl;

  out_file.close();

  delete[] acf;
  delete[] vacf;
  //  series.earse();
}
