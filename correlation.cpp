#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <numeric>

#define INTCUTOFF 0.05

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
using namespace std;

/* Allen, M.; Tildesley, D. Computer Simulation of Liquids; Oxford Science Publications, Clarendon Press: Oxford, 1989. */

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
  cout<< "velseries " << nSamples << endl;
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
      //      if(acf[i]<threshhold)
      //	break;
      
      I+=0.5*(acf[i]+acf[i+1])*timestep;
    }
  return(I);
}

// read time series from filename fname
// file is formed like NAMD colvar traj file
// store the number of points in numSamples
// each line should store one time step

vector<double> readSeriesNAMD(char *fname, int &numSamples, int field)
{
  ifstream datafile(fname);
  vector<double> series;
  double *timeSeries;
  int i=0;
  string line;
  istringstream iss;
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
	      string str2=line.substr(begin,23);
	      series.push_back(atof(str2.c_str()));
	      ++numSamples;
	    }
	  catch (exception& e)
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

vector<double> readSeriesRaw(char *fname, int &numSamples)
{
  ifstream datafile(fname);
  vector<double> series;
  double *timeSeries;
  int i=0;
  string line;
  istringstream iss;
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
	  catch (exception& e)
	    {
	      break;
	    }
	}
    }

  return(series);
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


// calculates intercept by linear interpolation

double limit(const vector<double>& x, const vector<double>& y)
{
    int n=x.size();
    double s_x  = accumulate(x.begin(), x.end(), 0.0);
    double s_y  = accumulate(y.begin(), y.end(), 0.0);
    double s_xx = inner_product(x.begin(), x.end(), x.begin(), 0.0);
    double s_xy = inner_product(x.begin(), x.end(), y.begin(), 0.0);
    double m=(n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    double b=(s_y-m*s_x)/n;
    
    cout << "#m " << m << endl;
    cout << "#b " << b << endl;
    return(b);
}

int main(int argc, char *argv[])
{
  int nCorr=1000;
  vector<double> series, seriesVel;
  double *velSeries;
  double *acf, *vacf, *timeSeries;
  double var, varVel, I;
  char *fname;
  double timestep=1.0;
  double var_m2;
  int field=1;
  int numSamples;
  double varAnalytical;
  double k=10.0*4.184*1000;
  double varVelAnalytical;
  
  varAnalytical=8.314*298.15/k;
  varVelAnalytical=8.314*298.15/(18.01/1000.0)*1E-10;
  
  cout << "varAnalytical " << varAnalytical << endl;
  cout << "varVelAnalytical " << varVelAnalytical << endl;
  
  if(argc<1)
    return(1);

  fname=argv[1];

  if(argc>1)
    field=atoi(argv[2]);

  series=readSeriesNAMD(fname, numSamples, field);
  timeSeries=&series[0];

  subtract_average(timeSeries, numSamples);
  velSeries=velocity_series(timeSeries, numSamples, timestep);
  
  numSamples=numSamples-2;

  acf=calcCorrelation(timeSeries, numSamples, nCorr);
  vacf=calcCorrelation(velSeries, numSamples, nCorr);

  var=variance(timeSeries, numSamples);
  varVel=variance(velSeries, numSamples);
  
  ofstream myfile;
  myfile.open(argv[3]);

  for(int i=0;i<nCorr;++i)
    myfile << i << " " << acf[i] << " " << vacf[i] << endl;
  myfile.close();

  I=integrateCorr(acf, nCorr, timestep);
  
  //  double s=0.003;
  
  vector<double> s_values;
  vector<double> intDs;

  cout << var << endl;
  cout << "#varVel " << varVel << endl;
  cout << "#varVelAnalytical " << varVelAnalytical << endl;

  double s=0.01;
  while(s<=0.04)
    {
      double laplaceVACF=laplace(vacf, timestep, s, nCorr);
      double Ds=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);

      cout << s << " " <<  Ds << endl;
      s_values.push_back(s);
      intDs.push_back(Ds);
      s=s+0.0001;
    }
  
  double limDs=limit(s_values, intDs);

  cout << "#I = " << I << endl;
  cout << "#var = " << var << endl;
  cout << "#D = " << var*var/I << " A2/fs " << endl;
  cout << "#D = " << var*var/I*0.1 << " cm2/s " << endl;
  cout << "#Ds= " << limDs << " A2/s " << endl;

  delete[] acf;
  delete[] vacf;
  //  series.earse();
}
