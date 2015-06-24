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
  //  cout << "variance << " << v2 << endl;
  v2/=nSamples;
  return(v2);
}

void subtract_average(double *y, int nSamples)
{
  double avg=0.0;

  for(int i=0;i<nSamples;++i)
    avg+=y[i];
 
  avg/=nSamples;

  for(int i=0;i<nSamples;++i)
    y[i]-=avg;
}
// calculate the velocity time series by numerically differentiating the position time series

double *velocity_series(double *y, int nSamples, double timestep)
{
  double *vel=new double[nSamples-1];
  double avg=0.0;

  // calculate velocity by finite difference approach and convert to m/s
  for(int i=0;i<nSamples-1;++i)
    vel[i]=(y[i+1]-y[i])/timestep*1E-10/1E-15;

  for(int i=0;i<nSamples-1;++i)
    avg+=vel[i];

  avg/=nSamples-1;

  for(int i=0;i<nSamples-1;++i)
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

double chi2(double *y, int nSamples, double var)
{
  long double chi2=0.0;
  for(int i=0;i<nSamples;++i)
    {
      //      chi2=(y[i]-
    }
}

vector<double> readSeries(char *fname, int &numSamples, int field)
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

double laplace(double *series, double timestep, double s, int length)
{
  double F=0.0;
  
  for(int i=0;i<length;++i)
    {
      F+=exp(-s*i*timestep)*series[i]*timestep;
    }
  return(F);
}

int main(int argc, char *argv[])
{
  int nCorr=5000;
  vector<double> series;
  double *vel;
  double *acf, *vacf, *timeSeries;
  double var, varVel, I;
  char *fname;
  double timestep=2.0;
  double timestep_s=timestep*1E-15;
  int field=1;
  int numSamples;

  if(argc<1)
    return(1);

  fname=argv[1];

  if(argc>1)
    field=atoi(argv[2]);

  series=readSeries(fname, numSamples, field);
  timeSeries=&series[0];

  numSamples=numSamples-1;

  subtract_average(timeSeries, numSamples);
  vel=velocity_series(timeSeries, numSamples, timestep);
  acf=calcCorrelation(timeSeries, numSamples, nCorr);
  vacf=calcCorrelation(vel, numSamples-1, nCorr);
  var=variance(timeSeries, numSamples);
  varVel=variance(vel, numSamples-1);
  
  ofstream myfile;
  myfile.open(argv[3]);
  
  for(int i=0;i<nCorr;++i)
    myfile << i << " " << acf[i] << " " << vacf[i] << endl;
  myfile.close();
  
  I=integrateCorr(acf, nCorr, timestep);

  double s=0.00001;
  while(s<0.01)
    {
      double laplaceVACF=laplace(vacf, timestep_s, s, nCorr);
      double Ds=-(laplaceVACF*var*varVel)/(laplaceVACF*(s*var+varVel/s)-var*varVel);
      cout << s << " " <<  laplaceVACF << " " << Ds << endl;
      s=s+0.00001;
    }
  cout << "#I = " << I << endl;
  cout << "#var = " << var << endl;
  cout << "#D = " << var*var/I << " A2/fs " << endl;
  cout << "#D = " << var*var/I*0.1 << " cm2/s " << endl;

  delete[] acf;
  delete[] vacf;
  delete[] timeseries;
  //  series.earse();
}
