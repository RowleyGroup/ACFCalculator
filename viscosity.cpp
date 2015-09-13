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

double **readSeriesNAMD(char *fname, int &numSamples, int field)
{
  ifstream datafile(fname);
  vector<vector<double> > series(9, vector<double>(10));
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
      if(line.compare(0, 9, "PRESSURE:")==0)
	{
	  istringstream iss(line);
	  // skip first two files

	  string sub;
	  iss >> sub;
	  iss >> sub;

	  for(i=0;i<9;++i)
	    {
	      iss >> sub;
	      string::size_type sz;
	      series[i].push_back(stod(sub,&sz));
	    }
	  ++numSamples;
	}
    }

  double **seriesArray=new double*[9];

  for(int i=0;i<9;++i)
    seriesArray[i]=new double[numSamples];

  for(int i=0;i<9;++i)
    for(int j=0;j<numSamples;++j)
      seriesArray[i][j]=series[i][j];

  ///  cout << numSamples << endl;

  return(seriesArray);
}

int main(int argc, char *argv[])
{
  int nCorr=1000;
  vector<double> series, seriesVel;
  double **seriesArray;
  double *velSeries;
  double *acf;
  double var, varVel, I;
  char *fname;
  double timestep=2.0;
  double var_m2;
  int field=1;
  int numSamples;
  vector<double> viscosityAll;

  double boxLength=31.0399376148;
  double kB=1.38064852E-23;
  double T=298.15;

  fname=argv[1];

  if(argc<1)
    return(1);

  seriesArray=readSeriesNAMD(fname, numSamples, field);

  for(int i=0;i<9;++i)
    {
      // use only off-diagonal pressure tensor series
      if(i%4==0)
	continue;

      acf=calcCorrelation(seriesArray[i], numSamples, nCorr);
      I=integrateCorr(acf, nCorr, timestep);
      cout << "#" << I << endl;
      double eta=boxLength*boxLength*boxLength/(kB*T)*I*1E-30*1E-15*1E10;
      cout << "#eta " << i << " " << eta << endl;
      viscosityAll.push_back(eta);
      for(int j=0;j<nCorr;++j)
	{
	  cout << j << " " << acf[j] << endl;
	}
      cout << endl;
      delete[] acf;
    }

  double sum = accumulate(viscosityAll.begin(), viscosityAll.end(), 0.0);
  double mean = sum / viscosityAll.size();

  double sq_sum = inner_product(viscosityAll.begin(), viscosityAll.end(), viscosityAll.begin(), 0.0);
  double stdev = sqrt(sq_sum / viscosityAll.size() - mean * mean);
  cout << "#eta = " << mean  << " " << stdev << endl;

}
