#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include <algorithm>

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

/*
double *histogram(double *y, int nSamples)
{
  double *hist=new double[100];
  int bin;
  
  for(int i=0;i<100;++i)
    hist[i]=0.0;
  
  for(int i=0;i<nSamples;++i)
    {
      bin=50.0*y[i]+50;
      if(bin>0 && bin<100)
	hist[bin]++;
    }
  return(hist);
}
*/

int main(int argc, char *argv[])
{
  int nCorr=5000;
  vector<double> series;
  double *acf, *timeSeries;
  double var, I;
  char *fname;
  double timestep=2.0;
  int field=1;
  int numSamples;

  if(argc<1)
    return(1);

  fname=argv[1];

  if(argc>1)
    field=atoi(argv[2]);
  //  int numSamples=countLines(fname)-1;

  series=readSeries(fname, numSamples, field);
  timeSeries=&series[0];

  numSamples=numSamples-1;

  subtract_average(timeSeries, numSamples);
  acf=calcCorrelation(timeSeries, numSamples, nCorr);
  var=variance(timeSeries, numSamples);

  ofstream myfile;
  myfile.open(argv[3]);
  for(int i=0;i<nCorr;++i)
    myfile << i << " " << acf[i] << endl;
  myfile.close();
  
  I=integrateCorr(acf, nCorr, timestep);

  cout << "#I = " << I << endl;
  cout << "#var = " << var << endl;
  cout << "#D = " << var*var/I << " A2/fs " << endl;
  cout << "#D = " << var*var/I*0.1 << " cm2/s " << endl;

  for(int i=0;i<nCorr;++i)
    cout << i << " " << acf[i] << endl;
  delete[] acf;
  //  series.earse();
}
