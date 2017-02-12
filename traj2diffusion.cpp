#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <numeric>

#include <boost/program_options.hpp>


#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define PI 3.1415926
#define RECSCALE32BIT 1
#define RECSCALE64BIT 2
#define RECSCALEMAX   2

#define BOLTZ 1.38065E0
//#define DEFAULT_T 373.0
#define DEFAULT_T 298.15
#define DEFAULT_DELTA 0.001
#define DEFAULT_TIMESTEP 1.0
#define DEFAULT_STRIDE 1
#define DEFAULT_MAX_S 1000

namespace po = boost::program_options;

double calcmsd(double *crd1, double *crd2)
{
  double dx, dy, dz;
  
  dx=crd2[0]-crd1[0];
  dy=crd2[1]-crd1[1];
  dz=crd2[2]-crd1[2];

  return( dx*dx + dy*dy + dz*dz);
}

// read time series from filename fname
// file is formed like NAMD colvar traj file
// store the number of points in numSamples
// each line should store one time step

double **readSeriesNAMD(char *fname, int &numSamples)
{
  std::ifstream datafile(fname);
  std::vector<double> seriesX, seriesY, seriesZ;
  std::string line;
  std::istringstream iss;
  double **crd;

  if(!datafile)
    {
      std::cerr << "Error: Cannot read " << fname << std::endl;
    }
  
  numSamples=0;

  while(getline(datafile,line))
    {
      if(line.at(0)!='#')
        {
          try
            {
	      std::string str2=line.substr(15,23);
              seriesX.push_back(atof(str2.c_str()));

	      std::string str3=line.substr(37,23);
              seriesY.push_back(atof(str3.c_str()));

	      std::string str4=line.substr(61,23);
              seriesZ.push_back(atof(str4.c_str()));
	      
              ++numSamples;
            }
          catch (std::exception& e)
            {
              break;
            }
        }
    }
  
  try
    {
      crd=new double*[numSamples];
      
      for(int i=0;i<numSamples;++i)
	{
	  crd[i]=new double[3];
	  crd[i][0]=seriesX[i];
	  crd[i][1]=seriesY[i];
	  crd[i][2]=seriesZ[i];
	}
    }
  catch(std::bad_alloc&)
    {
      std::cerr << "Could not allocate memory." << std::endl;
      exit(1);
    }

    seriesX.clear();
    seriesY.clear();
    seriesZ.clear();
    
    return(crd);
}

bool isNumeric(char in)
{
  if(isdigit(in))
    return(true);
  
  switch(in)
    {
    case('-'):
      return(true);
    case('+'):
      return(true);
    case('.'):
      return(true);
    case('e'):
      return(true);
    case('E'):
      return(true);
    }
  return(false);
}
 

double **readSeriesGeneral(char *fname, int &numSamples)
{
  std::ifstream datafile(fname, std::ifstream::in);
  std::vector<double> series;
  std::string line;
  std::vector<double> seriesX, seriesY, seriesZ;
  double **crd;
  
  numSamples=0;
  while(getline(datafile,line))
    {
      try
	{
	  //	  std::cout << line[0] << std::endl;
	  if(isNumeric(line[0]) || line[0]==' ')
	    {
	      int start=0;
	      int end=0;
	      
	      // skip first field (time step column)
	      while(start<line.length() && isNumeric(line[start])==false)
		++start;
	      
	      end=start;
	      
	      while(end<line.length() && isNumeric(line[end]))
		++end;
	      
	      start=end;
	      
	      // get x value
	      
	      // skip non numeric values leading string
	      
	      while(start<line.length() && isNumeric(line[start])==false)
		++start;
	      
	      // find end of string
	      end=start;
	      
	      while(end<line.length() && isNumeric(line[end]))
		++end;
	      
	      std::string str2=line.substr(start,end-start);
	      seriesX.push_back(atof(str2.c_str()));
	      
	      // get  y value
	      
	      // skip non numeric values leading string
	      
	      while(start<line.length() && isNumeric(line[start])==false)
		++start;

	      // find end of string
	      end=start;
	      
	      while(end<line.length() && isNumeric(line[end]))
		++end;
	      	      
	      std::string str3=line.substr(start,end-start);

	      seriesY.push_back(atof(str2.c_str()));

	      // get z value
	      
	      // skip non numeric values leading string
	      start=end;
	      
	      while(start<line.length() && isNumeric(line[start])==false)
		++start;
	      
	      // find end of string
	      end=start;
	      
	      while(end<line.length() && isNumeric(line[end]))
		++end;
	      
	      std::string str4=line.substr(start, end-start);

	      seriesZ.push_back(atof(str2.c_str()));
	      
	      ++numSamples;
	    }
	}
      catch (std::exception& e)
	{
	  break;	
	}
    }
  
  try
    {
      crd=new double*[numSamples];
      
      for(int i=0;i<numSamples;++i)
	{
	  crd[i]=new double[3];
	  crd[i][0]=seriesX[i];
	  crd[i][1]=seriesY[i];
	  crd[i][2]=seriesZ[i];
	  //	  std::cout << i << " " << crd[i][0] << " " << crd[i][1] << " " << crd[i][2] << std::endl;
	}
    }
  catch(std::bad_alloc&)
    {
      std::cerr << "Could not allocate memory." << std::endl;
      exit(1);
    }
  
  seriesX.clear();
  seriesY.clear();
  seriesZ.clear();
  
  return(crd);
}


double dist(double **coords[3], int a, int b)
{
    double rxab, ryab, rzab, r;

    rxab=coords[b][0]-coords[a][0];
    ryab=coords[b][1]-coords[a][1];
    rzab=coords[b][2]-coords[a][2];

    r=sqrt(rxab*rxab+ryab*ryab+rzab*rzab);
    
    return(r);
}

void image(double **crd, int nframes, double L)
{
  double **pbc_images;
  
  try
    {
      pbc_images=new double*[nframes];
      for(int i=0;i<nframes;++i)
	{
	  pbc_images[i]=new double[3];
	  pbc_images[i][0]=0.0;
	  pbc_images[i][1]=0.0;
	  pbc_images[i][2]=0.0;
	}
    }
  catch(std::bad_alloc&)
    {
      std::cerr << "Error: Could not allocate memory." << std::endl;
      exit(1);
    }
  
  #pragma omp parallel for
  for(int i=1;i<nframes;++i)
    {
      for(int j=0;j<3;++j)
	{
	  double delta=crd[i][j]-crd[i-1][j];

	  if(delta > L/2)
	    {
	      pbc_images[i][j]-=L;
	    }
	  
	  if(delta < -L/2)
	    {
	      pbc_images[i][j]+=L;
	    }
	}
    }

  // accumlate pbc_images
  
  for(int i=1;i<nframes;++i)
    {
      #pragma omp parallel for
      for(int j=0;j<3;++j)
	{
	  pbc_images[i][j]+=pbc_images[i-1][j];
	}
    }

  #pragma omp parallel for
  for(int i=1;i<nframes;++i)
    {
      for(int j=0;j<3;++j)
	{
	  crd[i][j]+=pbc_images[i][j];
	}
    }
  
  for(int i=0;i<nframes;++i)
    {
      delete[] pbc_images[i];
    }
  delete[] pbc_images;
}

double slope(double *msd, int n)
{
  int i;
  double sum=0.0;

  for(i=1;i<n;++i)
    {
      sum+=(msd[i]/i);
    }

  printf("#sum = %lf\n", sum); 
  return(sum/n);
}


int main(int argc, char **argv)
{
  enum InputType {namd, general};
  InputType type;
  
  int i, j, k, n;
  
  float *atomcoords=NULL;
  
  int natoms, natom_res, natom_firstres=0;
  
  double ave_volume;
  double total_volume=0.0;
  double total_mass=0.0;
  double density;
  
  int count;
  double T=DEFAULT_T;
  
  FILE *output_fh;
  
  int max_s=1000;
  
  double dt=DEFAULT_DELTA;
  
  double **crd;
  double *r_arr, *x_arr, *y_arr, *z_arr;
  
  int mflag=0;
  int stride=100;
  int nframes;
  bool imaging=false;
  double L=0.0;
  int c;
  
  int *msd_counter;
  double *t, *msd;
  double D;
  char *traj_fname=NULL;
  double timestep=1.0; // stores timestep in fs
    
  po::options_description desc("Allowed options", 1024, 512);
  
  desc.add_options()
    ("help,h", "produce help message")
    ("traj,t", po::value<std::string>()->required(), "file name containing the time series (e.g., the colvar.traj file")
    ("timestep,d", po::value<double>(&timestep)->default_value(DEFAULT_TIMESTEP), "time step of the simulation (fs)")
    ("type,f", po::value<std::string>(), "type of time series (namd, general)")
    ("stride,s", po::value<int>(&stride)->default_value(DEFAULT_STRIDE), "number of steps between the stride")
    ("cell,c", po::value<double>(&L), "length of unit cell")
    ("max,m", po::value<int>(&max_s)->default_value(DEFAULT_MAX_S), "maximum number of steps to calculate MSD over")
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
      
      if (vm.count("traj"))
	{
	  const std::string traj_fname_str=vm["traj"].as<std::string>();
	  traj_fname = new char[traj_fname_str.length()+1];
	  std::strcpy(traj_fname, traj_fname_str.c_str());
	}
      else
	{
	  std::cout << "Error: The number of the trajectory file must be provided." << std::endl;
	  exit(1);
	}
      
      // stride is the number of time steps between entries in the colvars file
      if (vm.count("stride"))
	{
	  stride=vm["stride"].as<int>(); 
	}
      else
	{
	  std::cout << "Stride of colvars file must be provided." << std::endl;
	  exit(1);
	}
      
      if (vm.count("timestep"))
	{
	  dt=vm["timestep"].as<double>();
	}
      else
	{
	  std::cout << "Time step of simulation must be provided" << std::endl;
	  exit(1);
	}
      
      if (vm.count("cell"))
	{
	  L=vm["cell"].as<double>();
	  imaging=true;
	}
      else
	{
	  std::cout << "#No imaging is performed. This assumes that the simulation code does NOT wrap positions from the periodic boundary conditions" << std::endl;
	}
      
      if(vm.count("type"))
	{
	  const std::string type_str=vm["type"].as<std::string>();
	  if(type_str.compare("namd")==0)
	    {
	      type=namd;
	    }
	}
      else
	{
	  std::cout << "#Type of file not provided, defaulting to general file reader" << std::endl;
	  type=general;
	}
    }
  catch ( const std::exception& e )
    {
      std::cerr << e.what() << std::endl;
      return 1;
    }
  
  if(traj_fname==NULL)
    {
      std::cerr << "Error: Time series file name must be provided" << std::endl;
      exit(0);
    }
  
  double spaceSeries=dt*stride;

  if(type==namd)
    {
      crd=readSeriesNAMD(traj_fname, nframes);
    }
  else
    {
      crd=readSeriesGeneral(traj_fname, nframes);
    }
      
  if(imaging)
    image(crd, nframes, L);
  
  std::cout << "#" << nframes << std::endl;
  
  if(nframes<0)
    {
      std::cerr << "Error: Skipping more frames than there are in the trajectory." << std::endl;
      exit(1);
    }
  
  count=0;
  
  try
    {
      msd=new double[max_s];
      msd_counter=new int[max_s];
    }
  catch(std::bad_alloc&)
    {
      std::cerr << "Error: Could not allocate memory." << std::endl;
      exit(1);	
    }
  
  bzero(msd, max_s*sizeof(double));
  bzero(msd_counter, max_s*sizeof(int));
  
  for(i=0;i<nframes;i+=stride)
    {
      for(j=0;j<max_s;j++)
	{
	  if( (i+j)>=nframes)
	    break;
	  
	  msd[j]+=calcmsd(crd[i], crd[i+j]);
	  msd_counter[j]=msd_counter[j]+1;
	}
    }
  
  for(i=0;i<max_s;++i)
    {
      msd[i]/=msd_counter[i];
      std::cout << i*dt*stride << " " << msd[i] << std::endl;
    }
  
  D=slope(msd, max_s)/spaceSeries/6.0;
  
  printf("#D %le (A2/fs)\n", D);
  printf("#D %le (cm2/s)\n", D*1E-1);
  printf("#D %le (m2/s)\n", D*1E-5);
  
  return(0);
}
