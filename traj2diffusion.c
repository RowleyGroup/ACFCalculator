#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

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

using namespace std;

int pattern_match(const char *str, const char *pattern) {
  enum State {
    Exact,      // exact match
    Any,        // ?
    AnyRepeat    // *
  };

  const char *s = str;
  const char *p = pattern;
  const char *q = 0;
  int state = 0;

  int match = 1;

  while (match && *p) 
    {
      if (*p == '*')
	{
	  state = AnyRepeat;
	  q = p+1;
	} 
      else if (*p == '?') state = Any;
      else state = Exact;

      if (*s == 0) break;

      switch (state) 
	{
	case Exact:
	  match = *s == *p;
	  s++;
	  p++;
	  break;
	  
	case Any:
	  match = 1;
	  s++;
	  p++;
	  break;

	case AnyRepeat:
	  match = 1;
	  s++;
	  
	  if (*s == *q) p++;
	  break;
	}
    }
  
  if (state == AnyRepeat) return (*s == *q);
  else if (state == Any) return (*s == *p);
  else return match && (*s == *p);
}

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
  ifstream datafile(fname);
  vector<double> seriesX, seriesY, seriesZ;
  double *timeSeries;
  int i=0;
  string line;
  istringstream iss;
  int begin;
  double **crd;

  if(!datafile)
    {
      cerr << "Error: Cannot read " << fname << endl;
    }
  
  numSamples=0;

  while(getline(datafile,line))
    {
      if(line.at(0)!='#')
        {
          try
            {
              string str2=line.substr(15,23);
              seriesX.push_back(atof(str2.c_str()));

              string str3=line.substr(37,23);
              seriesY.push_back(atof(str3.c_str()));

              string str4=line.substr(61,23);
              seriesZ.push_back(atof(str4.c_str()));
	      
              ++numSamples;
	      
            }
          catch (exception& e)
            {
              break;
            }
        }
    }
 
    if( (crd=(double **) calloc(numSamples, sizeof(double *)))==NULL)
      {
	fprintf(stderr, "Error: Could not allocate memory for crd.\n");
	exit(1);
      }
    
    for(i=0;i<numSamples;++i)
      {
	if( (crd[i]=(double *) calloc(3, sizeof(double)))==NULL)
	  {
	    fprintf(stderr, "Error: Could not allocate memory for crd.\n");
	    exit(1);
	  }
	crd[i][0]=seriesX[i];
	crd[i][1]=seriesY[i];
	crd[i][2]=seriesZ[i];
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
  int i, j, k;
  double shift[3];
  double delta;

  double **pbc_images;
  
  if( (pbc_images=(double **) calloc(nframes, sizeof(double *)))==NULL)
    {
      fprintf(stderr, "Error: Could not allocate memory for pbc_images.\n");
      exit(1);
    }
  
  for(i=0;i<nframes;++i)
    {
      if( (pbc_images[i]=(double *) calloc(3, sizeof(double)))==NULL)
	{
	  fprintf(stderr, "Error: Could not allocate memory for pbc_images.\n");
	  exit(1);
	}
      pbc_images[i][0]=0.0;
      pbc_images[i][1]=0.0;
      pbc_images[i][2]=0.0;
    }

  for(i=1;i<nframes;++i)
    {
      for(j=0;j<3;++j)
	{
	  delta=crd[i][j]-crd[i-1][j];
	  //	  cout << delta << endl;
	  if(delta > L/2)
	    {
	      pbc_images[i][j]-=L;
	    }
	  
	  if(delta<-L/2)
	    {
	      pbc_images[i][j]+=L;
	    }
	}
    }

  // accumlate pbc_images
  for(i=1;i<nframes;++i)
    {
      for(j=0;j<3;++j)
	{
	  pbc_images[i][j]+=pbc_images[i-1][j];
	}
    }
  for(i=1;i<nframes;++i)
    {
      for(j=0;j<3;++j)
	{
	  crd[i][j]+=pbc_images[i][j];
	}
    }
  
  for(i=0;i<nframes;++i)
    {
      delete[] pbc_images[i];
    }
  delete[] pbc_images;
}

double slope(double timestep, double *msd, int n)
{
  int i;
  double sum=0.0;

  for(i=1;i<n;++i)
    {
      sum+=(msd[i]/i);
    }

  printf("#sum = %lf timestep = %lf\n", sum, timestep); 
  return(sum/timestep/n);
}


int main(int argc, char **argv)
{
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

    int maxT=10000;

    double dt=DEFAULT_DELTA;

    double **crd;
    double *r_arr, *x_arr, *y_arr, *z_arr;
    
    int mflag=0;
    int stride=100;
    int nframes;
    int imaging=0;
    double L;
    int c;
    
    int *msd_counter;
    double *t, *msd;
    double D;
    char *traj_filename=NULL;
    
    opterr = 0;
    
    while ((c = getopt (argc, argv, "ime:T:t:p:s:d:l:")) != -1)
	switch (c)
	{
	    /* d, delta, timestep in ps */
	    case 'd':
	       dt=atof(optarg);
	       break;
	       /* i, turn on imaging for trajectories where this is off*/
	    case 'l':
	      L=atof(optarg);
	      break;
	    case 'i':
		imaging=1;
		break;
	    /* m, flag to include average cell dipole term (otherwise assumed to be zero) */
	    case 'm':
		mflag = 1;
		break;
	    /* T, temperature, temperature in K */
	    case 'T':
		T = atof(optarg);
		break;
	    case 's':
		stride=atoi(optarg);
		break;
	     case 't':
	       traj_filename = optarg;
	       break;
 	    case '?':
	     if (optopt == 't' || optopt == 'p' || optopt=='e')
		 fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	     else
	       {
		 fprintf (stderr, "Unknown option `-%c'.\n", optopt);
		 fprintf(stderr, "Syntax: %s -t traj_file [-d timestep (ps)] [-T temperature (K) \n", argv[0]);
		 
		 return 1;
	       }
	    default:
		abort ();
	}

    if(traj_filename==NULL)
      {
	cerr << "Error: Time series file name must be provided" << endl;
	exit(0);
      }
    
    crd=readSeriesNAMD(traj_filename, nframes);
    image(crd, nframes, L);
    
    cout << "#" << nframes << endl;

    if(nframes<0)
      {
	fprintf(stderr, "Error: Skipping more frames than there are in the trajectory.\n");
        exit(1);
      }

    count=0;

    if( (msd=(double *) calloc(maxT, sizeof(double)))==NULL)
      {
	fprintf(stderr, "Error: Could not allocate memory for mass array.\n");
	exit(1);
      }
    
    if( (msd_counter=(int *) calloc(maxT, sizeof(int)))==NULL)
      {
	fprintf(stderr, "Error: Could not allocate memory for msd_counter array.\n");
    	exit(1);
      }

    bzero(msd, maxT*sizeof(double));
    bzero(msd_counter, maxT*sizeof(int));

    for(i=0;i<nframes;i+=stride)
      {
	for(j=0;j<maxT;j++)
	  {
	    if( (i+j)>=nframes)
	      break;
	    
	    msd[j]+=calcmsd(crd[i], crd[i+j]);
	    msd_counter[j]=msd_counter[j]+1;
	  }
      }

    for(i=0;i<maxT;++i)
      {
	msd[i]/=msd_counter[i];
	cout << i << " " << msd[i] << endl;
		//, 
	//       pbc_images[i][0], pbc_images[i][1], pbc_images[i][2]);
      }

    
    D=slope(dt, msd, maxT)/stride/6.0;
    printf("#D %le %le (A2/fs)\n", D, slope(dt, msd, maxT));
    printf("#D %le (cm2/s)\n", D*1E-1);

    return(0);
}
