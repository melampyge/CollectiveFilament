#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
using namespace std; 

/*-----------------------------------------------------------------
Compute the mean square displacement based on a naive algorithm
use the input filename, the number of lines, and the limit of the
MSD as input parameters
-----------------------------------------------------------------*/
int main(int argc, char* argv[])
{ 
  // read in the information from the command line
  if (argc != 6) {
    cout << "Usage: " << argv[0] << " infile name      n read       n msr      limit for MSR      n every" << endl;
    return 0;
  }
  string fname = argv[1];
  int nall = atoi(argv[2]);
  int nmsr = atoi(argv[3]);
  int limit = atoi(argv[4]);
  int nevery = atoi(argv[5]);

  // declare variables, allocate arrays
  const double pi = 3.14159265359;
  double xi, yi, li, dphi;
  double x1,x2,y1,y2;
  double p1, p2;
  int nmax;

  double *x,*y,*phi,*phi2;
  double *msr, *t;
  int *count;
  int *logval;
  double di,dj;
  double lmin, lmax, dl;
  x = (double*) malloc(nall*sizeof(double));
  y = (double*) malloc(nall*sizeof(double));
  phi = (double*) malloc(nall*sizeof(double));
  phi2 = (double*) malloc(nall*sizeof(double));
  logval = (int*) malloc(nmsr*sizeof(int));
  msr = (double*) malloc(nmsr*sizeof(double));
  t = (double*) malloc(nall*sizeof(double));
  count = (int*) malloc(nmsr*sizeof(int));
  for (int i = 0; i < nmsr; i++) {
    msr[i] = 0.0;
    count[i] = 0;
  }
  lmin = log(1);
  lmax = log(limit);
  dl = (lmax - lmin)/(nmsr-1);
  for (int i = 0; i < nmsr; i++) logval[i] = (int) exp(lmin + i*dl);


  // read in the data
  cout << "  reading coordinates" <<  endl;
  ifstream infile;
  infile.open(fname.c_str());
  // error check
  if (!infile.is_open()) {
    cout << "Error, could not open file \n";
    return 0;
  }
  // skip the first two lines 
  string line;
  double tmp;
  getline(infile,line);
  getline(infile,line);
  // read in line by line;
  for (int i = 0; i < nall; i++) {
    infile >> t[i] >> x[i] >> y[i];
    for (int j = 0; j < nevery-1; j++) {
      infile >> tmp >> tmp;
    }
  }
  infile.close();

  // normalize x and y
  cout << "  normalizing x and y" << endl;
  for (int i = 0; i < nall; i++) {
    xi = x[i];
    yi = y[i];
    li = sqrt(xi*xi + yi*yi);
    x[i] = xi/li;
    y[i] = yi/li;
  }

  // compute local angle using the atan2 function
  cout << "  computing angles" << endl;
  for (int i = 0; i < nall; i++) {
    xi = x[i];
    yi = y[i];
    phi[i] = atan2(yi,xi);
  }

  // correct for periodicity when computing phi2
  cout << "  correcting periodicity" << endl;
  phi2[0] = phi[0];
  for (int i = 1; i < nall; i++) {
    dphi = phi[i] - phi[i-1];
    if (dphi < -pi) dphi += 2*pi;
    if (dphi > pi) dphi -= 2*pi;
    phi2[i] = phi2[i-1] + dphi;
  }

 

  // compute the mean square rotation
  cout << "  computing the msr" << endl;
  for (int i = 0; i < nall; i++) {
    p1 = phi2[i];
    for (int l = 0; l < nmsr; l++) {
      int j = logval[l];
      if (j + i < nall) { 
        p2 = phi2[i+j];
        dphi = (p2-p1)*(p2-p1);
        count[l] += 1;
        msr[l] += dphi;
      }
    }
  }
  cout << " Normalizing" << endl;
  // normalize
  for (int i = 0; i < nmsr; i++) {
    if (count[i] > 0) msr[i] /= count[i];
  }

  cout << " writing Results to file" << endl;
  // write results to files
  ofstream ofile;
  ofile.open("e2e_msr.data");
  ofile << "#time\tmsr\n";
  for (int i=0; i<nmsr; i++) {
    ofile << t[logval[i]] << "\t" <<  msr[i] << endl;
  }
  ofile.close();

  // free arrays
  free(x);
  free(y);
  free(phi);
  free(phi2);
  free(msr);
  free(t);
  free(count);
  free(logval);

  cout << "  DONE" << endl;
  return 0;
}
