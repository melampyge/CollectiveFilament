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
    cout << "Usage: " << argv[0] << " infile name      n read       limit for MSD      n every      outfile name" << endl;
    return 0;
  }
  string fname = argv[1];
  int nall = atoi(argv[2]);
  int limit = atoi(argv[3]);
  int nevery = atoi(argv[4]);
  string fname2 = argv[5];

  // declare variables, allocate arrays
  double *msd, *x, *y, *z,*t;
  int *count;
  double xi,yi,zi,xj,yj,zj;
  msd = (double*) malloc(limit*sizeof(double));
  count = (int*) malloc(limit*sizeof(int));
  for (int i = 0; i < limit; i++) {
    msd[i] = 0.0;
    count[i] = 0;
  }
  t = (double*) malloc(nall*sizeof(double));
  x = (double*) malloc(nall*sizeof(double));
  y = (double*) malloc(nall*sizeof(double));
  z = (double*) malloc(nall*sizeof(double));

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
    infile >> t[i] >> x[i] >> y[i] >> z[i];
    for (int j = 0; j < nevery-1; j++) {
      infile >> tmp >> tmp >> tmp >> tmp;
    }
  }
  infile.close();

  // compute the MSD
  cout << "  computing the MSD" << endl;
  int nprint = nall/1000;
  int k = 0;
  int nmax;
  for (int i = 0; i < nall; i++) {
    if (i % nprint == 0) {
      cout <<"    "<<  k << " / " << nall << endl;
      k = k + nprint;
    }
    xi = x[i];
    yi = y[i];
    zi = z[i];
    nmax = min(limit, nall - i);
    for (int j=1; j < nmax; j++) {
      xj = x[j+i];
      yj = y[j+i];
      zj = z[j+i];
      count[j] += 1;
      msd[j] += pow(xi-xj,2) + pow(yi-yj,2) + pow(zi-zj,2);
    }
  }
  // normalize
  for (int i = 1; i <limit; i++) {
    msd[i] /= count[i];
  }

  // write results to file
  cout << "  writing results to file" << endl;
  ofstream ofile;
  ofile.open(fname2.c_str());
  ofile << "#time\tMSD\n";
  for (int i=0; i<limit; i++) {
    ofile << t[i] << "\t" <<  msd[i] << endl;
  }
  ofile.close();

  // free arrays
  free(x);
  free(y);
  free(z);
  free(msd);
  free(t);
  free(count);

  cout << "  DONE" << endl;
  return 0;
}
