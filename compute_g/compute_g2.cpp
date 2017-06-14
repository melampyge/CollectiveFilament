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
  if (argc != 7) {
    cout << "Usage: " << argv[0] << " infile name      n read       limit for MSD      n atoms      selected atom     outfile name" << endl;
    return 0;
  }
  string fname = argv[1];
  int nall = atoi(argv[2]);
  int limit = atoi(argv[3]);
  int natoms = atoi(argv[4]);
  int nselected = atoi(argv[5]);
  string fname2 = argv[6];

  // declare variables, allocate arrays
  double *msd, *x, *y, *z,*t;
  double *comx, *comy, *comz;
  int *count;
  double xi,yi,zi,xj,yj,zj;
  msd = (double*) malloc(limit*sizeof(double));
  count = (int*) malloc(limit*sizeof(int));
  for (int i = 0; i < limit; i++) {
    count[i] = 0;
    msd[i] = 0.0;
  }
  t = (double*) malloc(nall*sizeof(double));
  x = (double*) malloc(nall*sizeof(double));
  y = (double*) malloc(nall*sizeof(double));
  z = (double*) malloc(nall*sizeof(double));
  comx = (double*) malloc(nall*sizeof(double));
  comy = (double*) malloc(nall*sizeof(double));
  comz = (double*) malloc(nall*sizeof(double));
  for (int i = 0; i < nall; i++) {
    comx[i] = 0.0;
    comy[i] = 0.0;
    comz[i] = 0.0;
  }
  

  // read in the data
  cout << "  reading coordinates" <<  endl;
  ifstream infile;
  infile.open(fname.c_str());
  // error check
  if (!infile.is_open()) {
    cout << "Error, could not open file \n";
    return 0;
  }
  // loop over all snapshots
  for (int i = 0; i < nall; i++) {
    // skip the first lines
    string line;
    double tmp;
    infile >> line;
    // read in the time
    infile >> line >> line >> t[i];
    for (int j = 0; j < natoms; j++) {
      infile >> tmp >> xi >> yi >> zi;
      if (j == nselected-1) {
	x[i] = xi;
	y[i] = yi;
	z[i] = zi;
      }
      comx[i] += xi;
      comy[i] += yi;
      comz[i] += zi;
    }
    comx[i] /= natoms;
    comy[i] /= natoms;
    comz[i] /= natoms;
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
    xi = x[i]-comx[i];
    yi = y[i]-comy[i];
    zi = z[i]-comz[i];
    nmax = min(limit, nall - i);
    for (int j=1; j < nmax; j++) {
      xj = x[j+i]-comx[j+i];
      yj = y[j+i]-comy[j+i];
      zj = z[j+i]-comz[j+i];
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
  free(count);
  free(comx);
  free(comy);
  free(comz);
  free(t);

  cout << "  DONE" << endl;
  return 0;
}
