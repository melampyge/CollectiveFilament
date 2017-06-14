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
    cout << "Usage: " << argv[0] << "       infile name      n read       n msd       limit for MSD      n every      outfile name" << endl;
    return 0;
  }
  string fname = argv[1];
  int nall = atoi(argv[2]);
  int nmsd = atoi(argv[3]);
  int limit = atoi(argv[4]);
  int nevery = atoi(argv[5]);
  string fname2 = argv[6];
  double xi,yi,zi,xj,yj,zj;

  // declare variables, allocate arrays
  double *msd, *x, *y, *z,*t;
  int *count;
  int *logval;
  double lmin, lmax, dl;
  msd = (double*) malloc(nmsd*sizeof(double));
  count = (int*) malloc(nmsd*sizeof(int));
  logval = (int*) malloc(nmsd*sizeof(int));
  lmin = log(1);
  lmax = log(limit);
  dl = (lmax - lmin)/(nmsd-1);
  for (int i = 0; i < nmsd; i++) logval[i] = (int) exp(lmin + i*dl);
  for (int i = 0; i < nmsd; i++) {
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
  int kprint = 0;
  int nmax;
  for (int i = 0; i < nall; i++) {
    //if (i % nprint == 0) {
    //  cout <<"    "<<  kprint << " / " << nall << endl;
    //  kprint = kprint + nprint;
    //}
    xi = x[i];
    yi = y[i];
    zi = z[i];
    for (int l = 0; l < nmsd; l++) {
      int j = logval[l];
      if (j + i < nall) {
        xj = x[j+i];
        yj = y[j+i];
        zj = z[j+i];
        count[l] += 1;
        msd[l] += pow(xi-xj,2) + pow(yi-yj,2) + pow(zi-zj,2);
      }
    }
  }
  // normalize
  for (int i = 0; i <nmsd; i++) {
    msd[i] /= count[i];
  }


  // write results to file
  cout << "  writing results to file" << endl;
  ofstream ofile;
  ofile.open(fname2.c_str());
  ofile << "#time\tMSD\n";
  for (int i=0; i<nmsd; i++) {
    ofile << t[logval[i]] << "\t" <<  msd[i] << endl;
  }
  ofile.close();


  // free arrays
  free(x);
  free(y);
  free(z);
  free(msd);
  free(t);
  free(count);
  free(logval);


  cout << "  DONE" << endl;
  return 0;
}
