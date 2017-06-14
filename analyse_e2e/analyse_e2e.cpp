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
  if (argc != 5) {
    cout << "Usage: " << argv[0] << " infile name      n read       limit for acf      n every" << endl;
    return 0;
  }
  string fname = argv[1];
  int nall = atoi(argv[2]);
  int limit = atoi(argv[3]);
  int nevery = atoi(argv[4]);

  // declare variables, allocate arrays
  const double pi = 3.14159265359;
  double xi, yi, li, dphi;
  double x1,x2,y1,y2;
  double p1, p2;
  int nmax;

  double *x,*y,*phi,*phi2;
  double *acf, *msr, *t;
  int *count;
  double di,dj;
  x = (double*) malloc(nall*sizeof(double));
  y = (double*) malloc(nall*sizeof(double));
  phi = (double*) malloc(nall*sizeof(double));
  phi2 = (double*) malloc(nall*sizeof(double));
  msr = (double*) malloc(limit*sizeof(double));
  t = (double*) malloc(limit*sizeof(double));
  acf = (double*) malloc(limit*sizeof(double));
  count = (int*) malloc(limit*sizeof(int));
  for (int i = 0; i < limit; i++) {
    msr[i] = 0.0;
    acf[i] = 0.0;
    count[i] = 0;
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

  
  // compute the acf
  cout << "  computing the acf" << endl;
  int nprint = nall/1000;
  int k = 0;
  for (int i = 0; i < nall; i++) {
    if (i % nprint == 0) {
      cout <<"    "<<  k << " / " << nall << endl;
      k = k + nprint;
    }
    x1 = x[i];
    y1 = y[i];
    nmax = min(limit, nall - i);
    for (int j=0; j < nmax; j++) {
      x2 = x[j+i];
      y2 = y[j+i];
      count[j] += 1;
      dphi = x1*x2 + y1*y2;
      acf[j] += dphi;
    }
  }
  // normalize
  for (int i = 0; i < limit; i++) {
    acf[i] /= count[i];
  }

  // reset counter
  for (int i = 0; i < limit; i++) {
    count[i] = 0;
  }

  // compute the mean square rotation
  cout << "  computing the msr" << endl;
  k = 0;
  for (int i = 0; i < nall; i++) {
    if (i % nprint == 0) {
      cout <<"    "<<  k << " / " << nall << endl;
      k = k + nprint;
    }
    p1 = phi2[i];
    nmax = min(limit, nall - i);
    for (int j=0; j < nmax; j++) {
      p2 = phi2[i+j];
      dphi = (p2-p1)*(p2-p1);
      count[j] += 1;
      msr[j] += dphi;
    }
  }
  // normalize
  for (int i = 0; i < limit; i++) {
    msr[i] /= count[i];
  }

  // write results to files
  cout << "  writing results to file" << endl;
  ofstream ofile;
  ofile.open("e2e_acf.data");
  ofile << "#time\tacf\n";
  for (int i=0; i<limit; i++) {
    ofile << t[i] << "\t" <<  acf[i] << endl;
  }
  ofile.close();
  ofstream ofile2;
  ofile2.open("e2e_msr.data");
  ofile2 << "#time\tmsr\n";
  for (int i=0; i<limit; i++) {
    ofile2 << t[i] << "\t" <<  msr[i] << endl;
  }
  ofile2.close();
  // for testing purposes, write phi2 and scaled x and y
  ofstream ofile3;
  ofile3.open("e2e_test.data");
  ofile3 << "#time\txs\tys\tphi\tphi2\n";
  for (int i=0; i<nall; i++) {
    ofile3 << t[i] << "\t" << x[i] << "\t" << y[i] << "\t" <<  phi[i]
           << "\t" << phi2[i] << endl;
  }
  ofile3.close();


  // free arrays
  free(x);
  free(y);
  free(phi);
  free(phi2);
  free(acf);
  free(msr);
  free(t);
  free(count);

  cout << "  DONE" << endl;
  return 0;
}
