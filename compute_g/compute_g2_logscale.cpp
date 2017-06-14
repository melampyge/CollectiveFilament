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
  if (argc != 10) {
    cout << "Usage: " << argv[0] << "      infile name       comfile      n read     n msd      limit for MSD      n atoms      selected atom     nevery      outfile name" << endl;
    return 0;
  }
  string fname = argv[1];
  string cname = argv[2];
  int nall = atoi(argv[3]);
  int nmsd = atoi(argv[4]);
  int limit = atoi(argv[5]);
  int natoms = atoi(argv[6]);
  int nselected = atoi(argv[7]);
  int nevery = atoi(argv[8]);
  string fname2 = argv[9];

  // declare variables, allocate arrays
  double *msd, *x, *y, *z,*t;
  double *comx,*comy,*comz,*comt;
  int *count;
  int *logval;
  float lmin,lmax,dl;
  double xi,yi,zi,xj,yj,zj;
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
  comx = (double*) malloc(nall*sizeof(double));
  comy = (double*) malloc(nall*sizeof(double));
  comz = (double*) malloc(nall*sizeof(double));
  comt = (double*) malloc(nall*sizeof(double));


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
      if (j == nselected-1) infile >> tmp >> x[i] >> y[i] >> z[i];
      else infile >> tmp >> tmp >> tmp >> tmp;
    }
  }
  infile.close();

  // read in the com
  ifstream infile2;
  infile2.open(cname.c_str());
  // error check
  if (!infile2.is_open()) {
    cout << "Error, could not open com file\n";
    return 0;
  }
  // read in all data, skip first two lines and nevery-1
  string line;
  double tmp;
  getline(infile2,line);
  getline(infile2,line);
  for (int i = 0; i < nall; i++) {
    infile2 >> comt[i] >> comx[i] >> comy[i] >> comz[i];
    for (int j = 0; j < nevery-1; j++) {
      infile2 >> tmp >> tmp >> tmp >> tmp;
    }
  }
  infile2.close();

  // compute the MSD
  cout << "  computing the MSD" << endl;
  int nprint = nall/1000;
  int kprint = 0;
  for (int i = 0; i < nall; i++) {
    if (i % nprint == 0) {
      cout <<"    "<<  kprint << " / " << nall << endl;
      kprint = kprint + nprint;
    }
    xi = x[i]-comx[i];
    yi = y[i]-comy[i];
    zi = z[i]-comz[i];
    for (int l=0; l < nmsd; l++) {
      int j = logval[l];
      if (j + i < nall) {
        xj = x[j+i]-comx[j+i];
        yj = y[j+i]-comy[j+i];
        zj = z[j+i]-comz[j+i];
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
  free(comx);
  free(comy);
  free(comz);
  free(comt);

  cout << "  DONE" << endl;
  return 0;
}
