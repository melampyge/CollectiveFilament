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
  if (argc != 9) {
    cout << "Usage: " << argv[0] << " infile name      n read       n atoms       n kvec      n time      nlimit     nevery      outfile name" << endl;
    return 0;
  }
  string fname = argv[1];
  int nall = atoi(argv[2]);
  int natoms = atoi(argv[3]);
  int nkvec = atoi(argv[4]);
  int nt = atoi(argv[5]);
  int limit = atoi(argv[6]);
  int nevery = atoi(argv[7]);
  string fname2 = argv[8];

  // allocate some variables
  double delx,dely,delz;
  double kmin, kmax;
  double dk;
  const double pi = 3.14159265359;
  
  // allocate arrays
  double *x,*y;            // positions of the taoms
  double *t;               // time 
  double *kvec;            // absolute value of k-vectors
  double *kxs, *kys;       // circular orientation of k-vectors
  double *sqtall, *sqtall2;          // structure factor as a function of all q and t
  double *sqt;             // structure factor as a function of average q and t
  int *logval;             // integer with local time values
  double lmin, lmax, dl;
  int nks = 24;
  double phi;
  double xi, yi, xj, yj, q1, q2, dx, dy;
  double expterm, term;
  double cterm, sterm;
  int i,j,k,l,m,n,idx,idx2;
 
  // allocate arrays
  x = (double*) malloc(natoms*nall*sizeof(double));
  y = (double*) malloc(natoms*nall*sizeof(double));
  t = (double*) malloc(nall*sizeof(double));
  kvec = (double*) malloc(nkvec*sizeof(double));
  kxs = (double*) malloc(nks*sizeof(double));
  kys = (double*) malloc(nks*sizeof(double));
  sqtall = (double*) malloc(nks*nkvec*nt*sizeof(double));
  sqtall2 = (double*) malloc(nks*nkvec*nt*sizeof(double));
  for (i = 0; i < nks*nkvec*nt; i++) {
    sqtall[i] = 0.0;
    sqtall2[i] = 0.0;
  }
  sqt = (double*) malloc(nkvec*nt*sizeof(double));
  for (i = 0; i < nkvec*nt; i++) sqt[i] = 0.0;
  logval = (int*) malloc(nt*sizeof(int));
  lmin = log(1);
  lmax = log(limit);
  dl = (lmax - lmin)/(nt-1);
  for (i = 0; i < nt; i++) logval[i] = int(exp(lmin + i*dl));
  logval[0] = 0;
  

  // reading in the coordinates and time
  cout << "  reading in all coordinates" << endl;
  ifstream infile;
  infile.open(fname.c_str());
  // error check
  if (!infile.is_open()) {
    cout << "Error, could not open file \n";
    return 0;
  }
  // loop over all snapshots
  string line;
  double tmp;
  for (i = 0; i < nall; i++) {
    infile >> line;
    infile >> line >> line >> t[i];
    for (j = 0; j < natoms; j++) {
      idx = i*natoms + j;
      infile >> tmp >> x[idx] >> y[idx] >> tmp;
    }
  }
  // close the file
  infile.close();

  // compute the fourier transform of the mean square displacements = S(q,t)
  cout << "  computing the structure factor" << endl;
  kmin = 0.01;
  kmax = 10.0;
  kmin = log(kmin);
  kmax = log(kmax);
  dk = (kmax - kmin)/nkvec;
  for (i = 0; i < nkvec; i++) kvec[i] = kmin + i*dk;
  for (i = 0; i < nkvec; i++) kvec[i] = exp(kvec[i]);
  // HARD CODE THE DESIRED VALUE FOR kvec:
  kvec[0] = 0.158489;
  for (i = 0; i < nks; i++) {
    phi = 2.0*pi/(nks)*i;
    kxs[i] = cos(phi);
    kys[i] = sin(phi);
  }
  // start the loops
  for (n = 0; n < nall-limit; n++) { // snapshots
    if (n % nevery != 0) continue;
    cout << "    n = " << n << endl;
    for (i = 0; i < natoms; i++) {  // atom i
      for (j = 0; j < natoms; j++) {  // atom j
        for (k = 0; k < nt; k++) { // time
          for (l = 0; l < nkvec; l++) { // absolute value of the k-vector
	    for (m = 0; m < nks; m++) { // orientation of the k-vector
	      idx = n*natoms + i;
              xi = x[idx];
	      yi = y[idx];
	      idx = (n+logval[k])*natoms + j;
	      xj = x[idx];
              yj = y[idx];
	      dx = xi - xj;
	      dy = yi - yj;
	      q1 = kvec[l]*kxs[m];
	      q2 = kvec[l]*kys[m];
	      expterm = dx*q1 + dy*q2;
	      cterm = cos(expterm);
	      //sterm = sin(expterm);
              idx = k*nkvec*nks + l*nks + m;
	      sqtall[idx] += cterm;
	      //sqtall2[idx] += sterm;
	    }
          }
        }
      }
    }
  }
  // include normalization factors
  for (i = 0; i < nt*nkvec*nks; i++) sqtall[i] /= natoms*((nall-limit)/nevery);
  // average to get sqt
  for (k = 0; k < nt; k++) {
    for (l = 0; l < nkvec; l++) {
      for (m = 0; m < nks; m++) {
	sqt[k*nkvec + l] += sqtall[k*nkvec*nks + l*nks + m];
      }
    }  
  }
  for (i = 0; i < nt*nkvec; i++) sqt[i] /= nks;


  // write results to file
  cout << "  writing results to file" << endl;
  ofstream ofile;
  ofile.open(fname2.c_str());
  ofile << "#k\tt\tS\n";
  for (i = 0; i < nkvec; i++) {
    for (j = 0; j < nt; j++) {
      ofile << kvec[i] << "\t" << t[logval[j]] << "\t" << sqt[j*nkvec + i] << endl;
    }
  }
  ofile.close();

  // free arrays
  free(x);
  free(y);
  free(t);
  free(kvec);
  free(kxs);
  free(kys);
  free(sqtall);
  free(sqt);
  free(logval);

  cout << "  DONE" << endl;
  return 0;
}
