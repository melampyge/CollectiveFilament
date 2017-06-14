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
    cout << "Usage: " << argv[0] << " infile name      n read       n atoms   n kvec        outfile name" << endl;
    return 0;
  }
  string fname = argv[1];
  int nall = atoi(argv[2]);
  nall = 400000; //MANUAL CORRECTION BECAUSE OF BAD DESIGN OF EXPERIMENT
  int natoms = atoi(argv[3]);
  int nkvec = atoi(argv[4]);
  string fname2 = argv[5];

  // allocate some variables
  double comx, comy, comz;
  double delx,dely,delz;
  double rgyrav;
  double kmin, kmax;
  double dk;
  const double pi = 3.14159265359;
  
  // allocate arrays
  double *x,*y,*z,*rgyr;
  double *kvec;
  double *kxs, *kys;
  double *sq;
  int nks = 24;
  double phi;
  double xi, yi, q1, q2;
  double expterm, term;
  double cterm, sterm;
 
  
  x = (double*) malloc(natoms*sizeof(double));
  y = (double*) malloc(natoms*sizeof(double));
  z = (double*) malloc(natoms*sizeof(double));
  rgyr = (double*) malloc(nall*sizeof(double));
  for (int i = 0; i < nall; i++) rgyr[i] = 0.0;
  kvec = (double*) malloc(nkvec*sizeof(double));
  kxs = (double*) malloc(nks*sizeof(double));
  kys = (double*) malloc(nks*sizeof(double));
  sq = (double*) malloc(nkvec*sizeof(double));
  for (int i = 0; i < nkvec; i++) sq[i] = 0.0;
  
  // compute the average radius of gyration
  // declare variables, allocate arrays
  cout << "  computing the radius of gyration" <<  endl;
  int nprint = nall/1000;
  int kprint = 0;
  int nmax;
  ifstream infile;
  infile.open(fname.c_str());
  // error check
  if (!infile.is_open()) {
    cout << "Error, could not open file \n";
    return 0;
  }
  // loop over all snapshots
  for (int i = 0; i < nall; i++) {
    //write progress to screen
    if (i % nprint == 0) {
      cout <<"    "<<  kprint << " / " << nall << endl;
      kprint = kprint + nprint;
    }
    // skip the first two lines
    string line;
    double tmp;
    infile >> line;
    infile >> line >> line >> line;
    for (int j = 0; j < natoms; j++) {
      infile >> tmp >> x[j] >> y[j] >> z[j];
    }
    // comute com
    comx = 0.0;
    comy = 0.0;
    comz = 0.0;
    for (int j = 0; j < natoms; j++) {
      comx += x[j];
      comy += y[j];
      comz += z[j];
    }
    comx /= natoms;
    comy /= natoms;
    comz /= natoms;
    // compute rgyr
    for (int j = 0; j < natoms; j++) {
      delx = x[j] - comx;
      dely = y[j] - comy;
      delz = z[j] - comz;
      rgyr[i] += delx*delx + dely*dely + delz*delz;
    }
    rgyr[i] /= natoms;
  }
  infile.close();

  // compute the average of rgyr**2
  rgyrav = 0.0;
  for (int i = 0; i < nall; i++) rgyrav += rgyr[i];
  rgyrav /= nall;
  rgyrav = sqrt(rgyrav);
  cout << "  rgyrav = " << rgyrav << endl;


  cout << "  computing the structure factor" <<  endl;
  // define the k vectors
  double incfac = 2.0; // how much more kvectors to use
  //kmin = 2*pi/rgyrav/incfac/100;
  //kmax = 2*pi*incfac;
  kmin = 0.01;
  kmax = 100.0;
  kmin = log(kmin);
  kmax = log(kmax);
  dk = (kmax - kmin)/nkvec;
  for (int i = 0; i < nkvec; i++) kvec[i] = kmin + i*dk;
  for (int i = 0; i < nkvec; i++) kvec[i] = exp(kvec[i]);
  for (int i = 0; i < nks; i++) {
    phi = 2*pi/(nks)*i;
    kxs[i] = cos(phi);
    kys[i] = sin(phi);
  }

  ifstream infile2;
  infile2.open(fname.c_str());
  // error check
  if (!infile2.is_open()) {
    cout << "Error, could not open file \n";
    return 0;
  }
  // loop over all snapshots to compute SQ
  kprint = 0;
  for (int i = 0; i < nall; i++) {
    // write down the progress to the screen
    if (i % nprint == 0) {
      cout <<"    "<<  kprint << " / " << nall << endl;
      kprint = kprint + nprint;
    }
    // skip the first two lines
    string line;
    double tmp;
    infile2 >> line;
    infile2 >> line >> line >> line;
    for (int j = 0; j < natoms; j++) {
      infile2 >> tmp >> x[j] >> y[j] >> z[j];
    }
    // compute the structure factor
    for (int k = 0; k < nkvec; k++) {
      term = 0.0;
      for (int l = 0; l < nks; l++) {
	q1 = kvec[k]*kxs[l];
	q2 = kvec[k]*kys[l];	
	cterm = 0.0;
	sterm = 0.0;
        for (int j = 0; j < natoms; j++) {
          xi = x[j];
          yi = y[j];

	  expterm = xi*q1 + yi*q2;
	  cterm += cos(expterm);
	  sterm += sin(expterm);
	}
	term += cterm*cterm + sterm*sterm;
      }
      term /= nks;
      sq[k] += term;
    }
  }
  infile2.close();
  for (int k = 0; k < nkvec; k++) sq[k] /= nall*natoms;

  // write results to file
  cout << "  writing results to file" << endl;
  ofstream ofile;
  ofile.open(fname2.c_str());
  ofile << "#k\tS\n";
  for (int i=0; i< nkvec; i++) {
    ofile << kvec[i] << "\t" <<  sq[i] << endl;
  }
  ofile.close();

  // free arrays
  free(x);
  free(y);
  free(z);
  free(rgyr);
  free(kvec);
  free(kxs);
  free(kys);
  free(sq);

  cout << "  DONE" << endl;
  return 0;
}
