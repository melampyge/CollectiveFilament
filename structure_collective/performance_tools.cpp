
#include <iostream>
#include <cmath>

#define pi M_PI

using namespace std;

class Performance_tools {
  public:

/***********************************************************************************/

  void compute_static_structure(double *S, double *x, double *y, int nsteps, int natoms, double box_size) {

    double delk = 2*pi/box_size;
    int N = static_cast<int>(box_size/2);
    int N2 = N*N;
   
    for (int step = 0; step < nsteps; step++) {

	for (int nx = 0; nx < N; nx++) {
	  double kx = delk*nx;

	    for (int ny = 0; ny < N; ny++) {
	      double ky = delk*ny;
	      int k = static_cast<int>(sqrt(kx*kx + ky*ky));
	      double costerm = 0;
	      double sinterm = 0;

	      for (int j = 0; j < natoms; j++) {
		double dotp = kx*x[step*natoms+j] + ky*y[step*natoms+j];
		costerm += cos(dotp);
		sinterm += sin(dotp);
	
	      } // particle loop

	    S[k] += costerm*costerm + sinterm*sinterm;

	    } // ky loop
	} // kx loop        
    } // time loop

    for (int j = 0; j < N2; j++) { S[j] /= (natoms*nsteps); }

  }
    
};

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/


extern "C" {

  Performance_tools* Performance_tools_new() { return new Performance_tools(); }

  void compute_static_structure(Performance_tools* performance_tools,
		   double *S, double *x, double *y,
		   int nsteps, int natoms, double box_size) {
    performance_tools->compute_static_structure(S, x, y, nsteps, natoms, box_size);
  }
}

