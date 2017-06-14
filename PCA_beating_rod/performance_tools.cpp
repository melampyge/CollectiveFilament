#include <iostream>
#include <math.h>

using namespace std;

class Performance_tools{
  public:

  void compute_curvature(double* curv, double* xall, double* yall,
			 int nsteps, int natoms) {
    // declare some required variables
    int idx, idx2;
    double x0, x1, x2, y0, y1, y2;
    double ax, ay, bx, by, cx, cy;
    double A, a, b, c;
    
    // compute the local curvature
    for (int i = 0; i < nsteps; i++) {
      for (int j = 0; j < natoms - 2; j++) {
        idx = natoms*i + j;
	idx2 = (natoms-2)*i + j;
        x0 = xall[idx];
	x1 = xall[idx+1];
	x2 = xall[idx+2];
	y0 = yall[idx];
	y1 = yall[idx+1];
	y2 = yall[idx+2];

	// vectors between the points
	ax = x0-x1;
	ay = y0-y1;
	bx = x2-x1;
	by = y2-y1;
	cx = x2-x0;
	cy = y2-y0;
	// size of the vectors
	a = sqrt(ax*ax + ay*ay);
	b = sqrt(bx*bx + by*by);
	c = sqrt(cx*cx + cy*cy);
	// area
	A = 0.5*(ax*by - ay*bx);
	// curvature
	curv[idx2] = 4.*A/(a*b*c);
      }
    }

    // return
    return;
  }
 
  /***********************************************************************************/

  void compute_curvature_single(double* curv, double* x, double* y, int natoms) {

    // declare some required variables
    double x0, x1, x2, y0, y1, y2;
    double ax, ay, bx, by, cx, cy;
    double A, a, b, c;
    
    // compute the local curvature
    for (int i = 0; i < natoms - 2; i++) {
      x0 = x[i];
      x1 = x[i+1];
      x2 = x[i+2];
      y0 = y[i];
      y1 = y[i+1];
      y2 = y[i+2];

      // vectors between the points
      ax = x0-x1;
      ay = y0-y1;
      bx = x2-x1;
      by = y2-y1;
      cx = x2-x0;
      cy = y2-y0;
      // size of the vectors
      a = sqrt(ax*ax + ay*ay);
      b = sqrt(bx*bx + by*by);
      c = sqrt(cx*cx + cy*cy);
      // area
      A = 0.5*(ax*by - ay*bx);
      // curvature
      curv[i] = 4.*A/(a*b*c);
    }

    // return
    return;
  }
 
  /***********************************************************************************/

  void compute_crosscorrelation(int* t_cc, int* time, double* a1, double* a2, double* cc1, double* cc2,
			       int* counter, int* linval,
			       int nsteps, int ncc, int limit, int dt) {
    // allocate variables
    double a1i, a1j, a2i, a2j;
    
    // zero cc and counter
    for (int i = 0; i < ncc; i++) cc1[i] = 0.0;
    for (int i = 0; i < ncc; i++) cc2[i] = 0.0;
    for (int i = 0; i < ncc; i++) counter[i] = 0;

    // fill linval array
    for (int i = 0; i < ncc; i++) linval[i] = i*dt;

    // fill t_cc array
    for (int i = 0; i < ncc; i++) t_cc[i] = time[linval[i]] - time[0];
  
    // compute the autocorrelation
    for (int i = 0; i < nsteps; i++) {
      a1i = a1[i];
      a2i = a2[i];
      for (int l = 0; l < ncc; l++) {
        int j = linval[l];
        if (j + i < nsteps) {
          a1j = a1[j+i];
	  a2j = a2[j+i];
          counter[l] += 1;
          cc1[l] += a1i*a2j;
	  cc2[l] += a2i*a1j;
        }
      }
    }
    
    // normalize autocorrelation
    for (int i = 0; i < ncc; i++) cc1[i] /= counter[i];
    for (int i = 0; i < ncc; i++) cc2[i] /= counter[i];

    // return
    return;
  }


  /***********************************************************************************/

  void blocking_method(double* blockdata, double* block_std, double* block_uncert, int n, double& acl, double& stdl) {

    // declare variables
    double c0;
    double smin, smax;
    double smin_old, smax_old;
    
    // compute the average
    acl = 0;
    for (int i = 0; i < n; i++) acl += blockdata[i];
    acl /= n;

    // perform the blocking method
    int j = 0;
    while (n >= 2) {
      // compute std for current block
      c0 = 0;
      for (int i = 0; i < n; i++) c0 += pow(blockdata[i] - acl, 2);
      c0 /= (n-1);
      // compute the standard deviation and fill it to blocking array
      c0 = sqrt(c0/(n-1));
      block_std[j] = c0;
      block_uncert[j] = c0/sqrt(2*(n-1));
      // perform block operation
      for (int i = 0; i < n/2; i++) blockdata[i] = 0.5*(blockdata[2*i] + blockdata[2*i+1]);
      n = n/2;
      // increase j
      j = j + 1;
    }
    
    // find the plateau value and fill it to stdl, set to -1 if no plateau is found
    // start from the lowermost value, ignore last nignore entries
    int nignore = 3;
    n = j - nignore -1;
    smin_old = block_std[0] - 2*block_uncert[0];
    smax_old = block_std[0] + 2*block_uncert[0];
    stdl = -1;
    for (int i = 1; i < n; i++) {
      smin = block_std[i] - 2*block_uncert[i];
      smax = block_std[i] + 2*block_uncert[i];
      if (smax_old > smin) {
	stdl = 0.5*(block_std[i] + block_std[i-1]);
	break;
      }
      smin_old = smin;
      smax_old = smax;
    }
    
    return;
  }
  
  /***********************************************************************************/

  void compute_crosscorrelation_with_errorbars(int* t_cc, int* time, double* x1, double* x2, double* cc,
					      double* std, double* blockdata, double* block_std,
					      double* block_uncert, int* linval,
					      int nsteps, int ncc, int limit, int dt) {
    // allocate variables
    double xi, xj;
    double ccl, stdl;
    
    // zero cc and std
    for (int i = 0; i < ncc; i++) cc[i] = 0.0;
    for (int i = 0; i < ncc; i++) std[i] = 0.0;

    // fill linval array
    for (int i = 0; i < ncc; i++) linval[i] = i*dt;

    // fill t_cc array
    for (int i = 0; i < ncc; i++) t_cc[i] = time[linval[i]] - time[0];

    // compute autocorrelation function with statistical uncertainties
    //   using the blocking method:
    //   have l as the outer loop
    //   compute all values for fixed l and store them in an array
    //   apply blocking method to this array; find true std automatically

    for (int l = 0; l < ncc; l++) {
      int j = linval[l];
      cout << "   " << l << " " << ncc << endl;
      for (int i = 0; i < nsteps-j; i++) {
        xi = x1[i];
	xj = x2[j+i];
	blockdata[i] = xi*xj;	
      }
      // compute mean and true standard deviation from the blocking method
      blocking_method(blockdata, block_std, block_uncert, nsteps-j, ccl, stdl);
      cc[l] = ccl;
      std[l] = stdl;
    }

    // return
    return;
  }

  
  
};

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/


extern "C" {
  Performance_tools* Performance_tools_new(){return new Performance_tools();}

  /***********************************************************************************/

  void compute_curvature(Performance_tools* performance_tools,
			 double* curv, double* xall, double* yall,
			 int nsteps, int natoms) {
    performance_tools->compute_curvature(curv, xall, yall, nsteps, natoms);
  }

  /***********************************************************************************/

  void compute_curvature_single(Performance_tools* performance_tools,
			 double* curv, double* x, double* y,
			 int natoms) {
    performance_tools->compute_curvature_single(curv, x, y, natoms);
  }

  /***********************************************************************************/

  void compute_crosscorrelation(Performance_tools* performance_tools,
			        int* t_cc, int* time, double* a1, double* a2, double* cc1, double* cc2,
       			        int* counter, int* linval,
			        int nsteps, int ncc, int limit, int dt) {
    performance_tools->compute_crosscorrelation(t_cc, time, a1, a2, cc1, cc2, counter, linval, nsteps, ncc, limit, dt);
  }

  /***********************************************************************************/

  void compute_crosscorrelation_with_errorbars(Performance_tools* performance_tools,
		 			       int* t_cc, int* time, double* x1, double* x2, double* cc, double* std,
		 			       double* blockdata, double* block_std, double* block_uncert, int* linval,
		 			       int nsteps, int ncc, int limit, int dt) {
    performance_tools->compute_crosscorrelation_with_errorbars(t_cc, time, x1, x2, cc, std, blockdata, block_std, block_uncert, linval, nsteps, ncc, limit, dt);
  }
  
}

