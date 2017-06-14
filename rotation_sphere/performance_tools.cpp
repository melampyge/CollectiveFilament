#include <iostream>
#include <math.h>

using namespace std;

class Performance_tools{
  public:

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

  void compute_ttac_with_errorbars(int* t_ac, int* time, double* tx, double*ty, double* ac,
					      double* std, double* blockdata, double* block_std,
					      double* block_uncert, int* linval,
					      int nsteps, int nac, int limit, int dt) {
    // allocate variables
    double txi, txj, tyi, tyj;
    double ccl, stdl;
    
    // zero cc and std
    for (int i = 0; i < nac; i++) ac[i] = 0.0;
    for (int i = 0; i < nac; i++) std[i] = 0.0;

    // fill linval array
    for (int i = 0; i < nac; i++) linval[i] = i*dt;

    // fill t_cc array
    for (int i = 0; i < nac; i++) t_ac[i] = time[linval[i]] - time[0];

    // compute autocorrelation function with statistical uncertainties
    //   using the blocking method:
    //   have l as the outer loop
    //   compute all values for fixed l and store them in an array
    //   apply blocking method to this array; find true std automatically

    for (int l = 0; l < nac; l++) {
      int j = linval[l];
      cout << "   " << l << " " << nac << endl;
      for (int i = 0; i < nsteps-j; i++) {
        txi = tx[i];
	txj = tx[j+i];
	tyi = ty[i];
	tyj = ty[j+i];
	blockdata[i] = txi*txj + tyi*tyj;	
      }
      // compute mean and true standard deviation from the blocking method
      blocking_method(blockdata, block_std, block_uncert, nsteps-j, ccl, stdl);
      ac[l] = ccl;
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

  void compute_ttac_with_errorbars(Performance_tools* performance_tools,
				   int* t_ac, int* time, double* tx, double* ty, double* ac, double* std,
		 		   double* blockdata, double* block_std, double* block_uncert, int* linval,
		 		   int nsteps, int nac, int limit, int dt) {
    performance_tools->compute_ttac_with_errorbars(t_ac, time, tx, ty, ac, std, blockdata, block_std, block_uncert, linval, nsteps, nac, limit, dt);
  }
  
}

