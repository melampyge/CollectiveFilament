#include <iostream>
#include <math.h>

using namespace std;

class Performance_tools{
  public:

  /***********************************************************************************/

  void compute_angle(double* ex, double* ey, double* phi, double* dp, int nsteps) {
    // allocate variables
    double dpi;
    const double pi = atan2(0,-1);
    // compute phi
    for (int i = 0; i < nsteps; i++) phi[i] = atan2(ey[i], ex[i]);

    // compute dp
    for (int i = 1; i < nsteps; i++) {
      dpi = phi[i] - phi[i-1];
      if (dpi < -pi) dpi += 2*pi;
      if (dpi > pi) dpi -= 2*pi;
      dp[i] = dpi;
    }

    // correct phi
    for (int i = 1; i < nsteps; i++) phi[i] = phi[i-1] + dp[i];

    // return
    return;
  }

  /***********************************************************************************/

  void compute_msr_loglog(int* t_msr, int* time, double* phi,
		   double* msr, int* counter, int* logval,
		   int nsteps, int nmsr, int limit) {
    // allocate variables
    double pi, pj;
    double lmin, lmax, dl;
    
    
    // zero msd and counter
    for (int i = 0; i < nmsr; i++) msr[i] = 0.0;
    for (int i = 0; i < nmsr; i++) counter[i] = 0;

    // fill logval array
    lmax = limit;
    lmin = log(1.0);
    lmax = log(lmax);
    dl = (lmax - lmin) / nmsr;
    for (int i = 0; i < nmsr; i++) logval[i] = int (exp(lmin + dl*i));

    // fill t_msr array
    for (int i = 0; i < nmsr; i++) t_msr[i] = time[logval[i]] - time[0];
  
    // compute the MSD
    for (int i = 0; i < nsteps; i++) {
      pi = phi[i];
      for (int l = 0; l < nmsr; l++) {
        int j = logval[l];
        if (j + i < nsteps) {
          pj = phi[j+i];
          counter[l] += 1;
          msr[l] += pow(pi-pj,2);
        }
      }
    }
    
    // normalize MSD
    for (int i = 0; i <nmsr; i++) msr[i] /= counter[i];

    // return
    return;
  }


  /***********************************************************************************/

  void compute_correlation(int* t_ee, int* time,
		           double* tx, double* ty, double* c_ee, int* counter, int* linval,
		           int nsteps, int limit) {
    // declare some variables
    double tx1, tx2, ty1, ty2;

    // zero c_ee and counter
    for (int i = 0; i < limit; i++) c_ee[i] = 0.0;
    for (int i = 0; i < limit; i++) counter[i] = 0;

    // fill linval array
    for (int i = 0; i < limit; i++) linval[i] = i;

    // fill t_ee array
    for (int i = 0; i < limit; i++) t_ee[i] = time[i];

    // compute the correlation function
    for (int i = 0; i < nsteps; i++) {
      tx1 = tx[i];
      ty1 = ty[i];
      for (int j = 0; j < limit; j++) {
        if (j + i < nsteps) {
          tx2 = tx[j+i];
	  ty2 = ty[j+i];
          counter[j] += 1;
          c_ee[j] += tx1*tx2 + ty1*ty2;
        }
      }
    }

    // normalize
    for (int i = 0; i < limit; i++) c_ee[i] /= counter[i];

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

  void compute_correlation_with_errorbars(int* t_ee, int* time,
		                          double* tx, double* ty, double* c_ee,
					  double* std, double* blockdata,
					  double* block_std, double* block_uncert,
					  int* linval, int nsteps, int limit) {
    // declare some variables
    double tx1, tx2, ty1, ty2;
    double ccj, stdj;

    // zero c_ee and std
    for (int i = 0; i < limit; i++) c_ee[i] = 0.0;
    for (int i = 0; i < limit; i++) std[i] = 0.0;

    // fill linval array
    for (int i = 0; i < limit; i++) linval[i] = i;

    // fill t_ee array
    for (int i = 0; i < limit; i++) t_ee[i] = time[i];

    // compute the correlation function with statistical uncertainties
    //   using the blocking methd:
    //   have j as the outer loop
    //   copute all values for fixed j and store them in an array
    //   apply blocking method to this array; find true std autmoatically

    for (int j = 0; j < limit; j++) {
      for (int i = 0; i < nsteps - j; i++) {
        tx1 = tx[i];
        ty1 = ty[i];
        tx2 = tx[i+j];
	ty2 = ty[i+j];
	blockdata[i] = tx1*tx2 + ty1*ty2;
        }
      // compute the mean and true standard deviation from the blocking method
      blocking_method(blockdata, block_std, block_uncert, nsteps-j, ccj, stdj);
      c_ee[j] = ccj;
      std[j] = stdj;
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

  void compute_angle(Performance_tools* performance_tools,
		     double* ex, double* ey, double* phi, double* dp, int nsteps) {
    performance_tools->compute_angle(ex, ey, phi, dp, nsteps);
  }

  /***********************************************************************************/

  void compute_msr_loglog(Performance_tools* performance_tools,
		   int* t_msr, int* time,
		   double* phi, double* msr, int* counter, int* logval,
		   int nsteps, int nmsr, int limit, int dt) {
    performance_tools->compute_msr_loglog(t_msr, time, phi, msr, counter, logval, nsteps, nmsr, limit);
  }

  /***********************************************************************************/

  void compute_correlation(Performance_tools* performance_tools,
		           int* t_ee, int* time,
		           double* tx, double* ty, double* c_ee, int* counter, int* linval,
		           int nsteps, int limit) {
    performance_tools->compute_correlation(t_ee, time, tx, ty, c_ee, counter, linval, nsteps, limit);
  }

  /***********************************************************************************/

  void compute_correlation_with_errorbars(Performance_tools* performance_tools,
		                          int* t_ee, int* time,
		                          double* tx, double* ty, double* c_ee,
					  double* std, double* blockdata,
					  double* block_std, double* block_uncert,
                                          int* linval, int nsteps, int limit) {
    performance_tools->compute_correlation_with_errorbars(t_ee, time, tx, ty, c_ee, std, blockdata, block_std, block_uncert, linval, nsteps, limit);
  }

}

