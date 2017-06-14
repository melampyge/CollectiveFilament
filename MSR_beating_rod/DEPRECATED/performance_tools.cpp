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

  void compute_msr(int* t_msr, int* time, double* phi,
		   double* msr, int* counter, int* linval,
		   int nsteps, int nmsr, int limit, int dt) {
    // allocate variables
    double pi, pj;
    
    // zero msd and counter
    for (int i = 0; i < nmsr; i++) msr[i] = 0.0;
    for (int i = 0; i < nmsr; i++) counter[i] = 0;

    // fill linval array
    for (int i = 0; i < nmsr; i++) linval[i] = i*dt;

    // fill t_msr array
    for (int i = 0; i < nmsr; i++) t_msr[i] = time[linval[i]] - time[0];
  
    // compute the MSD
    for (int i = 0; i < nsteps; i++) {
      pi = phi[i];
      for (int l = 0; l < nmsr; l++) {
        int j = linval[l];
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

  void compute_msr_with_derivatives(int* t_msr, int* time, double* phi, double* dphi, double* d2phi,
				    double* msr, double* dmsr, double* d2msr, double* d2msr_t1, double* d2msr_t2, int* counter, int* linval,
				    int nsteps, int nmsr, int limit, int dt) {
    // allocate variables
    double pi, pj;
    double dp, d2p;
    
    // zero msd and counter
    for (int i = 0; i < nmsr; i++) msr[i] = 0.0;
    for (int i = 0; i < nmsr; i++) dmsr[i] = 0.0;
    for (int i = 0; i < nmsr; i++) d2msr[i] = 0.0;
    for (int i = 0; i < nmsr; i++) d2msr_t1[i] = 0.0;
    for (int i = 0; i < nmsr; i++) d2msr_t2[i] = 0.0;

    for (int i = 0; i < nmsr; i++) counter[i] = 0;

    // fill linval array
    for (int i = 0; i < nmsr; i++) linval[i] = i*dt;

    // fill t_msr array
    for (int i = 0; i < nmsr; i++) t_msr[i] = time[linval[i]] - time[0];
  
    // compute the MSR and derivatives
    for (int i = 0; i < nsteps; i++) {
      pi = phi[i];
      for (int l = 0; l < nmsr; l++) {
        int j = linval[l];
        if (j + i < nsteps) {
          pj = phi[j+i];
	  dp = dphi[j+i];
	  d2p = d2phi[j+i];
          counter[l] += 1;
          msr[l] += pow(pi-pj,2);
	  dmsr[l] += 2*(pj-pi)*dp;
	  d2msr[l] += 2*(pj-pi)*d2p + 2*dp*dp;
	  d2msr_t1[l] += 2*(pj-pi)*d2p;
	  d2msr_t2[l] += 2*dp*dp;
        }
      }
    }
    
    // normalize MSR and derivatives
    for (int i = 0; i <nmsr; i++) msr[i] /= counter[i];
    for (int i = 0; i <nmsr; i++) dmsr[i] /= counter[i];
    for (int i = 0; i <nmsr; i++) d2msr[i] /= counter[i];
    for (int i = 0; i <nmsr; i++) d2msr_t1[i] /= counter[i];
    for (int i = 0; i <nmsr; i++) d2msr_t2[i] /= counter[i];

    // return
    return;
  }

  /***********************************************************************************/

  void compute_autocorrelation(int* t_ac, int* time, double* x, double* ac,
			       int* counter, int* linval,
			       int nsteps, int nac, int limit, int dt) {
    // allocate variables
    double xi, xj;
    
    // zero ac and counter
    for (int i = 0; i < nac; i++) ac[i] = 0.0;
    for (int i = 0; i < nac; i++) counter[i] = 0;

    // fill linval array
    for (int i = 0; i < nac; i++) linval[i] = i*dt;

    // fill t_msr array
    for (int i = 0; i < nac; i++) t_ac[i] = time[linval[i]] - time[0];
  
    // compute the autocorrelation
    for (int i = 0; i < nsteps; i++) {
      xi = x[i];
      for (int l = 0; l < nac; l++) {
        int j = linval[l];
        if (j + i < nsteps) {
          xj = x[j+i];
          counter[l] += 1;
          ac[l] += xi*xj;
        }
      }
    }
    
    // normalize autocorrelation
    for (int i = 0; i <nac; i++) ac[i] /= counter[i];

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

  void compute_autocorrelation_with_errorbars(int* t_ac, int* time, double* x, double* ac,
					      double* std, double* blockdata, double* block_std,
					      double* block_uncert, int* linval,
					      int nsteps, int nac, int limit, int dt) {
    // allocate variables
    double xi, xj;
    double acl, stdl;
    
    // zero ac, std and counter
    for (int i = 0; i < nac; i++) ac[i] = 0.0;
    for (int i = 0; i < nac; i++) std[i] = 0.0;

    // fill linval array
    for (int i = 0; i < nac; i++) linval[i] = i*dt;

    // fill t_msr array
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
        xi = x[i];
	xj = x[j+i];
	blockdata[i] = xi*xj;	
      }
      // compute mean and true standard deviation from the blocking method
      blocking_method(blockdata, block_std, block_uncert, nsteps-j, acl, stdl);
      ac[l] = acl;
      std[l] = stdl;
    }

    // return
    return;
  }

  /***********************************************************************************/

  void compute_msr_with_errorbars(int* t_msr, int* time, double* phi,
				  double* msr, double* std,
				  double* blockdata, double* block_std, double* block_uncert,
				  int* linval,
				  int nsteps, int nmsr, int limit, int dt) {
    // allocate variables
    double pi, pj;
    double msrl, stdl;
    
    // zero msd and counter
    for (int i = 0; i < nmsr; i++) msr[i] = 0.0;
    for (int i = 0; i < nmsr; i++) std[i] = 0.0;

    // fill linval array
    for (int i = 0; i < nmsr; i++) linval[i] = i*dt;

    // fill t_msr array
    for (int i = 0; i < nmsr; i++) t_msr[i] = time[linval[i]] - time[0];
  
    // compute the MSR with statistical uncertainties
    //   using the blocking method:
    //   have l as the outer loop
    //   compute all values for fixed l and store them in an array
    //   apply blocking method to this array; find true std automatically


    for (int l = 0; l < nmsr; l++) {
      int j = linval[l];
      cout << "   " << l << " " << nmsr << endl;
      for (int i = 0; i < nsteps - j; i++) {
        pi = phi[i];
	pj = phi[i+j];
	blockdata[i] = pow(pi-pj,2);
      }
      // compute mean and true standard deviation from the blocking method
      blocking_method(blockdata, block_std, block_uncert, nsteps-j, msrl, stdl);
      msr[l] = msrl;
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

  void compute_angle(Performance_tools* performance_tools,
		     double* ex, double* ey, double* phi, double* dp, int nsteps) {
    performance_tools->compute_angle(ex, ey, phi, dp, nsteps);
  }
  
  /***********************************************************************************/

  void compute_msr(Performance_tools* performance_tools,
		   int* t_msr, int* time,
		   double* phi, double* msr, int* counter, int* linval,
		   int nsteps, int nmsr, int limit, int dt) {
    performance_tools->compute_msr(t_msr, time, phi, msr, counter, linval, nsteps, nmsr, limit, dt);
  }

  /***********************************************************************************/

  void compute_msr_with_derivatives(Performance_tools* performance_tools,
				    int* t_msr, int* time,
				    double* phi, double* dphi, double* d2phi,
				    double* msr, double* dmsr, double* d2msr,
				    double* d2msr_t1, double* d2msr_t2,
				    int* counter, int* linval,
				    int nsteps, int nmsr, int limit, int dt) {
    performance_tools->compute_msr_with_derivatives(t_msr, time, phi, dphi, d2phi, msr, dmsr, d2msr, d2msr_t1, d2msr_t2, counter, linval, nsteps, nmsr, limit, dt);
  }

  /***********************************************************************************/

  void compute_autocorrelation(Performance_tools* performance_tools,
			       int* t_ac, int* time, double* x, double* ac,
       			       int* counter, int* linval,
			       int nsteps, int nac, int limit, int dt) {
    performance_tools->compute_autocorrelation(t_ac, time, x, ac, counter, linval, nsteps, nac, limit, dt);
  }

  /***********************************************************************************/

  void compute_autocorrelation_with_errorbars(Performance_tools* performance_tools,
					      int* t_ac, int* time, double* x, double* ac, double* std,
					      double* blockdata, double* block_std, double* block_uncert, int* linval,
					      int nsteps, int nac, int limit, int dt) {
    performance_tools->compute_autocorrelation_with_errorbars(t_ac, time, x, ac, std, blockdata, block_std, block_uncert, linval, nsteps, nac, limit, dt);
  }

  /***********************************************************************************/

  void compute_msr_with_errorbars(Performance_tools* performance_tools,
				  int* t_msr, int* time,
				  double* phi, double* msr, double* std,
				  double* blockdata, double* block_std, double* block_uncert,
				  int* linval, int nsteps, int nmsr, int limit, int dt) {
    performance_tools->compute_msr_with_errorbars(t_msr, time, phi, msr, std, blockdata, block_std, block_uncert, linval, nsteps, nmsr, limit, dt);
  }
  
}

