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
		   double* msr, int* counter, int* logval,
		   int nsteps, int nmsr, int limit, int nevery) {
    // allocate variables
    double pi, pj;
    double lmin, lmax, dl;
    
    // zero msd and counter
    for (int i = 0; i < nmsr; i++) msr[i] = 0.0;
    for (int i = 0; i < nmsr; i++) counter[i] = 0;

    // fill logval array
    lmin = log(1);
    lmax = log(limit-1);
    dl = (lmax - lmin)/(nmsr-1);
    for (int i = 0; i < nmsr; i++) logval[i] = (int) exp(lmin + i*dl);

    // fill t_msr array
    for (int i = 0; i < nmsr; i++) t_msr[i] = time[logval[i]] - time[0];

    // compute angle
    
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
		   double* phi, double* msr_tmp, int* counter, int* logval,
		   int nsteps, int nmsr, int limit, int nevery) {
    performance_tools->compute_msr(t_msr, time, phi, msr_tmp, counter, logval, nsteps, nmsr, limit, nevery);
  }
}

