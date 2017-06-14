#include <iostream>
#include <math.h>

using namespace std;

class Performance_tools{
  public:

  /***********************************************************************************/

  void compute_msd(int* t_msd, int* time, double* comx, double* comy,
		   double* msd, int* counter, int* linval,
		   int nsteps, int nmsd, int limit, int dt) {
    // allocate variables
    double xi,xj,yi,yj;
    
    // zero msd and counter
    for (int i = 0; i < nmsd; i++) msd[i] = 0.0;
    for (int i = 0; i < nmsd; i++) counter[i] = 0;

    // fill linval array
    for (int i = 0; i < nmsd; i++) linval[i] = i*dt;

    // fill t_msd array
    for (int i = 0; i < nmsd; i++) t_msd[i] = time[linval[i]] - time[0];
  
    // compute the MSD
    for (int i = 0; i < nsteps; i++) {
      xi = comx[i];
      yi = comy[i];
      for (int l = 0; l < nmsd; l++) {
        int j = linval[l];
        if (j + i < nsteps) {
          xj = comx[j+i];
	  yj = comy[j+i];
          counter[l] += 1;
          msd[l] += pow(xi-xj,2) + pow(yi-yj,2);
        }
      }
    }
    
    // normalize MSD
    for (int i = 0; i <nmsd; i++) msd[i] /= counter[i];

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

  void compute_msd_with_errorbars(int* t_msd, int* time, double* comx, double* comy,
				  double* msd, double* std, double* blockdata, double* block_std,
				  double* block_uncert, int* linval,
				  int nsteps, int nmsd, int limit, int dt) {
    // allocate variables
    double xi,xj,yi,yj;
    double msdl, stdl;
    
    // zero msd and std
    for (int i = 0; i < nmsd; i++) msd[i] = 0.0;
    for (int i = 0; i < nmsd; i++) std[i] = 0.0;

    // fill linval array
    for (int i = 0; i < nmsd; i++) linval[i] = i*dt;

    // fill t_msd array
    for (int i = 0; i < nmsd; i++) t_msd[i] = time[linval[i]] - time[0];
  
    // compute the MSD with statistical uncertainties
    //   using the blocking method:
    //   have l as the outer loop
    //   compute all values for fixed l and store them in an array
    //   apply blocking method to this array; find true std automatically

    for (int l = 0; l < nmsd; l++) {
      int j = linval[l];
      cout << "   " << l << " " << nmsd << endl;
      for (int i = 0; i < nsteps - j; i++) {
        xi = comx[i];
	yi = comy[i];
	xj = comx[i+j];
	yj = comy[i+j];
	blockdata[i] = pow(xi-xj,2) + pow(yi-yj,2);
      }
      // compute mean and true standard deviation from the blocking method
      blocking_method(blockdata, block_std, block_uncert, nsteps-j, msdl, stdl);
      msd[l] = msdl;
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

  void compute_msd(Performance_tools* performance_tools,
		   int* t_msd, int* time,
		   double* comx, double* comy, double* msd, int* counter, int* linval,
		   int nsteps, int nmsd, int limit, int dt) {
    performance_tools->compute_msd(t_msd, time, comx, comy, msd, counter, linval, nsteps, nmsd, limit, dt);
  }

  /***********************************************************************************/

  void compute_msd_with_errorbars(Performance_tools* performance_tools,
				  int* t_msd, int* time,
				  double* comx, double* comy, double* msd, double* std,
				  double* blockdata, double* block_std,
				  double* block_uncert, int* linval,
				  int nsteps, int nmsd, int limit, int dt) {
    performance_tools->compute_msd_with_errorbars(t_msd, time, comx, comy, msd, std, blockdata, block_std, block_uncert, linval, nsteps, nmsd, limit, dt);
  }

}

