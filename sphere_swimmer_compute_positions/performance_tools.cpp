#include <iostream>
#include <math.h>

using namespace std;

class Performance_tools{
  public:

  /***********************************************************************************/

  double neigh_min(double dx, double lx) {
    double dx1, dx2; 
    dx1 = dx + lx;
    dx2 = dx - lx;
    if (dx*dx < dx1*dx1 && dx*dx < dx2*dx2) return dx;
    if (dx1*dx1 < dx2*dx2) return dx1;
    return dx2;
  }

  /***********************************************************************************/

  void correct_pbc(double* comx, double* dx, double lx, int nsteps, int nmol) {
    // declare variables
    double dxj;
    // loop over all molecules
    for (int i = 0; i < nmol; i++) {
      // loop over all steps except the first
      for (int j = 1; j < nsteps; j++) {
	dxj = comx[i*nsteps + j] - comx[i*nsteps + j - 1];
	dx[j] = neigh_min(dxj,lx);
      }
      for (int j = 1; j < nsteps; j++) comx[i*nsteps + j] = comx[i*nsteps + j - 1] + dx[j];
    }
  }

  /***********************************************************************************/

  void compute_msd(int* t_msd, int* time, double* comx, double* comy,
		   double* msd, int* counter, int* logval,
		   int nsteps, int nmsd, int limit) {
    // allocate variables
    double xi, yi, xj, yj;
    double lmin, lmax, dl;
    
    // zero msd and counter
    for (int i = 0; i < nmsd; i++) msd[i] = 0.0;
    for (int i = 0; i < nmsd; i++) counter[i] = 0;

    // fill logval array
    lmin = log(1);
    lmax = log(limit-1);
    dl = (lmax - lmin)/(nmsd-1);
    for (int i = 0; i < nmsd; i++) logval[i] = (int) exp(lmin + i*dl);

    // fill t_msd array
    for (int i = 0; i < nmsd; i++) t_msd[i] = time[logval[i]] - time[0];

    // compute the MSD
    for (int i = 0; i < nsteps; i++) {
      xi = comx[i];
      yi = comy[i];
      for (int l = 0; l < nmsd; l++) {
        int j = logval[l];
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
    
  }

  /***********************************************************************************/

  void compute_msd_lin(int* t_msd, int* time, double* comx, double* comy,
		   double* msd, int* counter, int* linval,
		   int nsteps, int nmsd, int limit) {
    // allocate variables
    double xi, yi, xj, yj;
    
    // zero msd and counter
    for (int i = 0; i < nmsd; i++) msd[i] = 0.0;
    for (int i = 0; i < nmsd; i++) counter[i] = 0;

    // fill linval array
    for (int i = 0; i < nmsd; i++) linval[i] = i;

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
    
  }
  
};

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/


extern "C" {
  Performance_tools* Performance_tools_new(){return new Performance_tools();}

  void correct_pbc(Performance_tools* performance_tools,
		   double* comx, double* dx, double lx, int nsteps, int nmol) {

    performance_tools->correct_pbc(comx, dx, lx, nsteps, nmol);
  }

  /***********************************************************************************/

  void compute_msd(Performance_tools* performance_tools,
		   int* t_msd, int* time,
		   double* comx, double* comy, double* msd_tmp, int* counter, int* logval,
		   int nsteps, int nmsd, int limit) {
    performance_tools->compute_msd(t_msd, time, comx, comy, msd_tmp, counter, logval, nsteps, nmsd, limit);
  }

  /***********************************************************************************/

  void compute_msd_lin(Performance_tools* performance_tools,
		   int* t_msd, int* time,
		   double* comx, double* comy, double* msd_tmp, int* counter, int* linval,
		   int nsteps, int nmsd, int limit) {
    performance_tools->compute_msd_lin(t_msd, time, comx, comy, msd_tmp, counter, linval, nsteps, nmsd, limit);
  }
  
}

