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
}

