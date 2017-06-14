#include <iostream>
#include <math.h>

using namespace std;

class Rdf_site{
  public:
  void compute(int nsegx, int nsegy, int natoms, int* head, int* llist,
	       int* mol, double* x, double* y, double* rdf,
	       double lx, double ly, double rmax, int nbin) {
   
    // declare variable
    int sv1, sv2, i, j, val1, val2,i2,j2,a,b,seg;
    double dx, dy, rsq, r, x1, x2, y1, y2;

    // loop over all central cells
    for (i = 0; i < nsegx; i++) {
      for (j = 0; j < nsegy; j++) {
        // store header of current cell
	sv1 = head[i*nsegy + j];
        //loop over neighboring cells
        for (a = 0; a < 3; a++) {
	  i2 = (i-1+a+nsegx)%nsegx;
          for (b = 0; b < 3; b++) {
            j2 = (j-1+b+nsegy)%nsegy;
            // store header of neighbor cell
            sv2 = head[i2*nsegy + j2];
                    
            // restore head values at for each new cell
            val1 = sv1;
            val2 = sv2;
	    while (val1 != 0) {
              x1 = x[val1]/lx;
              y1 = y[val1]/ly;
              while (val2 != 0) {
                if (mol[val1] != mol[val2]) {
		  x2 = x[val2]/lx;
                  y2 = y[val2]/ly;

                  dx = x2-x1;
                  dx = dx - floor(dx + 0.5);
                  dx = dx*lx;
                  dy = y2-y1;
                  dy = dy - floor(dy + 0.5);
                  dy = dy*ly;
                  rsq = dx*dx + dy*dy;
                  if (rsq < rmax*rmax) {
		     r = sqrt(rsq);
                     seg = int (r/rmax*nbin);
                     rdf[seg] += 1;
		  }
		}
                val2 = llist[val2];
	      }
              val1 = llist[val1];
	      val2 = sv2;
	    }
	  }
	}
      }
    }
  }

};

extern "C" {
  Rdf_site* Rdf_site_new(){return new Rdf_site();}

  void compute(Rdf_site* rdf_site, int nsegx, int nsegy,
	                int natoms, int* head, int* llist,
	                int* mol, double* x, double* y, double* rdf,
	                double lx, double ly, double rmax, int nbin){
    rdf_site->compute(nsegx, nsegy, natoms, head, llist,
	              mol, x, y, rdf, lx, ly, rmax, nbin);
  }
}

