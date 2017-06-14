#include <iostream>
#include <math.h>

using namespace std;

class Performance_tools{
  public:

  void fill_neigh_matrix(int* neighs_mol, int* llist, int* head,
			 int nsegx, int nsegy, double* x, double* y,
			 double* tx, double* ty, int* mol, int nmol, double lx, double ly,
			 double dcrit, double ccrit) {
   
    // declare variable
    int sv1, sv2, i, j, val1, val2,i2,j2,a,b, mol1, mol2;
    double tx1, tx2, ty1, ty2, dx, dy, rsq, x1, x2, y1, y2;


   
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
	      tx1 = tx[val1];
	      ty1 = ty[val1];
	      mol1 = mol[val1];
              while (val2 != 0) {
                if (mol[val1] != mol[val2]) {
		  x2 = x[val2]/lx;
                  y2 = y[val2]/ly;
		  tx2 = tx[val2];
		  ty2 = ty[val2];
		  mol2 = mol[val2];

		  if (mol1 != mol2) {

                    dx = x2-x1;
                    dx = dx - floor(dx + 0.5);
                    dx = dx*lx;
                    dy = y2-y1;
                    dy = dy - floor(dy + 0.5);
                    dy = dy*ly;
                    rsq = dx*dx + dy*dy;
                    if (rsq < dcrit*dcrit) {
		      if (tx1*tx2 + ty1*ty2 >= ccrit) {
		        neighs_mol[mol1*nmol + mol2] += 1;
			neighs_mol[mol2*nmol + mol1] += 1;
		      }
		    }
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

  /*****************************************************************************************/

  void recursion(int* neighs, int* cl, int i, int nmol, double ncrit) {
    for (int j = 0; j < nmol; j++) {
      if (cl[j] == -1) {
        if (neighs[i*nmol + j] > ncrit) {
	  cl[j] = cl[i];
	  recursion(neighs,cl,j,nmol,ncrit);
	}
      }
    }
  }

  void cluster_search(int* neighs, int* cl, double lcrit, int nfil, int nmol) {
    // define critical overlap
    double ncrit = lcrit;
    //double ncrit = lcrit*nfil;
    int clmax = 0;
    // loop over all columns of the matrix:
    for (int i = 0; i < nmol; i++) {
      // assign a cluster id to the current molecule
      if (cl[i] == -1) {
        cl[i] = clmax;
	clmax++;
        recursion(neighs,cl,i,nmol,ncrit);
      }
    }
  }

  /*****************************************************************************************/

  void recurse_cluster(int* h, int* h2, int nx, int ny, int sx0, int sy0) {
    int sx;
    int sy;
    for (int i = -1; i < 2; i++) {
      for (int j = -1; j < 2; j++) {
        sx = sx0 + i;
        sy = sy0 + j;
        if (sx < 0 || sx >= nx || sy < 0 || sy >= ny) continue;
        if (h[sx*ny+sy] == 1 && h2[sx*ny+sy] == 0) {
          h2[sx*ny+sy] = 1;
          recurse_cluster(h,h2,nx,ny,sx,sy);
        }
      }
    }
  }
  
};

/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/


extern "C" {
  Performance_tools* Performance_tools_new(){return new Performance_tools();}

  void fill_neigh_matrix(Performance_tools* performance_tools,
			 int* neighs_mol, int* llist, int* head,
			 int nsegx, int nsegy, double* x, double* y,
			 double* tx, double* ty, int* mol, int nmol, double lx, double ly,
			 double dcrit, double ccrit) {

    performance_tools->fill_neigh_matrix(neighs_mol, llist, head,
					 nsegx, nsegy, x, y,
					 tx, ty, mol, nmol, lx, ly,
					 dcrit, ccrit);
  }

  /*****************************************************************************************/

  void cluster_search(Performance_tools* performance_tools,
		      int* neighs_mol, int* cl, double lcrit, int nfil, int nmol) {
    performance_tools->cluster_search(neighs_mol, cl, lcrit, nfil, nmol);
  }

    /*****************************************************************************************/

  void recurse_cluster(Performance_tools* performance_tools,
		       int* h, int* h2, int nx, int ny, int sx0, int sy0) {
    performance_tools->recurse_cluster(h, h2, nx, ny, sx0, sy0);
  }
  
}

