#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
using namespace std; 

/*-----------------------------------------------------------------
Compute the order parameter based on preparations made by the
python script
-----------------------------------------------------------------*/
int main()
{

  // declare variables
  int nall;
  double *x, *y, *cphi, *sphi;
  int *llist, *head;
  double lx,ly,rc;
  int nsegx, nsegy, ncells;
  int segx, segy, cell;
  int i,j,i2,j2,a,b;
  int sv1, sv2, val1, val2;
  double cos1, cos2, sin1, sin2;
  double x1,x2,y1,y2,dx,dy,dist;
  double Si = 0;
  int counter = 0;
  
  // read in values
  ifstream infile;
  infile.open("tmp.data");
  if (!infile.is_open()) {
    cout << "Error, could not open file \n";
  }
  // read in first line of variables
  infile >> nall >> rc >> lx >> ly;
  // allocate arrays
  x = (double*) malloc(nall*sizeof(double));
  y = (double*) malloc(nall*sizeof(double));
  cphi = (double*) malloc(nall*sizeof(double));
  sphi = (double*) malloc(nall*sizeof(double));
  // fill arrays
  for (i = 0; i < nall; i++)
    infile >> x[i] >> y[i] >> cphi[i] >> sphi[i];
  // close the file
  infile.close();

  // compute value for linked list
  nsegx = (int) lx/rc;
  nsegy = (int) ly/rc;
  ncells = nsegx*nsegy;
  // allocate memory for linked list
  llist = (int*) malloc(nall*sizeof(int));
  head =  (int*) malloc(ncells*sizeof(int));
  // set linked list to zero
  for (i = 0; i < nall; i++) llist[i] = 0;
  for (i = 0; i < ncells; i++) head[i] = 0;
  
  // generate the linked list
  for (i = 0; i < nall; i++) {
    segx = (int) x[i]/lx*nsegx;
    segy = (int) y[i]/ly*nsegy;
    cell = segx*nsegy + segy;
    llist[i] = head[cell];
    head[cell] = i;
  }


  
  // loop over the linked list
  for (i = 0; i < nsegx; i++) {
    for (j = 0; j < nsegy; j++) {
      sv1 = head[i*nsegy + j];
      for (a = - 1; a < 2; a++) {
	i2 = (i+a+nsegx)%nsegx;
	for (b = -1; b < 2; b++) {
	  j2 = (j+b+nsegy)%nsegy;
	  sv2 = head[i2*nsegy + j2];

	  // restore header values
	  val1 = sv1;
	  val2 = sv2;
	  while (val1 != 0) {
            cos1 = cphi[val1];
	    sin1 = sphi[val1];
	    x1 = x[val1]/lx;
	    y1 = y[val1]/ly;
            while (val2 != 0) {

	      if (val1 != val2) {
                x2 = x[val2]/lx;
		y2 = y[val2]/ly;
	        dx = x1-x2;
	        dx = dx - floor(dx + 0.5);
	        dx = dx*lx;
	        dy = y1-y2;
	        dy = dy - floor(dy + 0.5);
	        dy = dy*ly;
	        dist = dx*dx + dy*dy;
	        if (dist <= rc*rc) {
                  cos2 = cphi[val2];
	          sin2 = sphi[val2];
	          Si += cos1*cos2 + sin1*sin2;
	          counter += 1;
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

  Si = Si/counter;
  // write result to file
  ofstream ofile;
  ofile.open("tmp2.data");
  if (!ofile.is_open()) {
    cout << "Error, could not open output file \n";
  }
  ofile << Si << endl;
  ofile.close();
  
  // free arrays
  free(x);
  free(y);
  free(cphi);
  free(sphi);
  free(llist);
  free(head);

  return 0;
}
