
//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cmath>
#include "/usr/users/iff_th2/duman/hdf5/include/hdf5.h"
//#include "/homec/jiff26/jiff2610/hdf/include/hdf5.h"
#include "omp.h"

#define pi M_PI

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_sim_data (char *filename, int &nsteps, int &natoms, int &nmol, double &l) {
  /* get basic simulation data in hdf5 format */
  
  // open the file pointer
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  // create buffers to get single data
  
  int i_buffer[1];
  i_buffer[0] = 0;
  double d_buffer[1];
  d_buffer[0] = 0.;
  
  // number of steps
  
  hid_t dataset = H5Dopen(file, "info/nsteps", H5P_DEFAULT);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, i_buffer);
  nsteps = i_buffer[0];
  
  // number of atoms
  
  dataset = H5Dopen(file, "info/natoms", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, i_buffer);
  natoms = i_buffer[0];
  
  // number of molecules

  dataset = H5Dopen(file, "info/nmol", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, i_buffer);
  nmol = i_buffer[0];
  
  // box size
  
  dataset = H5Dopen(file, "info/box/lx", H5P_DEFAULT);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, d_buffer);
  l = d_buffer[0];   

  H5Dclose(dataset);
  H5Fclose(file);
  
  return;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void compute_static_structure (double *S, double *x, double *y, double *kxs, double *kys, double *kvec, int natoms, int nks, int Nmax, double bin_width, double delk) {
  /* calculate static structure per timestep with a running average */

  for (int k = 0; k < Nmax; k++) {
    double term = 0.0;

      for (int l = 0; l < nks; l++) {
	double kx = kxs[l]*kvec[k];
	double ky = kys[l]*kvec[k];
	double costerm = 0.;
	double sinterm = 0.;

	omp_set_num_threads(4);
	#pragma omp parallel for reduction(+:costerm,sinterm)
	for (int j = 0; j < natoms; j++) {
	  double dotp = kx*x[j] + ky*y[j];
	  costerm += cos(dotp);
	  sinterm += sin(dotp);
  
	} // particle loop

	term += costerm*costerm + sinterm*sinterm;

      } // orientation loop
      
      term /= nks;
      S[k] += term;
      
  } // absolute value loop        
  
  return;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {
                            
  // get the file name by parsing
  
  char *filename = argv[1];
  cout << "Doing analysis on the following file: \n" << filename << endl;

  // get simulation data
  
  int nsteps = 0;
  int natoms = 0;
  int nmol = 0;
  double l = 0.;

  get_sim_data(filename, nsteps, natoms, nmol, l);

  // timesteps

  double time[nsteps];
  for (int i = 0; i < nsteps; i++) time[i] = 0.;
  
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  hid_t dataset = H5Dopen(file, "time", H5P_DEFAULT); 
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time);
  
  /* calculation of static structure factor 
  * load position data of particles at every timestep
  * calculate the static structure factor on a per-timestep basis
  * make a running average along timesteps 
  */
  
  // set preliminary data
  
  double delk = 2*pi/l;
  cout << "delta k = " << delk << endl;
  double lamda_min = 0.4; 
  double lamda_max = l;
  double kmax = 2*pi/lamda_min;
  double kmin = 2*pi/lamda_max;
  int nks = 24;
  
  cout << "delta k = " << delk << " kmin = " << kmin << " kmax = " << kmax << endl;
  
  int njump = 10;
  int Nmax = 2500;
  kmax = log(kmax);
  kmin = log(kmin);  
  double wbin = (kmax-kmin)/Nmax;

  // populate log-spaced absolute value of the k vectors
  
  double *kvec = new double[Nmax];
  for (int j = 0; j < Nmax; j++) kvec[j] = kmin + j*wbin;
  for (int j = 0; j < Nmax; j++) kvec[j] = exp(kvec[j]);
  
  // populate circular orientation of the k vector 
  
  double *kxs = new double[nks];
  double *kys = new double[nks];
  for (int j = 0; j < nks; j++) {
    kxs[j] = cos(2*pi*j/nks);
    kys[j] = sin(2*pi*j/nks);
  }
    
  cout << "Nmax = " << Nmax << " wbin = " << wbin << endl;
    
  double *S = new double[Nmax];
  for (int j = 0; j < Nmax; j++) S[j] = 0.;    
  
  for (int step = 0; step < nsteps; step += njump) {
    
    cout << "step / nsteps : " << step << " / " << nsteps << endl;
  
    // allocate arrays
    
    double *xt = new double[natoms];
    for (int i = 0; i < natoms; i++) xt[i] = 0.; 

    double *yt = new double[natoms];
    for (int i = 0; i < natoms; i++) yt[i] = 0.; 

    /* load particle position data at every timestep */
    
    // file dataspace
    
    hid_t x_dataset = H5Dopen(file, "pos/x", H5P_DEFAULT);
    hid_t y_dataset = H5Dopen(file, "pos/y", H5P_DEFAULT);

    hid_t x_dataspace = H5Dget_space(x_dataset);
    hid_t y_dataspace = H5Dget_space(y_dataset);
    
    // hyperslab in dataset
    
    hsize_t offset[2];
    offset[0] = step; offset[1] = 0;
    
    hsize_t count[2];
    count[0] = 1; count[1] = natoms;
    
    status = H5Sselect_hyperslab(x_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Sselect_hyperslab(y_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    // memory dataspace
    
    hid_t x_memspace, y_memspace;
    hsize_t dimsm[2];
    dimsm[0] = 1; dimsm[1] = natoms;
    
    x_memspace = H5Screate_simple(2, dimsm, NULL);
    y_memspace = H5Screate_simple(2, dimsm, NULL);
    
    // memory hyperslab
    
    hsize_t offset_out[2];
    offset_out[0] = 0; offset_out[1] = 0;
    
    hsize_t count_out[2];
    count_out[0] = 1; count_out[1] = natoms;
    
    status = H5Sselect_hyperslab(x_memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
    status = H5Sselect_hyperslab(y_memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
    
    // read data from hyperslab in the file into the hyperslab in memory and display
    
    status = H5Dread(x_dataset, H5T_NATIVE_DOUBLE, x_memspace, x_dataspace, H5P_DEFAULT, xt);  
    status = H5Dread(y_dataset, H5T_NATIVE_DOUBLE, y_memspace, y_dataspace, H5P_DEFAULT, yt);  
    
    // calculate static structure factor at this timestep and add it to the average 
    
    compute_static_structure(S, xt, yt, kxs, kys, kvec, natoms, nks, Nmax, wbin, delk);
    
    // deallocate arrays
    
    delete [] xt;
    delete [] yt;
    H5Sclose(x_dataspace);
    H5Dclose(x_dataset);
    H5Sclose(y_dataspace);
    H5Dclose(y_dataset);  
    
  }
    
  H5Dclose(dataset);
  H5Fclose(file);
  
  // normalize and write the values
  
  ofstream outfile("STRUCTURE/static_structure_factor_cpp.data");

  double totalData = nsteps/njump;
  for (int j = 0; j < Nmax; j++)  { 
    S[j] /= (natoms*totalData); 
    outfile << kvec[j] << "\t\t" << S[j] << "\n";
  }

  outfile.close();
  
  return 0;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
