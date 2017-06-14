
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

void compute_static_structure (double *S, double *x, double *y, int natoms, int Nbins, int Nmax, double bin_width, double delk) {
  /* calculate static structure per timestep with a running average */

  for (int nx = 0; nx < Nmax; nx++) {
    double kx = delk*nx;

      for (int ny = 0; ny < Nmax; ny++) {
	double ky = delk*ny;
	double k = sqrt(kx*kx + ky*ky);
	int seg = static_cast<int>(k/bin_width);
	if (seg >= Nbins) seg = Nbins-1; 
	double costerm = 0.;
	double sinterm = 0.;

	omp_set_num_threads(8);
	#pragma omp parallel for reduction(+:costerm,sinterm)
	for (int j = 0; j < natoms; j++) {
	  double dotp = kx*x[j] + ky*y[j];
	  costerm += cos(dotp);
	  sinterm += sin(dotp);
  
	} // particle loop

	S[seg] += costerm*costerm + sinterm*sinterm;

      } // ky loop
  } // kx loop        
  
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
  double lamda_min = 3.325; 
  double lamda_max = l/2.;
  double kmax = 2*pi/lamda_min;
  double kmin = 2*pi/lamda_max;
  
  cout << "delta k = " << delk << " kmin = " << kmin << " kmax = " << kmax << endl;
  
  int Nbins = 100;
  int njump = 20;
  int Nmax = 300;
  double bin_width = (kmax-kmin)/Nbins;
  
  cout << "Nbins = " << Nbins << " bin_width = " << bin_width << endl;
  cout << "Nmax = " << Nmax << endl;
    
  double *S = new double[Nbins];
  for (int j = 0; j < Nbins; j++) S[j] = 0.;
    
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
    
    compute_static_structure(S, xt, yt, natoms, Nbins, Nmax, bin_width, delk);
    
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
  for (int j = 0; j < Nbins; j++)  { 
    S[j] /= (natoms*totalData); 
    outfile << j*bin_width << "\t\t" << S[j] << "\n";
  }

  outfile.close();
  
  return 0;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
