
import h5py
import numpy as np
import os.path


def main():
    
    path = '/local/duman/SIMULATIONS/many_polymers_5/density_0.08/kappa_2.5/fp_0.24/CLUSTER/cluster_sizes_23500000.txt'
    
    print path
    if os.path.exists(path):
        txtfile = np.loadtxt(path, dtype=int)
        print np.shape(txtfile)
        
#        f = h5py.File('test.hdf5', 'w')
#        size_grp = f.create_group('size')
#        size_grp.create_dataset('size', data=txtfile, compression='gzip')
        
        f = h5py.File('test.hdf5', 'r')
        size_grp = f['size']
        size = size_grp['size']
        print np.asarray(size)
        print size.shape
        
#        for el in data:
#            print el

        



if __name__ == '__main__':
    main()