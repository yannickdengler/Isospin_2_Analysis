using Pkg; Pkg.activate("Isospin_2_Analysis/I2julia")
using I2julia

hdf5dir  = "HDF5_logfiles/" 
hdf5list = readdir(hdf5dir,join=true)
filter!(isfile,hdf5list)
write_averaged_hdf5_files(hdf5list)
