using Pkg; Pkg.activate("analysis_code/I2julia")
using I2julia

hdf5dir  = "HDF5_logfiles/" 
hdf5list = readdir(hdf5dir,join=true)
filter!(isfile,hdf5list)
write_ensemble_list(hdf5list;filename="ensembles.csv")
