using Pkg; Pkg.activate("analysis_code/I2julia")
using I2julia

hdf5dir  = "/home/zierler_fabian/Nextcloud/isospin2/HDF5_files/" 
hdf5list = readdir(hdf5dir,join=true)
filter!(isfile,hdf5list)
write_ensemble_list(hdf5list;filename="ensembles.csv")
