using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia

hdf5dir  = "./output/HDF5_logfiles/" 
hdf5list = readdir(hdf5dir,join=true)
filter!(isfile,hdf5list)
write_ensemble_list(hdf5list;filename="output/ensembles.csv")
