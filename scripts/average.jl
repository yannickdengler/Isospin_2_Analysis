using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia

hdf5dir  = "./output/HDF5_logfiles/" 
hdf5list = readdir(hdf5dir,join=true)
filter!(isfile,hdf5list)
write_averaged_hdf5_files(hdf5list)
