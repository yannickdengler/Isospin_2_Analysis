using Pkg; Pkg.activate("analysis_code/I2julia")
using I2julia
using HDF5

hdf5dir  = "HDF5_logfiles/" 
hdf5list = readdir(hdf5dir,join=true)
filter!(isfile,hdf5list)

file_id = h5open(hdf5list[1])
corr = read(file_id,"correlators")

source_averaged, nsrc = average_sources(hdf5list[1]) 
avg_corr, Δavg_corr, ncfg = average_configurations(hdf5list[1])

avg_corr
Δavg_corr