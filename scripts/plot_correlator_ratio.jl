using Pkg; Pkg.activate("./scripts/src_jl")
using Plots; plotlyjs(ms=7,frame=:box)
using HDF5
using I2julia

hdf5dir  = "./output/HDF5_source_average"
hdf5list = readdir(hdf5dir,join=true) 
hdf5file = hdf5list[4]

fid = h5open(hdf5file,"r")

correlators = read(fid,"correlator")
corr_pipi = correlators[:,:,49]
corr_pi   = correlators[:,:,45]
@assert read(fid,"operators")[45] == "pi"
@assert read(fid,"operators")[49] == "pipi"

using Statistics
R = @. corr_pipi / corr_pi^2
N = size(R)[2]
R_mean = dropdims(mean(R,dims=2),dims=2)
R_sdev = dropdims(std(R,dims=2),dims=2)/sqrt(N)

using Plots
scatter(R_mean,yerr=R_sdev)