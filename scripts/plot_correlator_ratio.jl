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

BC = correlators[:,:,43]
AD = correlators[:,:,41]
corr_pipi = AD - BC

using Statistics
R = @. corr_pipi / corr_pi^2
T,N = size(R)
R_mean = dropdims(mean(R,dims=2),dims=2)
R_sdev = dropdims(std(R,dims=2),dims=2)/sqrt(N)

R_deriv = zeros(T-1,N)
for t in 1:T-1, i in 1:N
    R_deriv[t,i] = (corr_pipi[t,i] - corr_pipi[t+1,i]) / (corr_pi[t,i]^2 - corr_pi[t+1,i]^2)
end
R_deriv_mean = dropdims(mean(R_deriv,dims=2),dims=2)
R_deriv_sdev = dropdims(std(R_deriv,dims=2),dims=2)/sqrt(N)

using Plots
scatter(R_mean,yerr=R_sdev)
scatter(R_deriv_mean,yerr=R_deriv_sdev)