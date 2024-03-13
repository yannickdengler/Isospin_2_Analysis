using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia
using HDF5
using Plots

h5corrs = "output/correlators.hdf5"
h5fits  = "output/fitresults.hdf5"

corr_id = h5open(h5corrs,"r")
fits_id = h5open(h5fits,"r")
ensembles = keys(corr_id)
ensemble  = ensembles[9]

#=
ππo = 49
T   = read(corr_id,joinpath(ensemble,"N_T"))
cππ = read(corr_id,joinpath(ensemble,"correlator_deriv"))[:,:,ππo]
=#

ππo = 45
T   = read(corr_id,joinpath(ensemble,"N_T"))
cππ = read(corr_id,joinpath(ensemble,"correlator"))[:,:,ππo]

meff, Δmeff = implicit_meff_jackknife(cππ;sign=+1)

t = 2:T
scatter(t,meff[t],yerr=Δmeff[t])
