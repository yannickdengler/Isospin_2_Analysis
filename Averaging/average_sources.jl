using Pkg; Pkg.activate("Isospin_2_Analysis/I2julia")
using I2julia
using HDF5

hdf5dir  = "HDF5_logfiles/" 
hdf5list = readdir(hdf5dir,join=true)
filter!(isfile,hdf5list)

file_id = h5open(hdf5list[1])
corr = read(file_id,"correlators")
ops = read(file_id,"operators")
ops = read(file_id,"plaquette")

keys_hdf = keys(file_id)
keys_without_correlators = filter(!isequal("correlators"),keys(file_id))

# average over all sources and perform numerical derivative
source_averaged, nsrc = average_sources(hdf5list[1]) 
source_averaged_deriv = correlator_derivative(source_averaged)

# average over MC samples
corr, Δcorr, ncfg = average_configurations(source_averaged)
corr_deriv, Δcorr_deriv, ncfg = average_configurations(source_averaged_deriv)

source_averaged[:,:,1]
meff, Δmeff = implicit_meff_jackknife(source_averaged[:,:,1];sign=+1)
meff
Δmeff
using Plots
scatter(meff,yerr=Δmeff)

#=
# Plot a test
using Plots
function plot_corr(c,Δc,ops;op_ind=1)
    T, nops = size(c)
    @assert length(ops) == nops
    return scatter(c[:,op_ind],yerr=Δc[:,op_ind],label=ops[op_ind],yscale=:log10,legend=:bottomright)
end
plt1 = plot_corr(corr,Δcorr,ops,op_ind=1)
plt2 = plot_corr(abs.(corr_deriv),Δcorr_deriv,ops,op_ind=1)
=#
