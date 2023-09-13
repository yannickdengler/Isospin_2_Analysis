using Pkg; Pkg.activate("./scripts/src_jl")
using Plots; gr(ms=7)
using HDF5
using I2julia

h5dir = "./output/HDF5_corrfitter_results/"
beta=7.05
mass=-0.835
plt = plot()
E_min=1
E_max=2
plot_energy_levels!(plt,h5dir;beta,mass,group="pipi",marker=:rect,E_min,E_max)
#plot_energy_levels!(plt,h5dir;beta,mass,group="rho",marker=:diamond)
plot_energy_levels!(plt,h5dir;beta,mass,group="pi",marker=:circle)
#plot!(plt,ylims=(1.08,1.15))

