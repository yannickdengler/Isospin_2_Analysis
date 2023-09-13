using Pkg; Pkg.activate("./scripts/src_jl")
using Plots
using HDF5

h5dir = "./output/HDF5_corrfitter_results/"
plot_egery_levels(h5dir;beta=7.2,mass=-0.75,group="pipi")
