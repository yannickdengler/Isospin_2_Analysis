using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia

h5dir = "./output/HDF5_corrfitter_results/"
ensemble_sets, minf, Δminf = all_infinite_volume_goldstones(h5dir)