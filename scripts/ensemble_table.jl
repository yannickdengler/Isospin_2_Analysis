using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia

write_ensemble_list("output/logfiles.hdf5";outdir="output/tables",filename="pipi_fitintervals_default.csv")
