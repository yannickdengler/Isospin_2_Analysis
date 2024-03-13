loglistI2="./output/isospin_logfiles_I2_list"
loglist="./output/isospin_logfiles_list"

# activate correct julia environment and set up dependencies
julia -e 'using Pkg; Pkg.activate("./scripts/src_jl")'
#julia -e 'using Pkg; Pkg.activate("./scripts/src_jl"); Pkg.instantiate()'

#find ./input/  -name "out_scattering_I2" > $loglistI2
#find ./input/  -name "out_spectrum" > $loglist
#python3 scripts/HDF5.py $loglistI2

#julia scripts/ensemble_table.jl
#julia scripts/average.jl
#julia scripts/write_fpi_correlator.jl
#python3 scripts/fitting.py
julia scripts/energy_levels_table.jl
