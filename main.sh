filename="out_scattering_I2"
loglist="./input/isospin_logfiles_list"

# activate correct julia environment and set up dependencies
#julia -e 'using Pkg; Pkg.activate("./scripts/src_jl"); Pkg.instantiate()'

#find ./input/  -name $filename  > $loglist
#python3 scripts/HDF5.py $loglist
#julia scripts/ensemble_table.jl
#julia scripts/average.jl
python3 scripts/fitting.py
