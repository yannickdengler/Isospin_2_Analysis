filename="out_scattering_I2"
loglist="./input/isospin_logfiles_list"

find ./input/  -name $filename  > $loglist
python3 scripts/HDF5.py $loglist
julia scripts/ensemble_table.jl
julia scripts/average.jl
python3 scripts/fitting.py
