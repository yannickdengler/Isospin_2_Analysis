using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia
using HDF5

function write_g0g5_correlator_to_hdf5(hdf5correlators,logfile_list)

    files = readlines(logfile_list)
    fid   = h5open(hdf5correlators,"cw")

    # hdf5 group names in the putput logfile
    identifiers = keys(fid)

    # get the matching ensemble names from the VSC4 naming scheme
    first_cnfg_filenames = first.(read.(Ref(fid),joinpath.(identifiers,"filenames"))) 
    ensemble_names = getindex.(splitpath.(first_cnfg_filenames),8)


    for file in files
        fπ_correlator = parse_spectrum(file,"DEFAULT_SEMWALL TRIPLET")["g0g5"]
        fπ_correlator = permutedims(fπ_correlator) # flip to match vonvention of logfile

        # get the coorrspondng key for that particular measurement
        ensemble_id = findfirst(x->occursin(x,file),ensemble_names)
        group = identifiers[ensemble_id]

        # write correlator to file
        fid[joinpath(group,"g0g5_correlator")] = fπ_correlator

    end
    close(fid)
end

logfile_list = "output/isospin_logfiles_list"
hdf5correlators = "output/correlators.hdf5"
write_g0g5_correlator_to_hdf5(hdf5correlators,logfile_list)
