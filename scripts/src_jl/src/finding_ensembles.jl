"""
    find_matching_files(files,beta,mass[,group])

Given an interable containing filenames to hdf5 files, this function finds all 
files that match the specified values of the inverse coupling `beta` and the 
fermion mass. 
"""
function find_matching_files(files;beta,mass,group="pi")
    file_ids = h5open.(files) # opens all files 
    index = 1:length(files)   # indices of specific files
    # check elementwise for a match
    correct_beta = read.(file_ids,joinpath(group,"beta")) .== beta  
    correct_mass_m1 = read.(file_ids,joinpath(group,"m_1")) .== mass
    correct_mass_m2 = read.(file_ids,joinpath(group,"m_1")) .== mass
    # set every index that does not match to zero, and filter all vanishing 
    # entries. Only the indices that correspond to matching ensembles remain.
    mask = correct_beta.*correct_mass_m1.*correct_mass_m2
    indices = filter!(!iszero,mask.*index)
    # get the filenames that matched
    matched_files = getindex(files,indices)
    # close all files that have been opened initially
    close.(file_ids)
    return matched_files
end
"""
    find_largest_volume(h5dir,beta,mass[,group])

Given a directory `h5dir` that contains hdf5 files, this function finds 
the file with the largest spatial extent in lattice units that matches the 
inverse coupling `beta` and the fermion mass `mass`. 
"""
function find_largest_volume(h5dir;beta,mass,group="pi")
    files = readdir(h5dir,join=true)
    matched = find_matching_files(files;beta,mass,group)
    file_ids = h5open.(matched)
    Ls = read.(file_ids,joinpath(group,"N_L"))
    ind = findmax(Ls)[2]
    return matched[ind]
end
"""
    unique_ensemble_sets(h5dir[;group])

Finds all unique sets of the fermion mass, inverse coupling and gauge group in 
the hdf5 files in the directory h5dir. Returns an array of tuples in the form
of `(mass, beta, gauge_group)`.

"""
function unique_ensemble_sets(h5dir;group="pi")
    files = readdir(h5dir,join=true)
    fids   = h5open.(files)
    m1 = read.(fids,joinpath(group,"m_1"))
    m2 = read.(fids,joinpath(group,"m_2"))
    beta = read.(fids,joinpath(group,"beta"))
    gauge_group = read.(fids,joinpath(group,"gauge_group"))
    @assert m1 == m2
    ensemble_sets = [(m1[i],beta[i],gauge_group[i]) for i in eachindex(fids)]
    unique!(ensemble_sets)
    return ensemble_sets
end