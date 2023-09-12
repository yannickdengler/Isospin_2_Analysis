using Pkg; Pkg.activate("./scripts/src_jl")
using Plots
using HDF5

"""
    find_matching_files(files,beta,mass)

Given an interable containing filenames to hdf5 files, this function find all 
files that match the specified values of the inverse coupling `beta` and the 
fermion mass. 
"""
function find_matching_files(files;beta,mass)
    file_ids = h5open.(files) # opens all files 
    index = 1:length(files)   # indices of specific files
    # check elementwise for a match
    correct_beta = read.(file_ids,"pipi/beta") .== beta  
    correct_mass_m1 = read.(file_ids,"pipi/m_1") .== mass
    correct_mass_m2 = read.(file_ids,"pipi/m_1") .== mass
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

h5dir = "./output/HDF5_corrfitter_results/"
files = readdir(h5dir,join=true)

matched_files = find_matching_files(files,beta=7.2,mass=-0.78)

fid = h5open.(matched_files,"r")
E = read.(fid,"pipi/E")
ΔE = read.(fid,"pipi/Delta_E")
E_bs = read.(fid,"pipi/E_bs")
ΔE_bs = read.(fid,"pipi/Delta_E_bs")
T = read.(fid,"pipi/N_T")
L = read.(fid,"pipi/N_L")


n1=1
n2=2
plt = plot()
for i in eachindex(matched_files)
    energies = E_bs[i][n1:n2]
    Δenergies = ΔE_bs[i][n1:n2]
    Ls = L[i]*ones(n)
    scatter!(plt,Ls,energies,yerr=Δenergies,label="L=$(L[i]) T=$(T[i])") 
end
plot!(n=2,legend=:outerright)