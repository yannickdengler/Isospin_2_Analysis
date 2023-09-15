"""
    find_matching_files(files,beta,mass)

Given an interable containing filenames to hdf5 files, this function finds all 
files that match the specified values of the inverse coupling `beta` and the 
fermion mass. 
"""
function find_matching_files(files;beta,mass,group)
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
    plot_energy_levels(h5dir;beta,mass,group,E_min=1,E_max=1)
    plot_energy_levels!(plt,h5dir;kws...)

Plot the energy levels of all ensembles with the specified inverse coupling 
`beta` and bare mass `mass`. The data needs to be stored in the form of hdf5
files containing the corrfitter output in the directory `h5dir` in the group
`group`.

The plot will contain all energy levels between `E_min` and `E_max`. By default 
only the lowest energy state is plotted.

"""
plot_energy_levels(h5dir;kws...) = plot_energy_levels!(plot(),h5dir;kws...)
function plot_energy_levels!(plt,h5dir;beta,mass,group,E_min=1,E_max=1,marker=:circle,showinf=false)
    files = readdir(h5dir,join=true)
    matched_files = find_matching_files(files;beta,mass,group)

    # read the relevant data from the hdf files
    fid = h5open.(matched_files,"r")
    E  = read.(fid,joinpath(group,"E"))
    ΔE = read.(fid,joinpath(group,"Delta_E"))
    T = read.(fid,joinpath(group,"N_T"))
    L = read.(fid,joinpath(group,"N_L"))
    close.(fid)

    # sort data by spatial extent
    perm = sortperm(L)
    permute!(ΔE,perm)
    permute!(E,perm)
    permute!(T,perm)
    permute!(L,perm)

    # convert array of array into a 2-dimensional array, i.e. matrix
    E = hcat(E...)
    ΔE = hcat(ΔE...)
 
    # perform the plot
    for level in E_min:E_max
        label = "level $level ($group)"
        ylabel = "energy"
        xlabel = "1/L"
        linestyle = :dash
        xticks_val = inv.(L)
        xticks_lab = ["1/$l" for l in L]
        if showinf
            push!(xticks_val,0.0)
            push!(xticks_lab,"inf")
        end
        xticks = (xticks_val,xticks_lab)
        plot!(;xlabel,ylabel,linestyle,xticks)
        plot!(plt,inv.(L),E[level,:],yerr=ΔE[level,:];marker,label) 
        if showinf
            xl = xlims(plt)
            xlims!(plt,(0.0,xl[2]))
        end
    end
end
