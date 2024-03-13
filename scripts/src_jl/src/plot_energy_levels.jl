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
function plot_energy_levels!(plt,h5dir;beta,mass,group,E_min=1,E_max=1,marker=:circle,showinf=false,factor=1)
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
    E  = factor .* hcat(E...)
    ΔE = factor .* hcat(ΔE...)
 
    # perform the plot
    for level in E_min:E_max
        label = isone(factor) ? "level $level ($group)" : "$factor × level $level ($group)"
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
        scatter!(plt,inv.(L),E[level,:],yerr=ΔE[level,:];marker,label) 
        if showinf
            xl = xlims(plt)
            xlims!(plt,(0.0,xl[2]))
        end
    end
end
function errorbars_plotting(eigvals,Δeigvals,ylim_lower)
    errorbar_lower = similar(Δeigvals)
    for i in eachindex(Δeigvals)
        if eigvals[i] - Δeigvals[i] < 0
            errorbar_lower[i] = eigvals[i] - ylim_lower
        else
            errorbar_lower[i] = Δeigvals[i]
        end
    end
    return errorbar_lower
end
function plot_correlator!(plt,t,eigvals,Δeigvals;kws...)
    # remove negative entries 
    plot_eigvals = replace(x -> x < 0 ? NaN : x, eigvals)
    # but check if they are compatible with zeros
    #@assert minimum(eigvals + Δeigvals) > 0 
    # find smallest error bar that remains on the positive side 
    ylim_lower_v1 = minimum(filter(x -> x >0, eigvals - Δeigvals))  
    # alternatively, take smallest positive value and lower by one order of magnitude
    ylim_lower_v2 = minimum(filter(isfinite, plot_eigvals)) / 10
    ylim_lower = min(ylim_lower_v1, ylim_lower_v2)
    # set up error bars for plot
    errorbar_lower = errorbars_plotting(eigvals,Δeigvals,ylim_lower) 
    # add data to existing plot
    scatter!(plt,t,plot_eigvals,yerrors=(errorbar_lower, Δeigvals);kws...)
    return plt
end
function add_mass_band!(plt,m,Δm;label="",alpha=0.5,kws...)
    hspan!(plt,[m+Δm,m-Δm];label,alpha,kws...)
end
function add_fit_range!(plt,tmin,tmax,E,ΔE;label="",kws...)
    plot!(plt,tmin:tmax, E*ones(length(tmin:tmax)), ribbon = ΔE; label, kws...)
end