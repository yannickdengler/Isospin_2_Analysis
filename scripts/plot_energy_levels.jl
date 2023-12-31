using Pkg; Pkg.activate("./scripts/src_jl")
using Plots; gr(ms=7,frame=:box)
using Plots; plotlyjs(ms=7,frame=:box)
using HDF5
using LsqFit
using I2julia

function mass_on_largest_volume(h5dir;beta,mass,group)
    largets_volune = find_largest_volume(h5dir;beta,mass,group)
    fid = h5open(largets_volune)
    mπLmax  = first(read(fid,joinpath(group,"E")))
    ΔmπLmax = first(read(fid,joinpath(group,"Delta_E")))
    return mπLmax, ΔmπLmax
end
function add_pion_volume_extrapolation_plot!(plt,h5dir;beta,mass,factor=1,skip=0)
    # get an estimate for infinite volume pion mass
    E, ΔE, L = I2julia.finite_volume_data(h5dir;beta,mass,group="pi")
    # skip smalles lattices
    range = 1+skip:length(L)
    E, ΔE, L = E[range], ΔE[range], L[range] 
    fit, model = finitevolume_goldstone(L,E,ΔE)
    prm, std = abs.(fit.param), abs.(stderror(fit))
    @assert fit.converged "fit did not converge"

    # get ribbons
    x = minimum(L):0.05:maximum(L)
    middle = factor .* model(x,prm)
    upper  = factor .* model(x,prm.+std)
    lower  = factor .* model(x,prm.-std)

    plot!(plt,inv.(x),middle, ribbon=(upper-middle,middle-lower),label = "FV fit")
    return plt
end
function add_mass_band!(plt,m,Δm;label="",alpha=0.5)
    hspan!(plt,[m+Δm,m-Δm];label,alpha)
end
function all_finite_volume_plots(h5dir,ensemble_sets)
    for ind in eachindex(ensemble_sets)
        mass, beta, gauge_group = ensemble_sets[ind]

        plotname = "$(gauge_group)_beta$(beta)_mass$(mass)"
        E_min, E_max = 1, 1
        skip = 0
        
        #infinite volume extrapolation 
        E, ΔE, L = I2julia.finite_volume_data(h5dir;beta,mass,group="pi")
        length(E) - skip < 3 && continue
        fit, model = finitevolume_goldstone(L[skip+1:end],E[skip+1:end],ΔE[skip+1:end])
        mπinf, Δmπinf = fit.param[1], stderror(fit)[1]

        plt = plot(legend=:outertopright,title="β=$beta, mass=$mass, gauge group = $gauge_group")
        plot_energy_levels!(plt,h5dir;beta,mass,group="pipi",marker=:rect,E_min,E_max)
        plot_energy_levels!(plt,h5dir;beta,mass,group="pi",marker=:circle,showinf=false,factor=2)
        add_pion_volume_extrapolation_plot!(plt,h5dir;beta,mass,factor=2,skip)
        add_mass_band!(plt,2mπinf,2Δmπinf;label="2mpi(L → ∞)",alpha=0.5)
        ispath("./plots/finite_volume/") || mkpath("./plots/finite_volume/") 
        savefig(plt,"./plots/finite_volume/$plotname.pdf")
        display(plt)
    end
end

h5dir = "./output/HDF5_corrfitter_results_optim/"
ensemble_sets = unique_ensemble_sets(h5dir;group = "pi")
all_finite_volume_plots(h5dir,ensemble_sets)

#plot_energy_levels!(plt,h5dir;beta,mass,group="rho",marker=:diamond)
#plot_energy_levels!(plt,h5dir;beta,mass,group="pi",marker=:circle,showinf=false,factor=1)
#add_mass_band!(plt,mπinf,Δmπinf;label="mpi(L → ∞)",alpha=0.5)
#add_mass_band!(plt,4mπinf,4Δmπinf;label="4mpi(L → ∞)",alpha=0.5)