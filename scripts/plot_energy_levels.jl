using Pkg; Pkg.activate("./scripts/src_jl")
using Plots; gr(ms=7,frame=:box)
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
function add_pion_volume_extrapolation_plot!(h5dir;beta,mass)
    # get an estimate for infinite volume pion mass
    E, ΔE, L = I2julia.finite_volume_data(h5dir;beta,mass,group="pi")
    fit, model = finitevolume_goldstone(L,E,ΔE)
    prm, std = abs.(fit.param), abs.(stderror(fit))
    @assert fit.converged "fit did not converge"

    # get ribbons
    x = minimum(L):0.05:maximum(L)
    middle = model(x,prm)
    upper = model(x,prm.+std)
    lower = model(x,prm.-std)

    plot!(plt,inv.(x),middle, ribbon=(upper-middle,middle-lower),label = "FV fit")
    return plt
end
function add_mass_band!(plt,m,Δm;label="",alpha=0.5)
    hspan!(plt,[m+Δm,m-Δm];label,alpha)
end

h5dir = "./output/HDF5_corrfitter_results/"

ind = 15
ensemble_sets = unique_ensemble_sets(h5dir;group = "pi")
mass, beta, gauge_group = ensemble_sets[ind]

E_min, E_max = 1, 1
# naive proxy for infinite volume  mass
mπLmax, ΔmπLmax = mass_on_largest_volume(h5dir;beta,mass,group="pi")
EππLmax, ΔEππLmax = mass_on_largest_volume(h5dir;beta,mass,group="pipi")
# actual infinite volume extrapolation 
E, ΔE, L = I2julia.finite_volume_data(h5dir;beta,mass,group="pi")
fit, model = finitevolume_goldstone(L,E,ΔE)
mπinf, Δmπinf = fit.param[1], stderror(fit)[1]

plt = plot(legend=:outertopright,title="β=$beta, mass=$mass, gauge group = $gauge_group")
#plot_energy_levels!(plt,h5dir;beta,mass,group="pipi",marker=:rect,E_min,E_max)
#plot_energy_levels!(plt,h5dir;beta,mass,group="rho",marker=:diamond)
plot_energy_levels!(plt,h5dir;beta,mass,group="pi",marker=:circle,showinf=false)
add_pion_volume_extrapolation_plot!(h5dir;beta,mass)
add_mass_band!(plt,mπinf, Δmπinf;label="mpi(L → ∞)",alpha=0.5)
