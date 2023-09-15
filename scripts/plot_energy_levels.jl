using Pkg; Pkg.activate("./scripts/src_jl")
using Plots; gr(ms=7,frame=:box)
using HDF5
using I2julia

function find_largest_volume(h5dir;beta,mass,group)
    files = readdir(h5dir,join=true)
    matched = I2julia.find_matching_files(files;beta,mass,group)
    file_ids = h5open.(matched)
    Ls = read.(file_ids,joinpath(group,"N_L"))
    ind = findmax(Ls)[2]
    return matched[ind]
end
function mass_on_largest_volume(h5dir;beta,mass,group)
    largets_volune = find_largest_volume(h5dir;beta,mass,group)
    fid = h5open(largets_volune)
    mπLmax  = first(read(fid,joinpath(group,"E")))
    ΔmπLmax = first(read(fid,joinpath(group,"Delta_E")))
    return mπLmax, ΔmπLmax
end

h5dir = "./output/HDF5_corrfitter_results/"
h5dir = "./output/HDF5_corrfitter_results_ncg2/"
h5dir = "./output/HDF5_corrfitter_results_ncg3/"
beta=5.4
mass=-0.89
plt = plot(legend=:outertopright)
E_min=1
E_max=1
mπLmax, ΔmπLmax = mass_on_largest_volume(h5dir;beta,mass,group="pi")
EππLmax, ΔEππLmax = mass_on_largest_volume(h5dir;beta,mass,group="pipi")

plot_energy_levels!(plt,h5dir;beta,mass,group="rho",marker=:diamond)
plot_energy_levels!(plt,h5dir;beta,mass,group="pipi",marker=:rect,E_min,E_max)
plot_energy_levels!(plt,h5dir;beta,mass,group="pi",marker=:circle,showinf=false)

hspan!(plt,2 .* [mπLmax+ΔmπLmax,mπLmax-ΔmπLmax],label="2 mpi(L_max)")
hspan!(plt,4 .* [mπLmax+ΔmπLmax,mπLmax-ΔmπLmax],label="4 mpi(L_max)")
plot!(plt,ylims=(0.8*mπLmax,4.1*mπLmax))
