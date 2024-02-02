using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia
using HDF5

function mass_on_largest_volume(h5dir;beta,mass,group)
    largets_volune = find_largest_volume(h5dir;beta,mass,group)
    fid = h5open(largets_volune)
    mπLmax  = first(read(fid,joinpath(group,"E")))
    ΔmπLmax = first(read(fid,joinpath(group,"Delta_E")))
    return mπLmax, ΔmπLmax
end
Δratio(x,y,Δx,Δy) = sqrt((Δx / y)^2 + (Δy*x / y^2)^2)
Δproduct(x,y,Δx,Δy) = sqrt((Δx * y)^2 + (Δy*x)^2)

h5dir = "./output/HDF5_corrfitter_results/"
files = readdir(h5dir,join=true)

io = open("output/overview_table.csv", "w");
#io = Base.stdout

println(io,"group, beta, mass, L, T, N, m_pi/m_rho, E_pipi/m_pi_Lmax, E_pipi, E_pipi rel.error in %, E_pipi/m_pi_Lmax rel.error in %")
for file in files

    fid = h5open(file)

    # Get number of configurations used
    N = length(fid["pipi/montecarlotimes"])

    # overall lattice parameters
    group = fid["pipi/gauge_group"][]
    beta = fid["pipi/beta"][]
    mass = fid["pipi/m_1"][]
    L = fid["pipi/N_L"][]
    T = fid["pipi/N_T"][]

    # Energy levels
    Eππ = fid["pipi/E"][1]
    Eπ = fid["pi/E"][1]
    Eρ = fid["rho/E"][1]

    ΔEππ = fid["pipi/Delta_E"][1]
    ΔEπ = fid["pi/Delta_E"][1]
    ΔEρ = fid["rho/Delta_E"][1]

    # pion mass on largest volume
    mπLmax, ΔmπLmax = mass_on_largest_volume(h5dir;beta,mass,group="pi")

    # chiral limit parameter
    ratioχ  = Eπ/Eρ 
    Δratioχ = Δratio(Eπ,Eρ,ΔEπ,ΔEρ)

    # closeness to threshold
    ratiothresh  = Eππ/mπLmax 
    Δratiothresh = Δratio(Eππ,mπLmax,ΔEππ,ΔmπLmax)

    # relative error
    Δrelratiothresh  = Δratiothresh/ratiothresh
    ΔrelEππ = ΔEππ/Eππ
    ΔrelEρ  = ΔEρ/Eρ
    ΔrelEρ  = ΔEρ/Eρ

    # quality parameters
    χ2dof = fid["pipi/chi2"][] / fid["pipi/dof"][]

    str_rχ = errorstring(ratioχ, Δratioχ)
    str_th = errorstring(ratiothresh,Δratiothresh)
    str_rl = round(100ΔrelEππ,sigdigits=2)
    str_ππ = errorstring(Eππ,ΔEππ)
    str_rel_thresh = round(100Δrelratiothresh,sigdigits=2)
  
    # print to files
    println(io,"$group, $beta, $mass, $L, $T, $N, $str_rχ, $str_th, $str_ππ, $str_rl, $str_rel_thresh")

    close(fid)
end
close(io)