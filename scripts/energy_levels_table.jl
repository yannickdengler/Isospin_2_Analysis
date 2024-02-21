using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia
using HDF5
using LaTeXStrings

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

io1 = open("output/overview_table.csv", "w");
io2 = open("output/overview_table_machine_readable.csv", "w");
io3 = open("output/energy_levels.csv", "w");
io4 = open("output/energy_levels_tex.csv", "w");
#io = Base.stdout

println(io1,"group, beta, mass, L, T, N_conf, m_pi/m_rho, E_pipi/m_pi_Lmax, E_pipi, E_pipi rel.error in %, E_pipi/m_pi_Lmax rel.error in %")
println(io2,"group, beta, mass, L, T, N_conf, m_pi/m_rho, Delta_m_pi/m_rho, m_pi, Delta_m_pi, E_pipi, Delta_E_pipi")
println(io3,"group, beta, mass, L, T, N_conf, m_pi/m_rho, Delta_m_pi/m_rho, m_pi, Delta_m_pi, E_pipi, Delta_E_pipi")
println(io4,L"$\beta$ & $a m_{0}$ & $N_L$ & $N_T$ & $n_{\rm config}$ & $ m_\pi/m_\rho$ & $m_\pi$ & $E_{\pi\pi}$ \\\\")
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

    # unrenormalized pion decay constant and plaquette

    # pion mass on largest volume
    mπLmax, ΔmπLmax = mass_on_largest_volume(h5dir;beta,mass,group="pi")
    same_bare_parameters = find_matching_files(readdir(h5dir,join=true);beta,mass,group="pi")
    N_levels = length(same_bare_parameters)

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
    str_π  = errorstring(Eπ,ΔEπ)
    str_rel_thresh = round(100Δrelratiothresh,sigdigits=2)

    # print to files
    println(io1,"$group, $beta, $mass, $L, $T, $N, $str_rχ, $str_th, $str_ππ, $str_rl, $str_rel_thresh")
    println(io2,"$group, $beta, $mass, $L, $T, $N, $ratioχ, $Δratioχ, $Eπ, $ΔEπ, $Eππ, $ΔEππ")

    # write only if we have at least three levels and the gauge group is Sp(4)
    if group == "SP(4)" && N_levels > 2    
        println(io3,"$group, $beta, $mass, $L, $T, $N, $ratioχ, $Δratioχ, $Eπ, $ΔEπ, $Eππ, $ΔEππ")
        println(io4,"$beta & $mass & $L & $T & $N & $str_rχ & $str_π & $str_ππ \\\\")
    end

    close(fid)
end
close(io1)
close(io2)
close(io3)
close(io4)