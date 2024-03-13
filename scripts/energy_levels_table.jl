using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia
using HDF5
using LaTeXStrings
using Statistics

function write_table_1(h5file;tabledir="output/tables")
    ispath(tabledir) || mkpath(tabledir)

    Δratio(x,y,Δx,Δy) = sqrt((Δx / y)^2 + (Δy*x / y^2)^2)
    Δproduct(x,y,Δx,Δy) = sqrt((Δx * y)^2 + (Δy*x)^2)

    fid0  = h5open(h5file,"r")
    ensembles = keys(fid0)

    io1 = open("$tabledir/results.csv", "w");
    io2 = open("$tabledir/results_readable.csv", "w");
    io4 = open("$tabledir/energy_levels_tex.csv", "w");
    #io = Base.stdout

    println(io1,"group, beta, mass, L, T, N_conf, m_pi/m_rho, E_pipi, fpi, plaq")
    println(io2,"group, beta, mass, L, T, N_conf, m_pi/m_rho, Delta_m_pi/m_rho, m_pi, Delta_m_pi, E_pipi, Delta_E_pipi, fpi, Delta_fpi, plaq, Delta_plaq")
    println(io4,L"$\beta$ & $a m_{0}$ & $N_L$ & $N_T$ & $n_{\rm config}$ & $ m_\pi/m_\rho$ & $m_\pi$ & $E_{\pi\pi}$ & $f_{\pi}$ & $\langle P \rangle$ \\\\ \hline \hline")

    for ensemble in ensembles

        fid = fid0[ensemble]
        # Get number of configurations used
        N = length(fid["pipi/montecarlotimes"])
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

        # plaquette and pion decay constant
        fπ  = fid["pi/fpi_ren"][]
        Δfπ = fid["pi/Delta_fpi_ren"][]
        plaq= fid["pi/plaquette"][]
        p, Δp = mean(plaq), std(plaq)/sqrt(length(plaq))

        # chiral limit parameter
        ratioχ  = Eπ/Eρ 
        Δratioχ = Δratio(Eπ,Eρ,ΔEπ,ΔEρ)
        χ2dof = fid["pipi/chi2"][] / fid["pipi/dof"][]

        str_rχ = errorstring(ratioχ, Δratioχ)
        str_ππ = errorstring(Eππ,ΔEππ)
        str_π  = errorstring(Eπ,ΔEπ)
        str_fπ = errorstring(fπ,Δfπ)
        str_p  = errorstring(p,Δp)

        # print to files
        println(io1,"$group, $beta, $mass, $L, $T, $N, $str_rχ, $str_ππ, $str_fπ, $str_p")
        println(io2,"$group, $beta, $mass, $L, $T, $N, $ratioχ, $Δratioχ, $Eπ, $ΔEπ, $Eππ, $ΔEππ, $fπ, $Δfπ, $p, $Δp")  
        println(io4,"$beta & $mass & $L & $T & $N & $str_rχ & $str_π & $str_ππ  & $str_fπ & $str_p \\\\")
    end
    close(io1)
    close(io2)
    close(io4)
    close(fid0)
end
write_table_1("output/fitresults.hdf5")
#h5open("output/fitresults.hdf5","r")