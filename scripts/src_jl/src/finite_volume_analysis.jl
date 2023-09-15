function finitevolume_goldstone(L,m,Δm)
    @. model(L,p)  = p[1]*(1+abs(p[2])*exp(-L*p[1])/abs(L*p[1])^(3/2))
    fit  = curve_fit(model,L,m,inv.(Δm.^2),ones(2))
    return fit
end
function finitevolume(L,m,Δm,mGS_inf)
    @. model(L,p)  = p[1]*(1+abs(p[2])*exp(-L*mGS_inf)/abs(L*mGS_inf)^(3/2))
    fit  = curve_fit(model,L,m,inv.(Δm.^2),ones(2))
    return fit
end
function finite_volume_data(h5dir;beta,mass,group)
    files = readdir(h5dir,join=true)
    matched = I2julia.find_matching_files(files;beta,mass,group)
    fids = h5open.(matched)
    L = read.(fids,joinpath(group,"N_L"))
    E = first.(read.(fids,joinpath(group,"E")))
    ΔE = first.(read.(fids,joinpath(group,"Delta_E")))
    close.(fids)
    return E, ΔE, L
end
function all_infinite_volume_goldstones(h5dir)
    group = "pi"
    ensemble_sets = unique_ensemble_sets(h5dir;group = "pi")
    nsets = length(ensemble_sets)
    minf  = zeros(nsets)
    Δminf = zeros(nsets)
    for ind in eachindex(ensemble_sets)
        mass, beta, gauge_group = ensemble_sets[ind]
        E, ΔE, L = finite_volume_data(h5dir;beta,mass,group)
        length(E) < 2 && continue
        fit = finitevolume_goldstone(L,E,ΔE)
        minf[ind] = fit.param[1]
        Δminf[ind] = stderror(fit)[1]
    end
    ensemble_sets, minf, Δminf
end