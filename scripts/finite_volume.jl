using Pkg; Pkg.activate("./scripts/src_jl")
using I2julia
using HDF5

function finitevolume_goldstone(L,m)
    @. model(L,p)  = p[1]*(1+abs(p[2])*exp(-L*p[1])/abs(L*p[1])^(3/2))
    fit  = curve_fit(model,L,m,ones(2))
    return fit
end
function finitevolume(L,m,mGS_inf)
    @. model(L,p)  = p[1]*(1+abs(p[2])*exp(-L*mGS_inf)/abs(L*mGS_inf)^(3/2))
    fit  = curve_fit(model,L,m,ones(2))
    return fit
end