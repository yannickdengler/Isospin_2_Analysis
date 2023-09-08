"""
        average_sources(hdf_file) 

    Takes an hdf5 file according and extract the correlator information from it.
    Then an average over all stochastic sources is performed. Instead of the 
    filename and array containing the correlator data can also be provided.

    NOTE: This assumes that the correlator consist of a sum of single diagrams.
    In case of a product of multiple diagrams this averaging is not necessarily 
    correct. E.g. for disconnected meson contributions this does not yield the 
    correct Euclidean time dependence. 
"""
function average_sources(hdf_file::AbstractString) 
    file_id = h5open(hdf_file)
    corr = read(file_id,"correlators")
    nsrc = read(file_id,"N_hits")
    close(file_id)
    return average_sources(corr)
end
function average_sources(corr::AbstractArray) 
    # According to our specification this correlator needs to have 5 indices
    # corresponding to (t,tMC,nsrc,nops)
    @assert ndims(corr) == 4
    nsrc = size(corr)[3]
    source_averaged = mean(corr,dims=3)
    # this removes the source-dimension from the array
    source_averaged = dropdims(source_averaged,dims=3)
    return source_averaged, nsrc
end

"""
        average_configurations(corr::AbstractArray)

    Average over all configurations and calculate the estimated variance. This 
    assumes that the average over all stochastic sources has been performed.

        average_configurations(hdf_file::AbstractArray)

    Averages over the correlator data from a specified hdf5 files. Assumes that
    stochastic source average has been performed. 

    It is assumed that the configurations are thermalised and independent, i.e. 
    that no autocorrelation exists.
"""
function average_configurations(corr::AbstractArray) 
    # According to our specification this correlator needs to have 4 indices
    # corresponding to (t,tMC,nmom,nops)
    @assert ndims(corr) == 3
    ncfg = size(corr)[2]
    avg_corr  = mean(corr,dims=2)
    Δavg_corr = std(corr,dims=2)./sqrt(ncfg)
    # this removes the source-dimension from the array
    avg_corr  = dropdims(avg_corr,dims=2)
    Δavg_corr = dropdims(Δavg_corr,dims=2)
    return avg_corr, Δavg_corr, ncfg
end
function average_configurations(hdf_file::AbstractString) 
    file_id = h5open(hdf_file)
    corr = read(file_id,"correlators")
    close(file_id)
    # If the array is 4-dimensional we assume that the average over all sources 
    # has already been performed
    @assert ndims(corr) == 3 || ndims(corr) == 4
    if ndims(corr) == 3
        return average_configurations(corr)
    elseif ndims(corr) == 4
        source_averaged, nsrc = average_sources(corr)
        return average_configurations(source_averaged)
    end
end
"""
    write_averaged_hdf5_files(hdf5list)

Given a list of hdf files that contain the data extracted from the logfiles, 
theis function performs the averaging over all stochastic sources and the 
average over all Monte-Carlo samples. 

The data is saved in the hdf5-format in the directories `hdf5dir_src_averaged`
and `hdf5dir_MC_averaged`. The directories are created at the same level as the
directory containing the initial hdf5 file. 

The source-averaged file contains the following data:
    - `correlator`: average over all stochastic sources.
    - `correlator_deriv`: numerical derivative of the data after averaging

The Monte-Carlo-averaged file contains
    - `correlator`: averaged correlator
    - `correlator_deriv`: averaged first derivative of the correlator
    - `correlator_2nd_deriv`: averaged first derivative of the correlator
    - `meff`: effective mass using an implicit definition
    - `meff_deriv`: effective mass of the first derivative of the correlator

Uncertainties are provided using the empirical standard deviation of the mean 
for the correlators and jackknife resampling for the effective masses. They are
labelled by the prefix "Delta_", i.e. the uncertainties of the correlator are 
saved as `Delta_correlator`.

"""
function write_averaged_hdf5_files(hdf5list;hdf5dir_src_averaged = "HDF5_source_average", hdf5dir_MC_averaged  = "HDF5_MC_average")
    for logfile in hdf5list
        file_id = h5open(logfile)
        corr = read(file_id,"correlators")

        keys_hdf = keys(file_id)
        keys_without_correlators = filter(!isequal("correlators"),keys(file_id))

        ispath(hdf5dir_src_averaged) || mkpath(hdf5dir_src_averaged)
        ispath(hdf5dir_MC_averaged)  || mkpath(hdf5dir_MC_averaged)
        filename = splitpath(logfile)[end]
        file_src = h5open(joinpath(hdf5dir_src_averaged,filename),"w")
        file_MC  = h5open(joinpath(hdf5dir_MC_averaged,filename),"w")

        # save data that has not been changed
        for k in keys_without_correlators
            file_src[k] = read(file_id,k)
            file_MC[k] = read(file_id,k)
        end

        source_averaged, nsrc = average_sources(corr) 
        source_averaged_deriv = correlator_derivative(source_averaged)
        file_src["correlator"]       = source_averaged 
        file_src["correlator_deriv"] = source_averaged_deriv

        # 2nd derivative as a check for convexity of diagonal correlators
        source_averaged_2ndderiv = correlator_derivative(correlator_derivative(source_averaged))

        c, Δc, ncfg = average_configurations(source_averaged)
        c_deriv, Δc_deriv, ncfg = average_configurations(source_averaged_deriv)
        meff, Δmeff = implicit_meff_jackknife(source_averaged;sign=+1)
        meff_deriv, Δmeff_deriv = implicit_meff_jackknife(source_averaged_deriv;sign=-1)
        c_2nd_deriv, Δc_2nd_deriv, ncfg = average_configurations(source_averaged_2ndderiv)

        file_MC["correlator"] = c 
        file_MC["Delta_correlator"] = Δc 
        file_MC["correlator_deriv"] = c_deriv
        file_MC["Delta_correlator_deriv"] = Δc_deriv
        file_MC["correlator_2nd_deriv"] = c_2nd_deriv
        file_MC["Delta_correlator_2nd_deriv"] = Δc_2nd_deriv
        file_MC["meff"] = meff
        file_MC["Delta_meff"] = Δmeff
        file_MC["meff_deriv"] = meff_deriv
        file_MC["Delta_meff_deriv"] = Δmeff_deriv

        close(file_id)
        close(file_src)
        close(file_MC)
    end
end
