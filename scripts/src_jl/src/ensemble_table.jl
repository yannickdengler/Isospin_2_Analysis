"""
        write_ensemble_list(hdf5list) 

    Take an iterable `hdf5list` containing the absolute, full filenames of the
    hdf5 files and creates a csv table of all corresponding  ensembles. 
"""
function write_ensemble_list(hdf5list;filename="ensembles.csv") 
    # open file for writing list of ensembles
    io = open(filename,"w")
    write(io,"β;m;L;T;nsrc;ncfg;name\n");


    # loop over all relevant
    for hdf_file in hdf5list
        file_id = h5open(hdf_file)

        T = read(file_id,"N_T")[1]
        L = read(file_id,"N_L")[1]
        β = read(file_id,"beta")[1]
        m1 = read(file_id,"m_1")[1]
        m2 = read(file_id,"m_2")[1]

        corr = read(file_id,"correlators")
        nsrc = size(corr)[3]
        ncfg = size(corr)[2]

        # we always assume degenerate fermions
        @assert m1 == m2

        write(io,"$β;$m1;$L;$T;$nsrc;$ncfg;$(basename(hdf_file))\n");
        close(file_id)
    end
    close(io)
    sort_ensemble_file!(filename)
end
function sort_ensemble_file!(filename)
    data, header = readdlm(filename,';',Any,header=true)
    data = sortslices(data,dims=1)
    writedlm(filename,vcat(header,data),';')
end
