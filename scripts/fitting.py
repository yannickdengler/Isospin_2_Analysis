import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import os
import matplotlib.pyplot as plt
import csv

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()]

def make_models(T,tmin,tmax):
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=T, tmin=tmin, tmax=tmax, a='a', b='a', dE='dE')]

def make_prior(N):
    prior = gv.BufferDict()
    # NOTE: use a log-Gaussion distrubtion for forcing positive energies
    # NOTE: Even with this code they can be recovered by providing loose priors of 0.1(1) for both
    prior['log(a)']  = gv.log(gv.gvar(N * ['1(1)']))
    prior['log(dE)'] = gv.log(gv.gvar(N * ['1(1)']))
    return prior

def bootstrap_fit(fitter,dset,T,tmin,tmax,n=20,printing=False):
    pdatalist = (cf.process_dataset(ds, make_models(T,tmin,tmax)) for ds in gv.dataset.bootstrap_iter(dset, n=n))
    bs = gv.dataset.Dataset()
    for bsfit in fitter.bootstrapped_fit_iter(pdatalist=pdatalist):
        bs.append(E=np.cumsum(bsfit.pmean['dE']),a=bsfit.pmean['a'])
    bs = gv.dataset.avg_data(bs, bstrap=True)
    E = bs['E']
    a = bs['a']
    if printing:
        print('bootstrap: ',30 * '=')
        print('{:2}  {:15}  {:15}'.format('E', E[0], E[1]))
        print('{:2}  {:15}  {:15}'.format('a', a[0], a[1]))
    return E, a

def first_fit_parameters(fit):
    p = fit.p
    E = np.cumsum(p['dE'])
    a = p['a']
    chi2 = fit.chi2     
    dof = fit.dof
    return E, a, chi2, dof

def print_fit_param(fit):
    E, a, chi2, dof = first_fit_parameters(fit) 
    print('{:2}  {:15}  {:15}'.format('E', E[0], E[1]))
    print('{:2}  {:15}  {:15}'.format('a', a[0], a[1]))
    print('chi2/dof = ', chi2/dof, '\n')

def main(data,T,tmin,tmax,Nmax,plotname="test",plotdir="./plots/",antisymmetric=False,plotting=False,printing=False):
    T = - abs(T) if antisymmetric else abs(T) 
    fitter = cf.CorrFitter(models=make_models(T,tmin,tmax))
    avg = gv.dataset.avg_data(data)
    p0 = None
    # TODO: find good Nmax
    for N in range(2,Nmax+1):
        prior = make_prior(N)
        fit = fitter.lsqfit(data=avg, prior=prior, p0=p0)
        p0 = fit.pmean
        if printing:
            print('nterm =', N, 30 * '=')
            #print(fit)
            print_fit_param(fit)
    E, a, chi2, dof = first_fit_parameters(fit) 
    # NOTE: A bootstrap fit can only be performed if`the object `fitter` has 
    # already been used to perform a fit.
    # NOTE: The bootstrap analysis is performed using the priors and initial 
    # parameters used in the last invokation of the previous fit. 
    E_bs, a_bs = bootstrap_fit(fitter, data, T, tmin, tmax)
    # NOTE: From the lsqfit documentation
    # There are several different views available for each plot, specified by parameter view:os.
    #   'ratio': Data divided by fit (default).
    #   'diff': Data minus fit, divided by data’s standard deviation.
    #   'std': Data and fit.
    #   'log': 'std' with log scale on the vertical axis.
    #   'loglog': ‘std’` with log scale on both axes.
    if plotting:
        os.makedirs(plotdir+plotname, exist_ok=True)
        fit.show_plots(view='ratio',save=plotdir+plotname+'/ratio.pdf')
        fit.show_plots(view='log'  ,save=plotdir+plotname+'/data.pdf')
    return E, a, E_bs, a_bs, chi2, dof

def save_corrfitter_results(fid,outfile_id,ensemble,group,E,a,E_bs,a_bs,chi2,dof,antisymmetric,Nmax,tmin,tmax,binsize):
    f = outfile_id

    E_mean = [E_i.mean for E_i in E]
    E_sdev = [E_i.sdev for E_i in E]
    E_bs_mean = [E_i.mean for E_i in E_bs]
    E_bs_sdev = [E_i.sdev for E_i in E_bs]
    a_mean = [a_i.mean for a_i in a]
    a_sdev = [a_i.sdev for a_i in a]
    a_bs_mean = [a_i.mean for a_i in a_bs]
    a_bs_sdev = [a_i.sdev for a_i in a_bs]

    for key in fid.keys():
        if "correlator" not in key:
            f.create_dataset(ensemble+group+key, data = get_hdf5_value(fid,key))

    f.create_dataset(ensemble+group+"E", data = E_mean)
    f.create_dataset(ensemble+group+"E_bs", data = E_bs_mean)
    f.create_dataset(ensemble+group+"A", data = a_mean)
    f.create_dataset(ensemble+group+"A_bs", data = a_bs_mean)
    f.create_dataset(ensemble+group+"Delta_E", data = E_sdev)
    f.create_dataset(ensemble+group+"Delta_E_bs", data = E_bs_sdev)
    f.create_dataset(ensemble+group+"Delta_A", data = a_sdev)
    f.create_dataset(ensemble+group+"Delta_A_bs", data = a_bs_sdev)
    f.create_dataset(ensemble+group+"antisymmetric", data = antisymmetric)
    f.create_dataset(ensemble+group+"chi2", data = chi2)
    f.create_dataset(ensemble+group+"dof", data = dof)
    f.create_dataset(ensemble+group+"Nexp", data = Nmax)
    f.create_dataset(ensemble+group+"tmin", data = tmin)
    f.create_dataset(ensemble+group+"tmax", data = tmax)
    f.create_dataset(ensemble+group+"binsize", data = binsize)
    return

def fit_all_files(infile,outfile,groups,tmins,tmaxs,binsize=1):

    fid = h5py.File(infile,'r')
    oid = h5py.File(outfile, 'w')
    
    for i in range(0,len(groups)):

        group = groups[i] 
        tmin = tmins[i]
        tmax = tmaxs[i]

        # read the data from the hdf5 file
        T    = get_hdf5_value(fid,group+'N_T')
        L    = get_hdf5_value(fid,group+'N_L')
        m    = get_hdf5_value(fid,group+'m_1')
        beta = get_hdf5_value(fid,group+'beta')
        corr = get_hdf5_value(fid,group+'correlator')
        ops  = get_hdf5_value(fid,group+'operators')
        corr_deriv = get_hdf5_value(fid,group+'correlator_deriv')

        plotname = "beta{}_m{}_L{}_T{}".format(beta,m,L,T)
        print(plotname)

        # start with pipi correlator
        corr_pipi = -corr_deriv[48,:,:]
        corr_pipi = dict(Gab=corr_pipi)
        dset = gv.dataset.Dataset(corr_pipi,binsize=binsize)

        antisymmetric = True
        plotdir = "./plots/"
        Nmax = 10

        E, a, E_bs, a_bs, chi2, dof = main(dset,T,tmin,tmax,Nmax,plotname,plotdir,antisymmetric)
        save_corrfitter_results(fid[group],oid,group,"pipi/",E,a,E_bs,a_bs,chi2,dof,antisymmetric,Nmax,tmin,tmax,binsize)

        # then fit both the pion and the vector meson
        corr_pi = corr[44,:,:]
        corr_rho = corr[46,:,:]
        corr_pi = dict(Gab=corr_pi)
        dset_pi = gv.dataset.Dataset(corr_pi,binsize=binsize)
        corr_rho = dict(Gab=corr_rho)
        dset_rho = gv.dataset.Dataset(corr_rho,binsize=binsize)

        antisymmetric = False
        plotdir = "./plots/"
        Nmax = 10
        tmin = 1
        tmax = T/2

        E, a, E_bs, a_bs, chi2, dof = main(dset_pi,T,tmin,tmax,Nmax,plotname,plotdir,antisymmetric)
        save_corrfitter_results(fid[group],oid,group,"pi/",E,a,E_bs,a_bs,chi2,dof,antisymmetric,Nmax,tmin,tmax,binsize)
        E, a, E_bs, a_bs, chi2, dof = main(dset_rho,T,tmin,tmax,Nmax,plotname,plotdir,antisymmetric)
        save_corrfitter_results(fid[group],oid,group,"rho/",E,a,E_bs,a_bs,chi2,dof,antisymmetric,Nmax,tmin,tmax,binsize)


def read_filelist_fitparam(parameterfile):
    reader = csv.reader(open(parameterfile))
    # create list that contain the fitting information
    tmins  = []
    tmaxs  = []
    groups = []
    # skip line containing headers
    next(reader, None)
    for row in reader:
        groups.append(row[6])
        tmins.append(int(row[7]))
        tmaxs.append(int(row[8]))

    return groups, tmins, tmaxs

parameterfile  = './input/pipi_fitintervals.csv'
groups, tmins, tmaxs = read_filelist_fitparam(parameterfile)

infile  = './output/correlators.hdf5'
outfile = './output/fitresults.hdf5'
fit_all_files(infile,outfile,groups,tmins,tmaxs)