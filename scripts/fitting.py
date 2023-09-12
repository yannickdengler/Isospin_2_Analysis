import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import os
import matplotlib.pyplot as plt

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

def bootstrap_fit(fitter,dset,T,tmin,tmax,n=20):
    pdatalist = (cf.process_dataset(ds, make_models(T,tmin,tmax)) for ds in gv.dataset.bootstrap_iter(dset, n=n))
    bs = gv.dataset.Dataset()
    for bsfit in fitter.bootstrapped_fit_iter(pdatalist=pdatalist):
        bs.append(E=np.cumsum(bsfit.pmean['dE']),a=bsfit.pmean['a'])
    bs = gv.dataset.avg_data(bs, bstrap=True)
    E = bs['E']
    a = bs['a']
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

def main(data,T,Nmax=5,plotname="test",plotdir="./plots/",antisymmetric=False):
    tmin = 3
    tmax = abs(T/2) - 1
    T = - abs(T) if antisymmetric else abs(T) 
    fitter = cf.CorrFitter(models=make_models(T,tmin,tmax))
    avg = gv.dataset.avg_data(data)
    p0 = None
    # TODO: find good Nmax
    for N in range(2,Nmax+1):
        print('nterm =', N, 30 * '=')
        prior = make_prior(N)
        fit = fitter.lsqfit(data=avg, prior=prior, p0=p0)
        p0 = fit.pmean
        #print(fit)
        print_fit_param(fit)
        #if N == 2:
    print_fit_param(fit)
    E, a, chi2, dof = first_fit_parameters(fit) 
    print('bootstrap: ',30 * '=')
    # NOTE: A bootstrap fit can only be performed if`the object `fitter` has 
    # already been used to perform a fit.
    # NOTE: The bootstrap analysis is performed using the priors and initial 
    # parameters used in the last invokation of the previous fit. 
    E_bs, a_bs = bootstrap_fit(fitter, data, T, tmin, tmax)
    # NOTE: From the lsqfit documentation
    # There are several different views available for each plot, specified by parameter view:
    #   'ratio': Data divided by fit (default).
    #   'diff': Data minus fit, divided by data’s standard deviation.
    #   'std': Data and fit.
    #   'log': 'std' with log scale on the vertical axis.
    #   'loglog': ‘std’` with log scale on both axes.
    #fit.show_plots(view='log')
    os.makedirs(plotdir+plotname, exist_ok=True)
    fit.show_plots(view='ratio',save=plotdir+plotname+'/ratio.pdf')
    plt.close()
    fit.show_plots(view='log'  ,save=plotdir+plotname+'/data.pdf')
    plt.close()

filedir  = './output/HDF5_source_average/'
filelist = os.listdir(filedir)

for i in range(0,len(filelist)):
    filesrc  = filedir+filelist[i]
    fid = h5py.File(filesrc,'r')

    # read the data from the hdf5 file
    T = get_hdf5_value(fid,'N_T')
    L = get_hdf5_value(fid,'N_L')
    m = get_hdf5_value(fid,'m_1')
    beta = get_hdf5_value(fid,'beta')
    corr = get_hdf5_value(fid,'correlator')
    corr_deriv = get_hdf5_value(fid,'correlator_deriv')
    ops = get_hdf5_value(fid,'operators')

    plotname = "beta{}_m{}_L{}_T{}".format(beta,m,L,T)
    print(plotname)

    # pion correlator data
    # TODO: Assert that operators match
    corr_pi   = corr[44,:,:]
    corr_rho  = corr[46,:,:]
    corr_pipi = -corr_deriv[48,:,:]
    
    corr = dict(Gab=corr_pipi)
    dset = gv.dataset.Dataset(corr)
    main(dset,T,Nmax=10,plotname=plotname,antisymmetric=True)