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

def main(data,T,tmin,tmax,Nmax,plotting=False,printing=False):
    T = abs(T) 
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
        fit.show_plots(view='ratio')
        fit.show_plots(view='log'  )
    return E, a, E_bs, a_bs, chi2, dof

def fit_single_file(tmin=1,delta_tmax=0,Nmax=5):
    hdf5file = './input/data_connected.hdf5'
    f = h5py.File(hdf5file)
    for run in f.keys():
        for ens in f[run].keys():
            T    = get_hdf5_value(f[run][ens],"lattice")[0]
            L    = get_hdf5_value(f[run][ens],"lattice")[1]
            corr = get_hdf5_value(f[run][ens],"g0g5")*L**3/2
            plaq = get_hdf5_value(f[run][ens],"plaquette")
            beta = get_hdf5_value(f[run][ens],"beta")
            mass = get_hdf5_value(f[run][ens],"quarkmasses")[0]

            dset = gv.dataset.Dataset(dict(Gab=corr))
            p = gv.dataset.avg_data(plaq)

            # renormalization from lattice perturbation theory 
            ZA = 1 + (5/4)*(-12.82-3)*8/(16*np.pi**2)/(beta*p)

            tmax = T/2 - delta_tmax
            
            E, a, E_bs, a_bs, chi2, dof = main(dset,T,tmin,tmax,Nmax,plotting=False)
            fpi     = a_bs[0]*np.sqrt(2/E[0])
            fpi_ren = ZA*a_bs[0]*np.sqrt(2/E[0])

            print(ens,",",T,",",L,",",mass,",",beta,",",E[0].mean,",",E[0].sdev,",",fpi_ren.mean,",",fpi_ren.sdev)


print("ensemble,T,L,m0,beta,mpi,mpi_err,fpi,fpi_err")
fit_single_file()