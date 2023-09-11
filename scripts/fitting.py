import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import os

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()]

def make_models(T,tmin=2):
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=T, tmin=tmin, tmax=abs(T)/2, a='a', b='a', dE='dE')]

def make_prior(N):
    prior = gv.BufferDict()
    # NOTE: use a log-Gaussion distrubtion for forcing positive energies
    # NOTE: Even with this code they can be recovered by providing loose priors of 0.1(1) for both
    prior['log(a)']  = gv.log(gv.gvar(N * ['1(1)']))
    prior['log(dE)'] = gv.log(gv.gvar(N * ['1(1)']))
    return prior

def bootstrap_fit(fitter,dset,T,n=20):
    pdatalist = (cf.process_dataset(ds, make_models(T)) for ds in gv.dataset.bootstrap_iter(dset, n=n))
    bs = gv.dataset.Dataset()
    for bsfit in fitter.bootstrapped_fit_iter(pdatalist=pdatalist):
        bs.append(E=np.cumsum(bsfit.pmean['dE']),a=bsfit.pmean['a'])
    bs = gv.dataset.avg_data(bs, bstrap=True)
    E = bs['E']
    a = bs['a']
    print('{:2}  {:15}  {:15}'.format('E', E[0], E[1]))
    print('{:2}  {:15}  {:15}'.format('a', a[0], a[1]))

def print_fit_param(fit):
        p = fit.p
        E = np.cumsum(p['dE'])
        a = p['a']
        logGBF = fit.logGBF # Gaussian Bayes Factor
        chi2 = fit.chi2     
        dof = fit.dof
        Q = fit.Q
        print('{:2}  {:15}  {:15}'.format('E', E[0], E[1]))
        print('{:2}  {:15}  {:15}'.format('a', a[0], a[1]))
        # TODO: Better printing
        #print('log(BGF) = ', logGBF)
        #print('Q = ', Q)
        print('chi2/dof = ', chi2/dof, '\n')

def main(data,T,Nmax=5):
    fitter = cf.CorrFitter(models=make_models(T))
    avg = gv.dataset.avg_data(data)
    p0 = None
    # TODO: find good Nmax
    for N in range(2,Nmax+1):
        print('nterm =', N, 30 * '=')
        prior = make_prior(N)
        fit = fitter.lsqfit(data=avg, prior=prior, p0=p0)
        print_fit_param(fit)
        p0 = fit.pmean
        #print(fit)
    print('bootstrap: ',30 * '=')
    # NOTE: A bootstrap fit can only be performed if`the object `fitter` has 
    # already been used to perform a fit.
    # NOTE: The bootstrap analysis is performed using the priors and initial 
    # parameters used in the last invokation of the previous fit. 
    bootstrap_fit(fitter, data, T)
    # NOTE: From the lsqfit documentation
    # There are several different views available for each plot, specified by parameter view:
    #   'ratio': Data divided by fit (default).
    #   'diff': Data minus fit, divided by data’s standard deviation.
    #   'std': Data and fit.
    #   'log': 'std' with log scale on the vertical axis.
    #   'loglog': ‘std’` with log scale on both axes.
    fit.show_plots(view='log')

filedir  = './output/HDF5_source_average/'
filelist = os.listdir(filedir)
filesrc  = filedir+filelist[0]
fid = h5py.File(filesrc,'r')

# read the data from the hdf5 file
T = get_hdf5_value(fid,'N_T')
corr = get_hdf5_value(fid,'correlator')
corr_deriv = get_hdf5_value(fid,'correlator_deriv')
ops = get_hdf5_value(fid,'operators')

# pion correlator data
# TODO: Assert that operators match
corr_pi   = corr[44,:,:]
corr_rho  = corr[46,:,:]
corr_pipi = -corr_deriv[48,:,:]
print(filesrc)

corr = dict(Gab=corr_pipi)
dset = gv.dataset.Dataset(corr)
main(dset,-T,Nmax=10)