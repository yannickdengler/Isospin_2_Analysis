import gvar as gv
import corrfitter as cf
import h5py
import numpy as np

def get_hdf5_value(hdf5file,key):
    return hdf5file[key][()
                         ]

def make_models(T):
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=T, a='a', b='a', dE='dE')]

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
        print('chi2/dof = ', chi2/dof)

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


fn = '/home/zierler_fabian/Nextcloud/Isospin_2_Analysis/output/HDF5_source_average/Scattering_I2_SP(4)_beta6.900_m1-0.870_m2-0.870_T20_L10_logfile.hdf5'
f  = h5py.File(fn,'r')

# read the data from the hdf5 file
T = get_hdf5_value(f,'N_T')
corr = get_hdf5_value(f,'correlator')
corr_deriv = get_hdf5_value(f,'correlator_deriv')
ops = get_hdf5_value(f,'operators')

# pion correlator data
corr_pi   = corr[44,:,:]
corr_rho  = corr[46,:,:]
corr_pipi = corr_deriv[48,:,:]
# TODO: Assert that operators match
print(ops[44])
print(ops[46])
print(ops[48])
print(T)

corr = dict(Gab=corr_pi)
dset = gv.dataset.Dataset(corr)

main(dset,T)