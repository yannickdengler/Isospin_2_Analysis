import gvar as gv
import corrfitter as cf
import h5py
import numpy as np
import lsqfit
import matplotlib

fn = '/home/zierler_fabian/Nextcloud/Singlets_Data/data.h5'
fn = '/home/fabian/Downloads/data.h5'
f  = h5py.File(fn,'r')

h5key = 'runsSp4/Lt24Ls12beta6.9m1-0.89m2-0.89/out_spectrum/DEFAULT_SEMWALL TRIPLET_g5'
corr = dict(Gab=f[h5key])
dset = gv.dataset.Dataset(corr)

print(f[h5key].shape[1])

def make_models():
    """ Create corrfitter model for G(t). """
    return [cf.Corr2(datatag='Gab', tp=24, a='a', b='a', dE='dE')]

def make_prior(N):
    prior = gv.BufferDict()
    # NOTE: use a log-Gaussion distrubtion for forcing positive energies
    # NOTE: Even with this code they can be recovered by providing loose priors of 0.1(1) for both
    prior['log(a)']  = gv.log(gv.gvar(N * ['1(1)']))
    prior['log(dE)'] = gv.log(gv.gvar(N * ['1(1)']))
    return prior

def bootstrap_fit(fitter,dset,n=20):
    pdatalist = (cf.process_dataset(ds, make_models()) for ds in gv.dataset.bootstrap_iter(dset, n=n))
    bs = gv.dataset.Dataset()
    for bsfit in fitter.bootstrapped_fit_iter(pdatalist=pdatalist):
        bs.append(E=np.cumsum(bsfit.pmean['dE']),a=bsfit.pmean['a'])
    bs = gv.dataset.avg_data(bs, bstrap=True)
    E = bs['E']
    a = bs['a']
    print('{:2}  {:15}  {:15}'.format('E', E[0], E[1]))
    print('{:2}  {:15}  {:15}'.format('a', a[0], a[1]))

def main(data):
    fitter = cf.CorrFitter(models=make_models())
    avg = gv.dataset.avg_data(data)
    p0 = None
    for N in [2, 3, 4, 5]:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)
        fit = fitter.lsqfit(data=avg, prior=prior, p0=p0)
        print(fit)
        p0 = fit.pmean
    print(30 * '=', 'bootstrap')
    # NOTE: A bootstrap fit can only be performed if`the object `fitter` has 
    # already been used to perform a fit.
    bootstrap_fit(fitter, data)

main(dset)