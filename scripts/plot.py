"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-6
    @Last Modified by: Yannick Dengler
    
    Plot Correlator and effective masses etc.
 """
 
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect                                            # should be done somewhere else in the future
import src_py.read_HDF5_logfile as read_HDF
import src_py.error_classes as err
import src_py.read_HDF5_logfile as HDF_log 
# print(plt.rcParams.keys())

def set_errorbar_settings():
    plt.rcParams["errorbar.capsize"] = 5
    # plt.solid_capstyle = "projecting"                                     # not needed!! dumb
    plt.rcParams["lines.linestyle"] = ""
    plt.rcParams["lines.markersize"] = 10

def get_mean_std_corr(Corr):
    Corr_src = np.mean(Corr, axis=1)
    return np.mean(Corr_src, axis=1), np.std(Corr_src, axis=1)

def Corr_settings():
    plt.yscale("log")
    plt.grid()
    plt.xlabel("$n_t$")
    plt.ylabel("C")

def eff_mass_settings():
    plt.grid()
    plt.xlabel("$n_t$")
    plt.ylabel("$m_{eff}$")


def calc_eff_mass_impl_deri(Corr):   
    m_eff = np.zeros(len(Corr)-3)
    def zero_eff_mass(eff_mass, ratio, index):
        return np.sinh(eff_mass*(T_2-index))/np.sinh(eff_mass*(T_2-(index+1))) - ratio
    for i in range(len(Corr)-3):
        ratio = (Corr[i]-Corr[i+2])/(Corr[i+1]-Corr[i+3])
        T_2 = (len(Corr)-2)//2
        if (T_2-i) == 0 or (T_2-(i+1)) == 0:                    # Only happens for sinh. For cosh both values are well-defined
            m_eff[i] = 0
        else:
            res = bisect(f=zero_eff_mass, a=1e-30, b=1000, args = (ratio,i))
            if np.isnan(res):
                m_eff[i] = 0
            else:
                m_eff[i] = bisect(f=zero_eff_mass, a=1e-30, b=1000, args = (ratio,i))
    return m_eff

# def plot_Corr(filename):
#     Corr = read_HDF.get_corr_from_HDF5_logfile(filename)
#     Operators = read_HDF.get_ops_from_HDF5_logfile(filename)
#     N_L, N_T, gauge_group, beta, m_1, m_2 = read_HDF.get_info_from_HDF5_logfile(filename)

#     (Corr, Corr_err) = get_mean_std_corr(Corr)
#     Corr_settings()
#     plt.title(r"%s, $\beta$=%1.1e, $m_1$=%1.1e, $m_2$=%1.1e"%(gauge_group, beta, m_1, m_2))
#     for ind in (0,44):
#         Corr_Op = Corr[ind]
#         Corr_err_Op = Corr_err[ind]
#         xarr = np.arange(len(Corr_Op))
#         plt.errorbar(xarr, Corr_Op, yerr = Corr_err_Op)
#     plt.show()

def plot_eff_mass(filename):
    (corr, ops, info) = HDF_log.get_corr_ops_info_from_HDF5_logfile(filename)
    N_L, N_T, gauge_group, beta, m_1, m_2 = info
    test = err.measurement("basic_%s_beta_%1.3f_m1_%1.3f_m2_%1.3f_T%i_L%i"%(gauge_group, beta, m_1, m_2, N_T, N_L), measure_func=None, sampling_args = ("JK_SAMEDIM", 0, 0))
    test.read_from_HDF()
    test.print_everything()
    eff_mass_settings()
    xarr = np.arange(1.5,N_T-1.5)
    print(len(xarr), len(test.results["m_eff_impl_deri_"+"pipi"].median))
    for op in ("pi", "rho", "pipi"):
        plt.errorbar(xarr, test.results["m_eff_impl_deri_"+op].median, (test.results["m_eff_impl_deri_"+op].ep,test.results["m_eff_impl_deri_"+op].em), label = op)
    plt.legend()
    plt.show()

def plot_corr(filename):
    (corr, ops, info) = HDF_log.get_corr_ops_info_from_HDF5_logfile(filename)
    N_L, N_T, gauge_group, beta, m_1, m_2 = info
    test = err.measurement("basic_%s_beta_%1.3f_m1_%1.3f_m2_%1.3f_T%i_L%i"%(gauge_group, beta, m_1, m_2, N_T, N_L), measure_func=None, sampling_args = ("JK_SAMEDIM", 0, 0))
    test.read_from_HDF()
    test.print_everything()
    eff_mass_settings()
    xarr = np.arange(N_T)
    # print(len(xarr), len(test.results["m_eff_impl_deri_"+"pipi"].median))
    plt.yscale("log")
    for op in ("pi", "rho", "pipi"):
        plt.errorbar(xarr, test.results["Corr_"+op].median, (test.results["Corr_"+op].ep,test.results["Corr_"+op].em), label = op)
    plt.legend()
    plt.show()



    # Corr = read_HDF.get_corr_from_HDF5_logfile(filename)
    # Operators = read_HDF.get_ops_from_HDF5_logfile(filename)
    # N_L, N_T, gauge_group, beta, m_1, m_2 = info

    # (Corr, Corr_err) = get_mean_std_corr(Corr)
    # eff_mass_settings()
    # plt.title(r"%s, $\beta$=%1.1e, $m_1$=%1.1e, $m_2$=%1.1e"%(gauge_group, beta, m_1, m_2))
    # for ind in (0,44):
    #     Corr_Op = Corr[ind]
    #     xarr = np.arange(len(Corr_Op)-3)
    #     m_eff = calc_eff_mass_impl_deri(Corr_Op)
    #     m_eff_err = m_eff/10.                                            # TOY MODEL
    #     plt.errorbar(xarr, m_eff, m_eff_err)
    # plt.show()

# def plot_inf_volume_extrapolation(filename):                    # work in Progress (plot inf volume extrapolation from E(L) and fit params with error)



# set_errorbar_settings()

# plot_Corr("output/HDF5_logfiles/Scattering_I2_SP(4)_beta7.200_m1-0.780_m2-0.780_T24_L12_logfile.hdf5")
# plot_eff_mass("../output/HDF5_logfiles/Scattering_I2_SP(4)_beta6.900_m1-0.900_m2-0.900_T24_L12_logfile.hdf5")
