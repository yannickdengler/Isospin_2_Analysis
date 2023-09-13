"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-7
    @Last Modified by: Yannick Dengler
    
    Calculates energy levels from correlators
 """

import src_py.read_HDF5_logfile as read_HDF5_logfile
from scipy.optimize import curve_fit
from scipy.optimize import bisect
import numpy as np
import matplotlib.pyplot as plt


def calc_corr(data, args):
    result = {}
    result["Corr"] = []
    for i in range(len(data)):
        result["Corr"].append(data[i])
    return result

def calc_corr_tilde(data, args):
    result = {}
    result["Corr_tilde"] = []
    for i in range(len(data)-2):
        result["Corr_tilde"].append(data[i]-data[i+2])
    return result

def calc_eff_mass_impl_deri(data, args):
    def zero_eff_mass(eff_mass, ratio, index):
        return np.sinh(eff_mass*(T_2-index))/np.sinh(eff_mass*(T_2-(index+1))) - ratio    
    result = {}
    result["m_eff_impl_deri"] = []
    for i in range(len(data)-3):
        ratio = (data[i]-data[i+2])/(data[i+1]-data[i+3])
        T_2 = (len(data)-2)//2
        if (T_2-i) == 0 or (T_2-(i+1)) == 0:     
            result["m_eff_impl_deri"].append(0)
        else:
            res = bisect(f=zero_eff_mass, a=1e-30, b=1000, args = (ratio,i))
            if np.isnan(res):
                result["m_eff_impl_deri"].append(0)
            else:
                result["m_eff_impl_deri"].append(res)
    return result

def calc_convexity(data, args):
    result = {}
    result["convexity"] = []
    for i in range(len(data)-2):
        result["convexity"].append((data[i]-2*data[i+1]+data[i+2])/4)
    return result

def calc_eff_mass_log(data, args):
    result = {}
    result["eff_mass_log"] = []
    for i in range(len(data)-1):
        m_eff = np.log(data[i])/np.log(data[i+1])
        if np.isnan(m_eff) or np.isinf(m_eff):
            result["eff_mass_log"].append(0)
        else:
            result["eff_mass_log"].append(np.log(data[i])/np.log(data[i+1]))
    return result

def basic_analysis(data, args):
    result = {}
    for func in [calc_corr, calc_corr_tilde, calc_eff_mass_impl_deri, calc_convexity, calc_eff_mass_log]:
        result_tmp = func(data, args)
        for key in result_tmp:
            result[key] = result_tmp[key]
    return result


def basic_analysis_ops(data, args):
    Ops = args[0]
    N_T = args[1]
    result = {}
    for i, Op in zip(range(len(Ops)), Ops):
        result_tmp = basic_analysis(data[i*N_T:(i+1)*N_T], args = [])
        for key in result_tmp:
            result[key+"_%s"%Op] = result_tmp[key]
    return result



def fit_corr(data, args):
    result = {}
    corr_sinh = calc_corr_tilde(data, args)["Corr_tilde"]
    print(corr_sinh)
    filename = args[0]
    Op_ind = args[1]
    N_L, N_T, gauge_group, beta, m_1, m_2 = read_HDF5_logfile.get_info_from_HDF5_logfile(filename)
    fit_limits = np.genfromtxt("../input/fit_limits/Scattering_I2_%s_beta%1.2e_m1%1.3e_m2%1.3e_%i_%i"%(gauge_group, beta, m_1, m_2, N_T, N_L), dtype=int)
    fit_range = fit_limits[Op_ind]
    if not N_T == len(data):
        print("Warning! N_T is diffent in fit_corr")
    def corr_fit_func(n_t, E_0, A_0):
        return A_0*np.sinh((N_T/2.-n_t)*E_0)
    xdata = np.arange(1,N_T-1)
    popt, pcov = curve_fit(f=corr_fit_func, xdata=xdata[fit_range[0]:fit_range[1]+1], ydata=corr_sinh[fit_range[0]:fit_range[1]+1])
    # plt.scatter(xdata, corr_sinh)
    # plt.plot(xdata, corr_fit_func(corr_sinh, popt[0], popt[1]))
    # plt.scatter(xdata, [abs(ele) for ele in corr_sinh])
    # plt.plot(xdata, [abs(ele) for ele in corr_fit_func(xdata, popt[0], popt[1])])
    # plt.yscale("log")
    # plt.axvline(fit_range[0])
    # plt.axvline(fit_range[1]+1)
    # plt.show()
    # print(popt)
    result["E_0"] = [popt[0],]
    result["A_0"] = [popt[1],]
    return result