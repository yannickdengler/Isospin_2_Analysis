"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-7
    @Last Modified by: Yannick Dengler
    
    Calculates m_inf from a given set of E(L)
 """
 
import sys
from scipy.optimize import curve_fit
import numpy as np

def inf_mass_fit_Goldstone(N_L, m_inf, A_M):
    """
    Function that describes finite volume effects for the extrapolation of energy levels of Goldstone bosons to the infinite volume.    
    """
    return m_inf*(1+A_M*np.exp(-m_inf*N_L)/(m_inf*N_L)**(3/2.))
def inf_mass_fit_meson(N_L, m_inf, A_M):
    """
    Function that describes finite volume effects for the extrapolation of energy levels of non-Goldstone bosons to the infinite volume.    
    """
    return m_inf*(1+A_M*np.exp(-mass_Goldstone*N_L)/(mass_Goldstone*N_L)**(3/2.))

def extrapolate_to_inf_mass(data, args):
    """
    Function does the inifinite volume extrapolation for a given set of E(L). Returns: m_inf, A_M  
    """
    N_Ls = args[0]
    mass_Goldstone = args[1]
    E_0s = data
    results = {}
    if (len(N_Ls) != len(E_0s)):
        sys.exit('"extrapolate_to_inf_mass: len of arrays N_L and E_0 is not the same!"')
    if mass_Goldstone == 0:
        popt, pcov = curve_fit(f=inf_mass_fit_Goldstone, xdata=N_Ls, ydata=E_0s)
    else:
        popt, pcov = curve_fit(f=inf_mass_fit_meson, xdata=N_Ls, ydata=E_0s)
    results["m_inf"] = [popt[0],]
    results["A_M"] = [popt[1],]
    return results  
