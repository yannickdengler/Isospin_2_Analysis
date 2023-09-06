"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-5
    @Last Modified by: Yannick Dengler
    
    Calculate phase shifts from energy levels and m_inf
 """

def generalized_mom(data, args):
    result = {}
    E = data[0]
    m_meson = data[1]
    if E < 2*m_meson:
        print("Energy below 2m threshhold!!")
        result["generalized_mom"] = float("inf")
    else:
        result["generalized_mom"] =  2*np.arcsin(np.sqrt(0.5*(np.cosh(E/2.)-np.cosh(m_meson))))
    return result

def tan_phase_shift(data, args):
    result = {}
    E = data[0]
    m_meson = data[1]
    N_L = args[0]
    P = generalized_mom(E, m_PS, N_L)
    q = P*N_L/(2*np.pi)
    tan_PS = np.pi**(3/2.)*q/Zeta(q**2)
    if E < 2*m_PS:
        result["tan_PS"] = float("inf")
        result["P"] = float("inf")
    else:
        result["tan_PS"] = tan_PS
        result["P"] = P
    return result





