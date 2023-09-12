"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-8
    @Last Modified by: Yannick Dengler
    
    Takes data and parses it to the resampler and the functions to get an estimation of the error
"""

import src_py.error_classes as ec
import src_py.resampling as rs
import h5py
import numpy as np

def return_corr(data, args):
    """
    data contains the correlators, args contains 
    """
    result = {}
    for i in range(len(data)):
        result["corr_%i"%(i)] = [data[i],]
    return result

def get_corr_from_HDF5_logfile(filename):
    with h5py.File(filename,"r") as file:
        return file["correlators"][()]

def get_ops_from_HDF5_logfile(filename):
    with h5py.File(filename,"r") as file:
        return file["operators"][()]

def get_info_from_HDF5_logfile(filename):
    with h5py.File(filename,"r") as file:
        # for key in file:
        #     print(key, file[key][()])
        return file["N_L"][()], file["N_T"][()], file["gauge_group"][()].decode(), file["beta"][()], file["m_1"][()], file["m_2"][()]
     
def get_corr_ops_info_from_HDF5_logfile(filename):
    return (get_corr_from_HDF5_logfile(filename),get_ops_from_HDF5_logfile(filename),get_info_from_HDF5_logfile(filename))

def calc_and_save_corr_with_error_test(filename):
    (corr, ops, info) = get_corr_ops_info_from_HDF5_logfile(filename)
    test = ec.measurement("corr", measure_func=return_corr, sampling_args = ("JK_SAMEDIM", 0, 0))
    test.measure(orig_sample=np.swapaxes(np.mean(corr[0], axis=0),0,1), args=[])
    test.print_to_HDF()

# calc_and_save_corr_with_error_test("/home/dengler_yannick/Documents/Isospin_2_analysis/HDF5_logfiles/Scattering_I2_SP(4)_beta6.900_m1-0.900_m2-0.900_T20_L10_logfile.hdf5")


# test = ec.measurement("corr", measure_func=return_corr, sampling_args = ("JK_SAMEDIM", 0, 0))
# test.read_from_HDF()
# test.print_everything()

