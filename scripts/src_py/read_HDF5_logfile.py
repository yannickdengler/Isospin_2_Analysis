"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-7
    @Last Modified by: Yannick Dengler
    
    Functions that help with reading the HDF5_logfile
 """

import h5py

def get_corr_from_HDF5_logfile(filename):
    with h5py.File(filename,"r") as file:
        return file["correlators"][()]

def get_ops_from_HDF5_logfile(filename):
    with h5py.File(filename,"r") as file:
        return file["operators"][()]

def get_info_from_HDF5_logfile(filename):
    with h5py.File(filename,"r") as file:
        return file["N_L"][()], file["N_T"][()], file["gauge_group"][()].decode(), file["beta"][()], file["m_1"][()], file["m_2"][()]

def get_corr_ops_info_from_HDF5_logfile(filename):
    return (get_corr_from_HDF5_logfile(filename),get_ops_from_HDF5_logfile(filename),get_info_from_HDF5_logfile(filename))
     