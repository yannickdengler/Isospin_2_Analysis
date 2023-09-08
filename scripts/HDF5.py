"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-5
    @Last Modified by: Fabian Zierler 2023-Sep-6
    
    This file contains scripts to create a HDF5 from an output file from the HiRep Scattering Code
    Execute with: "python3 HDF5.py PATH_TO_LISTFILE"
 """

import numpy as np
import h5py
import os
import sys

def add_pi_rho_pipi_average_I2(Operators, Correlators, N_L):
    """
    Calculates the meaned and normalized correlation functions of the pi, rho and pipi Operators in the isospin-2 channel and adds them to the "Operator" list and the "Correlators" array
    """
    pi_rho_pipi = ("pi", "rho", "pipi")
    for Op in pi_rho_pipi:
        Operators.append(Op)
        Operators.append(Op+"_im")
    rho_tmp = []
    rho_im_tmp = []
    for Op in ("rho1_11", "rho1_12", "rho1_13", "rho1_21", "rho1_22", "rho1_23", "rho1_31", "rho1_32", "rho1_33", "rho2_11", "rho2_12", "rho2_13", "rho2_21", "rho2_22", "rho2_23", "rho2_31", "rho2_32", "rho2_33"):
        rho_tmp.append(Correlators[Operators.index(Op)]/(2*pow(N_L,3))) 
        rho_im_tmp.append(Correlators[Operators.index(Op+"_im")]/(2*pow(N_L,3))) 
    rho = np.mean(rho_tmp,axis=0)
    rho_im = np.mean(rho_tmp,axis=0)
    pi_tmp = []
    pi_im_tmp = []
    for Op in ("pi1", "pi2"):
        pi_tmp.append(Correlators[Operators.index(Op)]/(2*pow(N_L,3))) 
        pi_im_tmp.append(Correlators[Operators.index(Op+"_im")]/(2*pow(N_L,3))) 
    pi = np.mean(pi_tmp,axis=0)
    pi_im = np.mean(pi_im_tmp,axis=0)
    pipi = (Correlators[Operators.index("AD")]+Correlators[Operators.index("BC")])/(0.5*4*pow(N_L,6))
    pipi_im = (Correlators[Operators.index("AD_im")]+Correlators[Operators.index("BC_im")])/(0.5*4*pow(N_L,6))
    for Corr in (pi, pi_im, rho, rho_im, pipi, pipi_im):
        Correlators = np.append(Correlators, np.expand_dims(Corr, axis=0), axis = 0)
    return Operators, Correlators

def create_scattering(filename,hdfpath="../output/HDF5_logfiles/"):
    """
    Function that converts a logfile from a scattering measurement of HiRep to a HDF file with the relevant information
    """
    print("create_scattering: ", filename)

    logfile_name = filename
    for i in range(len(filename)-1):
        if filename[i] == "/":
            logfile_name = filename[i+1:]
    gauge_group = ""
    isospin_channel = int(999)
    beta = 0
    m_1 = 0
    m_2 = 0
    N_L = 0
    N_T = 0
    Acceptance = []                                             # ??
    Filenames = []                                              # Vector of Strings, filenames including the montecarlo time
    Operators = []                                              # The measured Operators (pi1, rho1, AD etc. )
    Montecarlotimes = []
    Plaquette = []
    num_src = 0

    with open(filename.strip()) as fi:                          # .strip() removes a potential "\n" in the end of the string
        data = fi.readlines()

        current_Operator_index = -1
        current_Montecarlotime = -1
        current_src = -1

        for lines in data:
            words = lines.split()
            if num_src == 0:
                if words[0] == "[CORR][0]Number":
                    if words[5] == "nhits":
                        num_src = int(words[7])
            if gauge_group == "":
                if words[0] == "[SYSTEM][0]Gauge":
                    gauge_group = words[2]
            if isospin_channel == 999:
                if words[0] == "[MAIN][0]Isospin":
                    isospin_channel = int(words[2])

            if beta == 0 and m_1 == 0 and m_2 == 0 and N_L == 0 and N_T == 0:
                if words[0] == "[MAIN][0]Configuration":
                    beta_index = 0
                    m1_index = 0
                    m2_index = 0
                    Lt_index = 0
                    Ls_index = 0
                    end_index = 0
                    for i in range(len(words[2])):
                        if words[2][i:i+2] == "Lt":    
                            Lt_index = i
                        if words[2][i:i+2] == "Ls":     
                            Ls_index = i
                        if words[2][i:i+4] == "beta":
                            beta_index = i
                        if words[2][i:i+2] == "m1" and m1_index == 0:
                            m1_index = i
                        if words[2][i:i+2] == "m2":
                            m2_index = i
                        if words[2][i:i+8] == "/configs":
                            end_index = i
                    print(words[2][Lt_index+2:Ls_index], words[2][Ls_index+2:beta_index], words[2][beta_index+4:m1_index], words[2][m1_index+2:m2_index], words[2][m2_index+2:end_index])
                    N_T = int(words[2][Lt_index+2:Ls_index])
                    N_L = int(words[2][Ls_index+2:beta_index])
                    beta = float(words[2][beta_index+4:m1_index])
                    m_1 = float(words[2][m1_index+2:m2_index])
                    m_2 = float(words[2][m2_index+2:end_index])
            if words[0] == "[MAIN][0]Configuration":
                if words[1] == "from":
                    Filenames.append(words[2])
                    for i in range(1,10):
                        if words[2][len(words[2])-i] == "n":
                            Montecarlotimes.append(int(words[2][len(words[2])-i+1:]))
            if words[0] == "[IO][0]Configuration":
                if words[7][:10] == "Plaquette=":
                    Plaquette.append(float(words[7][10:]))
            if words[0][:7] == "[IO][0]":
                for i in range(5,20):
                    if words[0][i:i+5] == "_src_":
                        Operator = words[0][7:i]
                        if Operator not in Operators:
                            Operators.append(Operator)

        Operators_w_im = []

        for Operator in Operators:
            Operators_w_im.append(Operator)
            Operators_w_im.append(Operator+"_im")
        
        print("Number of Motecarlo steps: ", len(Montecarlotimes))
        print("Number of sources: ", num_src)
        print("Isospin channel: ", isospin_channel)

        print("Writing Correlators...")

        Correlators = np.zeros((len(Operators_w_im), num_src, len(Montecarlotimes),N_T))
        for lines in data:
            words = lines.split()
            if words[0][:7] == "[IO][0]":
                for i in range(5,20):
                    if words[0][i:i+5] == "_src_":
                        current_Operator_index = Operators_w_im.index(words[0][7:i])
                        for j in range(4):
                            if words[0][i+j+5:i+j+9] == "_run":
                                current_src_index = int(words[0][i+5:i+j+5])
                        for i in range(1,10):
                            if words[0][len(words[0])-i] == "n":
                                current_Montecarlotime_index = Montecarlotimes.index(int(words[0][len(words[0])-i+1:]))
            if words[0][:7] == "[IO][0]":
                if len(words[0]) == 8:
                    Correlators[current_Operator_index][current_src_index][current_Montecarlotime_index][int(words[3])] = float(words[4])     #max(float(words[4]),1)
                    Correlators[current_Operator_index+1][current_src_index][current_Montecarlotime_index][int(words[3])] = float(words[5])   #max(float(words[5]),1)
            
        isospin_str = ""
        if isospin_channel == 2:
            isospin_str = "_I2"
            (Operators_w_im, Correlators) = add_pi_rho_pipi_average_I2(Operators_w_im, Correlators, N_L)

        print(len(Operators_w_im), " operators (w/ imag): ", Operators_w_im)
        print("Size of Correlator array [num_Operators][num_soruces][num_Montecarlotimes][N_T]:[%i][%i][%i][%i]"%(len(Correlators),len(Correlators[0]),len(Correlators[0][0]),len(Correlators[0][0][0])) )
        num_Montecarlotimes = len(Correlators[0][0])
        num_src = len(Correlators[0])


        os.makedirs(hdfpath, exist_ok=True)
        f = h5py.File(hdfpath+"Scattering%s_%s_beta%1.3f_m1%1.3f_m2%1.3f_T%i_L%i_logfile.hdf5"%(isospin_str, gauge_group, beta, m_1, m_2, N_T, N_L),"w")

        f.create_dataset("logfile name", data=logfile_name)
        f.create_dataset("isospin_channel", data=isospin_channel)
        f.create_dataset("N_mont", data = num_Montecarlotimes)
        f.create_dataset("N_hits", data = num_src)
        f.create_dataset("filenames", data = Filenames)
        f.create_dataset("plaquette", data = Plaquette)
        f.create_dataset("operators", data=Operators_w_im)
        f.create_dataset("montecarlotimes", data = Montecarlotimes)
        f.create_dataset("gauge_group", data = gauge_group)
        f.create_dataset("beta", data = beta)
        f.create_dataset("m_1", data = m_1)
        f.create_dataset("m_2", data = m_2)
        f.create_dataset("N_L", data = N_L)
        f.create_dataset("N_T", data = N_T)
        f.create_dataset("correlators", data = Correlators)
        print()


fi = open(sys.argv[1])
filelist = fi.read().splitlines()

for files in filelist:
    create_scattering(files)