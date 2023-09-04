"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-4
    @Last Modified by: Yannick Dengler
    
    This file contains scripts to create a HDF5 from an output file from the HiRep Scattering Code
    Execute with: "python3 HDF5.py PATH_TO_LOGFILE"
 """


import numpy as np
import h5py
import os
import sys

def create_scattering_momentum(filename):
    print("create_scattering: ", filename)

    gauge_group = ""
    beta = 0
    m_1 = 0
    m_2 = 0
    N_L = 0
    N_T = 0
    ptotx = 0                                                 
    ptoty = 0
    ptotz = 0
    nmom = 0
    Correlators = []                                            # 4-array, Operator (+Semwall etc.), src, Montecarlo-time, Lattice-time
    Acceptance = []                                             # ??
    Filenames = []                                              # Vector of Strings, filenames including the montecarlo time
    Operators = []                                              # The measured Operators (pi1, rho1, AD etc. )
    Momenta = []                                                # The measured Momenta
    Montecarlotimes = []
    Plaquette = []
    num_src = 0

    fi = open(filename)
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
        if nmom == 0:
            if words[0] == "[MAIN][0]Number":
                if words[1] == "of":
                    if words[2] == "momenta":
                        nmom = int(words[3])
        if ptotx == 0 and ptoty == 0 and ptotz == 0:
            if words[0] == "[MAIN][0]Total":
                if words[1] == "momentum:":
                    ptotx = int(words[2])
                    ptoty = int(words[3])
                    ptotz = int(words[4])

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
                        # print(Montecarlotimes[len(Montecarlotimes)-1])
                        # print(len(Montecarlotimes))
        if words[0] == "[IO][0]Configuration":
            if words[7][:10] == "Plaquette=":
                Plaquette.append(float(words[7][10:]))
        if words[0][:7] == "[IO][0]":
            for i in range(5,20):
                if words[0][i:i+5] == "_src_":
                    Operator = words[0][7:i]
                    if Operator not in Operators:
                        Operators.append(Operator)
            if len(words[0]) == 7:
                P = [int(words[1]), int(words[2]), int(words[3])]
                if P not in Momenta:
                    Momenta.append(P)
            if len(words[0]) == 8:
                P = [int(words[0][7]), int(words[1]), int(words[2])]
                if P not in Momenta:
                    Momenta.append(P)

    Operators_w_im = []

    for Operator in Operators:
        Operators_w_im.append(Operator)
        Operators_w_im.append(Operator+"_im")
    
    print(Operators_w_im)
    print(len(Operators_w_im))
    print(len(Montecarlotimes))
    print(num_src)
    print(Momenta)

    print("Writing Correlators:")
    Correlators = np.zeros((len(Operators_w_im), len(Momenta), num_src, len(Montecarlotimes),N_T))
    for lines in data:
        words = lines.split()
        if words[0][:7] == "[IO][0]":
            for i in range(5,20):
                if words[0][i:i+5] == "_src_":
                    current_Operator_index = Operators_w_im.index(words[0][7:i])
                    # print(current_Operator_index, words[0][7:i])
                    for j in range(4):
                        if words[0][i+j+5:i+j+9] == "_run":
                            current_src_index = int(words[0][i+5:i+j+5])
                    for i in range(1,10):
                        if words[0][len(words[0])-i] == "n":
                            current_Montecarlotime_index = Montecarlotimes.index(int(words[0][len(words[0])-i+1:]))
        if words[0][:7] == "[IO][0]":
            if len(words[0]) == 7:
                P = [int(words[1]), int(words[2]), int(words[3])]
                Correlators[current_Operator_index][Momenta.index(P)][current_src_index][current_Montecarlotime_index][int(words[4])] = float(words[5])     #max(float(words[4]),1)
                Correlators[current_Operator_index+1][Momenta.index(P)][current_src_index][current_Montecarlotime_index][int(words[4])] = float(words[6])   #max(float(words[5]),1)
            if len(words[0]) == 8:
                P = [int(words[0][7]), int(words[1]), int(words[2])]
                Correlators[current_Operator_index][Momenta.index(P)][current_src_index][current_Montecarlotime_index][int(words[3])] = float(words[4])     #max(float(words[4]),1)
                Correlators[current_Operator_index+1][Momenta.index(P)][current_src_index][current_Montecarlotime_index][int(words[3])] = float(words[5])   #max(float(words[5]),1)
        
    print(len(Correlators))
    print(len(Correlators[0]))
    print(len(Correlators[0][0]))
    print(len(Correlators[0][0][0]))
    print(len(Correlators[0][0][0][0]))

    f = h5py.File("./HDF5_files/Scattering_mom_%s_%1.2e_%1.3e_%1.3e_%i_%i_%i_(%i%i%i).hdf5"%(gauge_group, beta,m_1,m_2,N_T,N_L,nmom,ptotx,ptoty,ptotz),"w")

    Operator_set = f.create_dataset("Operators", data=Operators_w_im)
    ptot_set = f.create_dataset("P_tot", data=[ptotx,ptoty,ptotz])
    Momenta_set = f.create_dataset("Momenta", data=Momenta)
    Montecarlotimes_set = f.create_dataset("Montecarlotimes", data = Montecarlotimes)
    Plaquette_set = f.create_dataset("Plaquette", data = Plaquette)
    gauge_group_set = f.create_dataset("gauge_group", data = gauge_group)
    beta_set = f.create_dataset("beta",(1,), data = beta)
    m_1_set = f.create_dataset("m_1",(1,), data = m_1)
    m_2_set = f.create_dataset("m_2",(1,), data = m_2)
    N_L_set = f.create_dataset("L",(1,), data = N_L)
    N_T_set = f.create_dataset("T",(1,), data = N_T)
    nmom_set = f.create_dataset("nmom",(1,), data = nmom)
    ptot_set = f.create_dataset("ptot", data = (ptotx,ptoty,ptotz))
    Filename_set = f.create_dataset("Filenames", data = Filenames)
    Correlator_set = f.create_dataset("Correlators", data = Correlators)

    f.close()

create_scattering_momentum(sys.argv[1])