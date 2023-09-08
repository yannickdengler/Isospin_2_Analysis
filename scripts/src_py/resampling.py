"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-7
    @Last Modified by: Yannick Dengler
    
    Resamples along Montecarlotime for different methods. "JK" only works with one obervable (or "JK_SAMEDIM" if they come from the same gauge configurations). Otherwise "BS_FIX" (Bootstrapping) works works for everything but is not deterministic ("BS_SAMEDIM" is similar to JK_SAMEDIM but with Bootstrapping). 
 """

import sys
import numpy as np
import random

def resampling(OG_Sample, sampling_args):                # OG_sample = [observable][mont], returns [resample][observable], [observable] means different observable
    """
    Takes a list of shape [num_observable][T_mont] and resamples it arcording to the sampling args. The result is a list of shape [num_resampling][num_bservables]
    """    
    if sampling_args[0] == "JK":
        if len(OG_Sample) > 1:
            sys.exit("resampling: Number of Corr can not be larger than 1 for JK, use JK_SAMEDIM instead!")
        return resampling_JK(OG_Sample, sampling_args)
    elif sampling_args[0] == "JK_SAMEDIM":
        for i in range(len(OG_Sample)):
            if (len(OG_Sample[i]) != len(OG_Sample[0])):
                sys.exit("Samples for JK_SAMEDIM are not of same DIM, use Bootstrap instead")
        return resampling_JK_SAMEDIM(OG_Sample, sampling_args)
    elif sampling_args[0] == "BS_SAMEDIM":
        for i in range(len(OG_Sample)):
            if len(OG_Sample[i]) is not len(OG_Sample[0]):
                sys.exit("Samples for BS_SAMEDIM are not of same DIM, use BS instead")
        return resampling_BS_SAMEDIM(OG_Sample, sampling_args)
    elif sampling_args[0] == "BS_DIFFDIM":
        return resampling_BS_DIFFDIM(OG_Sample, sampling_args)
    elif sampling_args[0] == "None":
        return [np.mean(OG_Sample,axis=1),]
    else:
        sys.exit("No valid sampling method given!")



def resampling_JK(OG_Sample, sampling_args):
    """
    Resamples an OG_Sample of one Observable with the cut-1 Jackknife method
    """
    Resamples = np.zeros((len(OG_Sample[0]), 1))
    num_JK = len(Resamples)
    for i in range(num_JK):
        for j in range(num_JK):
            if i != j:
                Resamples[i][0] += OG_Sample[0][i]/num_JK
    return Resamples

def resampling_JK_SAMEDIM(OG_Sample, sampling_args):
    """
    Resamples an OG_Sample of one Observable with the cut-1 Jackknife method with more than one Observable coming from the same gauge configurations
    """
    num_JK = len(OG_Sample[0])
    num_obs = len(OG_Sample)
    Resamples = np.zeros((num_JK, num_obs))
    for i in range(num_obs):
        for j in range(num_JK):
            for k in range(num_JK):
                if j != k:
                    Resamples[j][i] += OG_Sample[i][j]/num_JK
    return Resamples



def resampling_BS_SAMEDIM(OG_Sample, sampling_args):
    """
    Resamples an OG_Sample of Observables by bootstrapping num times, 
    """
    num_BS = sampling_args[1]
    num_obs = len(OG_Sample)
    num_mont = len(OG_Sample[0])

    Resamples = np.zeros((num_BS, num_obs))

    for i in range(num_mont):
        for j in range(num_BS):
            randint = random.randint(0,num_mont-1)
            for k in range(num_obs):
                Resamples[j][k] += OG_Sample[k][randint]/num_mont
    return Resamples

def resampling_BS_DIFFDIM(OG_Sample, sampling_args):
    """
    Resamples an OG_Sample of Observables by bootstrapping num times, 
    """
    num_BS = sampling_args[1]
    num_obs = len(OG_Sample)
    # num_mont = len(OG_Sample[0])                                # kann nicht generell definiert werden!

    Resamples = np.zeros((num_BS, num_obs))

    for i in range(num_BS):
        for j in range(num_obs):
            num_mont = len(OG_Sample[j])
            for k in range(num_mont):
                randint = random.randint(0,num_mont-1)
                Resamples[i][j] += OG_Sample[j][randint]/num_mont
    return Resamples




# example_sample = np.zeros(shape = (10, 89))

# # for sampling in ["JK_SAMEDIM", "JK", "BS_SAMEDIM", "BS_DIFFDIM", "None"]:
# for sampling in ["JK_SAMEDIM", "BS_SAMEDIM", "BS_DIFFDIM", "None"]:
#     print(sampling+": ",len(resampling(example_sample,[sampling, 100])))