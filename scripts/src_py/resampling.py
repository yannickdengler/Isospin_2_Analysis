"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-7
    @Last Modified by: Yannick Dengler
    
    Resamples along Montecarlotime for different methods. "JK" only works with one obervable (or "JK_SAMEDIM" if they come from the same gauge configurations). Otherwise "BS_FIX" (Bootstrapping) works works for everything but is not deterministic ("BS_SAMEDIM" is similar to JK_SAMEDIM but with Bootstrapping). 
 """

import sys
import numpy as np
import random

def resampling(orig_sample, sampling_args):                # orig_sample = [observable][mont], returns [resample][observable], [observable] means different observable
    """
    Takes a list (orig_sample) of shape [num_observable][T_mont] and resamples it arcording to the sampling args. The result is a list (Resamples) of shape [num_resampling][num_observables]
    """    
    if sampling_args[0] == "JK":
        if len(orig_sample) > 1:
            sys.exit("resampling: Number of Corr can not be larger than 1 for JK, use JK_SAMEDIM instead!")
        return resampling_JK(orig_sample, sampling_args)
    elif sampling_args[0] == "JK_SAMEDIM":
        for i in range(len(orig_sample)):
            if (len(orig_sample[i]) != len(orig_sample[0])):
                sys.exit("Samples for JK_SAMEDIM are not of same DIM, use Bootstrap instead")
        return resampling_JK_SAMEDIM(orig_sample, sampling_args)
    elif sampling_args[0] == "BS_SAMEDIM":
        for i in range(len(orig_sample)):
            if len(orig_sample[i]) is not len(orig_sample[0]):
                sys.exit("Samples for BS_SAMEDIM are not of same DIM, use BS instead")
        return resampling_BS_SAMEDIM(orig_sample, sampling_args)
    elif sampling_args[0] == "BS_DIFFDIM":
        return resampling_BS_DIFFDIM(orig_sample, sampling_args)
    elif sampling_args[0] == "GAUSSIAN":
        for i in range(len(orig_sample)):
            if len(orig_sample[i]) != 1:
                sys.exit("Sample for GAUSSIAN is not 1 dimensional")
        return resampling_GAUSSIAN(orig_sample, sampling_args)
    elif sampling_args[0] == "None":
        return [np.mean(orig_sample,axis=1),]
    else:
        sys.exit("No valid sampling method given!")



def resampling_JK(orig_sample, sampling_args):
    """
    Resamples an orig_sample of one Observable with the cut-1 Jackknife method
    """
    Resamples = np.zeros((len(orig_sample[0]), 1))
    num_JK = len(Resamples)
    for i in range(num_JK):
        for j in range(num_JK):
            if i != j:
                Resamples[i][0] += orig_sample[0][i]/num_JK
    return Resamples

def resampling_JK_SAMEDIM(orig_sample, sampling_args):
    """
    Resamples an orig_sample of one Observable with the cut-1 Jackknife method with more than one Observable coming from the same gauge configurations
    """
    num_JK = len(orig_sample[0])
    num_obs = len(orig_sample)
    Resamples = np.zeros((num_JK, num_obs))
    for i in range(num_obs):
        for j in range(num_JK):
            for k in range(num_JK):
                if j != k:
                    Resamples[j][i] += orig_sample[i][j]/num_JK
    return Resamples



def resampling_BS_SAMEDIM(orig_sample, sampling_args):
    """
    Resamples an orig_sample of Observables by bootstrapping num times, 
    """
    num_BS = sampling_args[1]
    num_obs = len(orig_sample)
    num_mont = len(orig_sample[0])

    Resamples = np.zeros((num_BS, num_obs))

    for i in range(num_mont):
        for j in range(num_BS):
            randint = random.randint(0,num_mont-1)
            for k in range(num_obs):
                Resamples[j][k] += orig_sample[k][randint]/num_mont
    return Resamples

def resampling_BS_DIFFDIM(orig_sample, sampling_args):
    """
    Resamples an orig_sample of Observables by bootstrapping num times, 
    """
    num_BS = sampling_args[1]
    num_obs = len(orig_sample)
    # num_mont = len(orig_sample[0])                                # kann nicht generell definiert werden!

    Resamples = np.zeros((num_BS, num_obs))

    for i in range(num_BS):
        for j in range(num_obs):
            num_mont = len(orig_sample[j])
            for k in range(num_mont):
                randint = random.randint(0,num_mont-1)
                Resamples[i][j] += orig_sample[j][randint]/num_mont
    return Resamples

def resampling_GAUSSIAN(orig_sample, sampling_args):
    """
    Resamples an orig_sample of Observables distributing it with a gaussian distribution
    """
    std = sampling_args[1]                                  # array: one for every observable
    num_gauss = sampling_args[2]
    num_obs = len(orig_sample)

    Resamples = np.zeros((num_gauss, num_obs))
    for i in range(num_gauss):
        for j in range(num_obs):
            # Resamples[i][j] = random.normal(orig_sample[j][0], std[j], 1)
            Resamples[i][j] = np.random.normal(orig_sample[j][0], std[j], 1)
    return Resamples

