"""
    @Author: Yannick Dengler
    @Date:   2023-Sep-8
    @Last Modified by: Yannick Dengler
    
    Takes data and parses it to the resampler and the functions to get an estimation of the error
"""
import sys
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import src_py.resampling as rs

def mean_OG_sample(OG_Sample):
    mean = []
    for Corr_mont in OG_Sample:
        mean.append(np.mean(Corr_mont, axis=0))
    return mean

class measurement:
    def __init__(self, name, measure_func, sampling_args = ("JK", 0, 0)):
        self.name = name
        self.measure_func = measure_func
        self.sampling_args = sampling_args
        self.results = {}
    def measure(self, OG_Sample, args, check_for_bad_results = True):
        mean_res = self.measure_func(mean_OG_sample(OG_Sample), args)     
        for res_key, res in mean_res.items():
            # print(res_key, res)
            for val in res:
                if np.isnan(val) or np.isinf(val):
                    sys.exit("NaN or inf in Mean: %s"%(res_key))
        result_samples = {}                                                             
        for key in mean_res:
            result_samples[key] = []
        Resamples = rs.resampling(OG_Sample, self.sampling_args)
        tot = len(Resamples)                            
        for i, Resample in zip(range(len(Resamples)), Resamples):
            print(i*100./tot, "%")
            if check_for_bad_results:
                bad_res = True
                while bad_res:
                    res_tmp = self.measure_func(Resample, args)                               
                    for res in res_tmp.values():
                        if not(any(np.isnan(res)) or any(np.isinf(res))):
                            bad_res = False
            else:
                res_tmp = self.measure_func(Resample, args)                                   
            for key in res_tmp:
                result_samples[key].append(res_tmp[key])
        self.result_names = []
        for key in result_samples:
            self.result_names.append(key)
            self.results[key] = result(key, mean_res[key], result_samples[key], self.sampling_args)
    def visit_print(self):
        with h5py.File("../output/HDF5_logfiles"+self.name+".hdf5","r") as f:
            f.visit(print)
    def print_to_HDF(self, hdfpath = "../output/HDF5_resultfiles/"):
        os.makedirs(hdfpath, exist_ok=True)
        with h5py.File(hdfpath+self.name+".hdf5","w") as f:
            f.create_dataset(str(self.name), data = self.name)
            f.create_dataset("result_names", data = self.result_names)
            for i in range(len(self.sampling_args)):
                f.create_dataset("sampling_args"+str(i), data = self.sampling_args[i])
            for result in self.results.values():
                f.create_dataset(str(result.name)+"_sample", data = result.sample)
                f.create_dataset(str(result.name)+"_result", data = result.result)
            f.visit(print)
    def read_from_HDF(self, hdfpath = "../output/HDF5_resultfiles/"):
        with h5py.File(hdfpath+self.name+".hdf5","r") as f:
            f.visit(print)
            self.result_names = []
            for name in f["result_names"]:
                self.result_names.append(name.decode())
            for result_name in self.result_names:
                self.results[result_name] = result(result_name, np.asarray(f[result_name+"_result"]), np.asarray(f[result_name+"_sample"]), self.sampling_args)
    def print_everything(self):
        print(self.name)
        for result in self.results.values():
            print(result.name, result.median, result.e)

def JK_err(Sample, res):                                                               # Sample [sample][results]  # res [result]
    size = len(Sample)
    err = np.zeros(len(Sample[0]))                                                                                   # [result]
    tmp = np.swapaxes(Sample, 0, 1)
    for i in range(len(tmp)):
        for val in tmp[i]:
            err[i] += (res[i]-val)**2
        err[i] = np.sqrt(err[i]*(size-1)/size)
    return err

def BS_err(Sample):                                                                                 # Sample [sample][results]
    size = len(Sample)
    mean = np.mean(Sample, axis = 0)
    err = np.zeros(len(Sample[0]))                                                                                   # [result]
    tmp = np.swapaxes(Sample, 0, 1)
    for i in range(len(tmp)):
        for val in tmp[i]:
            err[i] += (mean[i]-val)**2
        err[i] = np.sqrt(err[i]/size)
    return err

class result:
    def __init__(self, name, result, res_sample, sampling_args = ("JK", 0, 0)):
        if not sampling_args[0] == "None":
            self.name = name                                          
            self.sampling_args = sampling_args
            self.result = result
            self.sample = res_sample                                            # [num_results]
            self.mean = np.mean(res_sample, axis=1)
            self.median = []
            self.e = []
            self.ep = []
            self.em = []
            for tmp in np.swapaxes(res_sample,0,1):
                num = len(tmp)                                                  # num_results
                tmp_2 = np.sort(tmp, kind="mergesort")
                median = tmp_2[num//2]
                self.median.append(median)
                ep_tmp = abs(tmp_2[math.ceil(num*2/3.)]-median)
                em_tmp = abs(tmp_2[math.floor(num/3.)]-median)
                self.ep.append(ep_tmp)
                self.em.append(em_tmp)
                self.e.append(max(ep_tmp,em_tmp))
            if sampling_args[0] == "JK" or sampling_args[0] == "JK_SAMEDIM":
                self.e_JK = JK_err(res_sample, result)
            elif sampling_args[0] == "BS" or sampling_args[0] == "BS_SAMEDIM" or sampling_args[0] == "BS_FIX":
                self.e_BS = BS_err(res_sample)

    def draw_histogram(self, ind=0, num_bins = 20):
        data = np.swapaxes(self.sample,0,1)[ind]
        start = 0
        stop = len(data)
        fig = plt.figure(figsize=(6, 6))
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                            left=0.1, right=0.9, bottom=0.1, top=0.9,
                            wspace=0.05, hspace=0.05)
        ax = fig.add_subplot(gs[1, 0])
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
        bins = np.linspace(min(data),max(data), num_bins)
        x_axis = []
        for i in range(len(data)):
            x_axis.append(start + i)
        ax.plot(x_axis, data)
        plt.hist(data,bins = bins, orientation="horizontal")
        plt.show()
    def percent_nan(self):
        nans = 0
        tot = len(self.sample)
        for res in self.sample:
            if any(np.isnan(res)):
                nans += 1
        print(nans*100./tot,"%")
