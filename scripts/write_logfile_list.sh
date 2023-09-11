#!/usr/bin/env bash

# Generates a list of all isospin=2 scattering logfiles within the directory './input'
# This list cann subsequently be used for generating the hdf binary files with 'HDF5.py'
find ../input/  -name "out_scattering_I2" > ../input/isospin_logfiles_list
