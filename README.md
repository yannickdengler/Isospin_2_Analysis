# Scattering Analysis

This is a markdownfile

## Logfile Parsing

The code "HDF5.py" in "Parsing" extracts the necessary informations from the logfile from HiRep scattering code and parses it to a HDF5 file for simpler use later on. The Code extracts the following information:

- "logfile name": The name of the parsed logfile 
- "isospin_chanel": The isospin chanel input used for the scattering code
- "N_hits": The number of used sources
- "montecarlotimes": A list of the analysed montecarlotimes. Not every gauge configuration is saved/analysed to assure statistically correlators
- "filenames": A list of the analysed gauge configuaration filenames
- "plaquette": A list of the avarage plaquette. Can be used to check for thermalization. The length is the number of montecarlo times
- "gauge_group": The gauge group
- "beta": The inverse bare coupling
- "m_1": The mass of the first fundamental fermion
- "m_2": The mass of the second fundamental fermion
- "N_L": The spatial extend of the lattice
- "N_T": The temporal extend of the lattice
- "operators": A list of the operators extracted from the logfile. Each operator has an entry for the real and the imaginary part. The set of operators depends on the isospin channel and includes a first analysis of the correlators. For the isospin-2 channel the Pion-, Rho- and Two-Pion correlation function is calculated and normalized from the other entries.
- "correlators": A array that contains the correlation function of each measured operator. The array has 4 indices and has the size: num_Operators x num_soruces x num_Montecarlotimes x N_T.