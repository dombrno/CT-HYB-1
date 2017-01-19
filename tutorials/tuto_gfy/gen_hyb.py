import os
import numpy as np
import matplotlib.pyplot as plt
from numpy import array
from pylab import figure, xticks, yticks
import pylab
import socket
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from math import pi
import sys
from numpy import array
import numpy as np
from math import pi
from matplotlib.colors import LogNorm
from matplotlib.ticker import AutoMinorLocator
import h5py as h5
from pathlib import Path, PurePath

target_location = PurePath(os.environ['HOME'], 'Code', 'Linear_response', 'Source', 'Misc')

print Path.cwd()

if str(target_location) not in sys.path:
    sys.path.append(str(target_location))
print sys.path

import analyze_output as an

data_root = PurePath(os.environ['HOME'], 'Code', 'Linear_response', 'Data',
                     'Alps2', 'lr0.4118', '1site', 'Delta_3.6',
                     'ch0', 'dope1.6', 'b22')

all_orbitals = ['a_up','b_down','a_down','b_up']
nb_sites = 1
nb_orbital_couples = len(all_orbitals)**2
nb_orbitals = len(all_orbitals * nb_sites)
n_orbitals = len(all_orbitals)


configs = ['real']

print "data root:\n", data_root
print
print "configs: \n", configs



config_index = -1
beta = 60.0
config = configs[config_index]
current_path = PurePath(data_root)

gtau2 = {}
print str(Path(current_path, 'c_delta_0.h5'))

print "trying to open", current_path

corresp = {}
corresp[0] = 0
corresp[1] = 1
corresp[2] = 2
corresp[3] = 3


with h5.File(str(Path(current_path, 'c_delta_0.h5')), 'r') as input_file:
    for i in range(n_orbitals * n_orbitals):
        row_index = i / n_orbitals
        col_index = i % n_orbitals
        target_h5_name = '/'
        my_group = input_file.get(target_h5_name)
        temp_data = (np.array(my_group.get('Delta_' + str(i))[:, 0] +
                        my_group.get('Delta_' + str(i))[:, 1] * 1j))
        # Minus sign: the od_hyb framework has this minus
        # inherited from Fortran
        gtau2[i] = -temp_data

output_file = str(Path(Path.cwd(), 'ec_delta.txt'))
with open(output_file, 'w') as output_file:
    for tau_index in np.arange(gtau2[i].shape[0]):
        for i in np.arange(n_orbitals):
            for j in np.arange(n_orbitals):
                orb_index = corresp[i] * n_orbitals + corresp[j]
                out_str = (str(tau_index) + '  ' +
                           str(i) + '  ' +
                           str(j) + '  ' +
                           str(np.real(gtau2[orb_index][tau_index])) + '  ' +
                           str(np.imag(gtau2[orb_index][tau_index])) + '\n')
                output_file.write(out_str)
