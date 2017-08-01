import sys
import numpy as np
from scipy.integrate import simps
from math import pi

# convention: a_up, b_down, a_down, b_up
Lambda = 0.0
nf=4
norb=nf/2
Uval = 8.0
Jval = 2.0
crystal_field = 3.4
shifts = np.array([0.0, -crystal_field,
                   0.0, -crystal_field], dtype=complex)


Himp = -0.5*Lambda*np.array( [
        [ 0,  0, -1J,  0,  0,  1],
        [ 0,  0,  0,  1J, -1,  0],
        [1J,  0,  0,   0,  0,-1J],
        [ 0,-1J,  0,   0,-1J,  0],
        [ 0, -1,  0,  1J,  0,  0],
        [ 1,  0, 1J,   0,  0,  0]], dtype=complex)

# from ch0/dope1.6/b22/real
mu = 1.1553024
mu = 1.1603603
mu = 2.2486

#2.5*Uval-5*Jval

for flavor in xrange(nf):
    Himp[flavor,flavor] -= mu
    Himp[flavor,flavor] += shifts[flavor]

f = open('hopping.txt','w')
for iorb in xrange(nf):
    for jorb in xrange(nf):
        print>>f, iorb, jorb, Himp[iorb,jorb].real, Himp[iorb,jorb].imag
f.close()

