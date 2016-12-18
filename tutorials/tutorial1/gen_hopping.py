import sys
import numpy as np
from scipy.integrate import simps
from math import pi

#Lambda = 1.0
# convention: a_up, a_down, b_up, b_down
Lambda = 0.0
nf=4
norb=nf/2
Uval = 8.0
Jval = 0.25*Uval
shifts = np.array([0.0, 0.0, -3.4, -3.4], dtype=complex)


Himp = -0.5*Lambda*np.array( [
        [ 0,  0, -1J,  0,  0,  1],
        [ 0,  0,  0,  1J, -1,  0],
        [1J,  0,  0,   0,  0,-1J],
        [ 0,-1J,  0,   0,-1J,  0],
        [ 0, -1,  0,  1J,  0,  0],
        [ 1,  0, 1J,   0,  0,  0]], dtype=complex)

# from b60/enforce_real
mu = 1.3870789578665168

#2.5*Uval-5*Jval

for flavor in xrange(nf):
    Himp[flavor,flavor] -= mu
    Himp[flavor,flavor] += shifts[flavor]

f = open('hopping.txt','w')
for iorb in xrange(nf):
    for jorb in xrange(nf):
        print>>f, iorb, jorb, Himp[iorb,jorb].real, Himp[iorb,jorb].imag
f.close()

