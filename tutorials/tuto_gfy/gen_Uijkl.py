import numpy as np
import sys

up = 0
down = 1


def complex_to_str(z):
    return "(%e,%e)"%(z.real,z.imag)

#See Shinaoka (2015): negative sign problem paper
def generate_U_tensor_SK(n_orb, U, JH):
    U_tensor = np.zeros((n_orb,2,n_orb,2,n_orb,2,n_orb,2),dtype=complex)

    num_elem = 0
    for iorb1 in xrange(n_orb):
        for iorb2 in xrange(n_orb):
            for iorb3 in xrange(n_orb):
                for iorb4 in xrange(n_orb):
                    coeff = 0.0
                    if iorb1==iorb2 and iorb2==iorb3 and iorb3==iorb4:
                        coeff = U
                    elif iorb1==iorb4 and iorb2==iorb3 and iorb1!=iorb2:
                        # this is U', i.e. dd interaction between different orbs
                        coeff = U-2*JH
                    elif iorb1==iorb3 and iorb2==iorb4 and iorb1!=iorb2:
                        coeff = JH
                    elif iorb1==iorb2 and iorb3==iorb4 and iorb1!=iorb3:
                        print iorb1, iorb2, iorb3, iorb4
                        coeff = JH

                    for isp in xrange(2):
                        for isp2 in xrange(2):
                            U_tensor[iorb1,isp,    iorb2,isp2,    iorb3,isp2,  iorb4,isp] += coeff
                            if coeff != 0.0:
                                num_elem += 1
    print U_tensor[0, 0, 0, 1, 1, 0, 1, 1]
    return U_tensor, num_elem

def generate_U_tensor_kunes(n_orb, U, JH, Jprime):
    U_tensor = np.zeros((n_orb,2,n_orb,2,n_orb,2,n_orb,2),dtype=complex)

    num_elem = 0
    for iorb1 in xrange(n_orb):
        for iorb2 in xrange(n_orb):
            for iorb3 in xrange(n_orb):
                for iorb4 in xrange(n_orb):
                    coeff = 0.0
                    if iorb1==iorb2 and iorb2==iorb3 and iorb3==iorb4:
                        coeff = U / 2
                        for isp in xrange(2):
                            isp2 = 1 - isp
                            U_tensor[iorb1, isp,    iorb2, isp2,
                                     iorb3, isp2,  iorb4, isp] = coeff
                            if np.abs(coeff)>0:
                                num_elem += 1
                    elif iorb1==iorb4 and iorb2==iorb3 and iorb1!=iorb2:
                        # this is U', i.e. dd interaction between different orbs
                        coeff = (U - 2 * JH) / 2.0  # divide by 2 for ab neq ba
                        for isp in xrange(2):
                            for isp2 in xrange(2):
                                U_tensor[iorb1, isp,    iorb2, isp2,
                                        iorb3, isp2,  iorb4, isp] = coeff
                                if np.abs(coeff)>0:
                                    num_elem += 1                                    
                        for isp in xrange(2):
                            U_tensor[iorb1, isp,    iorb2, isp,
                                        iorb3, isp,  iorb4, isp] -= JH / 2.0
                    elif iorb1==iorb3 and iorb2==iorb4 and iorb1!=iorb2:
                        coeff = JH / 2.0  # acount for exchange b <-> a
                        for isp in xrange(2):
                            isp2 = 1 - isp
                            U_tensor[iorb1, isp, iorb2, isp2,
                                     iorb3, isp2,  iorb4, isp] = coeff
                            if np.abs(coeff)>0:
                                num_elem += 1
                    elif iorb1==iorb2 and iorb3==iorb4 and iorb1!=iorb3:
                        coeff = Jprime / 2.0
                        for isp in xrange(2):
                            isp2 = 1 - isp
                            U_tensor[iorb1, isp, iorb2, isp2,
                                     iorb3, isp2,  iorb4, isp] = coeff
                            if np.abs(coeff)>0:
                                num_elem += 1
    return U_tensor, num_elem

def generate_dd_tensor(n_orb, U, JH):
    # generates the tensor with the proper crystal field splitting,
    # Using the aup bdown adown bup = (0, 1, 2, 3) convention
    n_spins = 2
    U_prime = U - 2.0 * JH
    U_tensor = np.zeros((n_orb,2,n_orb,2,n_orb,2,n_orb,2),dtype=complex)

    num_elem = 0
    for iorb1 in xrange(n_orb):
        iorb2 = 1 - iorb1
        for isp1 in xrange(n_spins):
            isp2 = 1 - isp1
            U_tensor[iorb1, isp1, iorb2, isp2,
                     iorb1, isp1, iorb2, isp2] = -0.5 * U_prime
            U_tensor[iorb1, isp1, iorb2, isp2,
                     iorb2, isp2, iorb1, isp1 ] = 0.5 * U_prime
            
            num_elem += 2

    for iorb1 in xrange(n_orb):
        iorb2 = 1 - iorb1
        for isp1 in xrange(n_spins):
            isp2 = isp1
            U_tensor[iorb1, isp1, iorb2, isp2,
                     iorb1, isp1, iorb2, isp2] = -0.5 * U_prime + 0.5 * JH
            U_tensor[iorb1, isp1, iorb2, isp2,
                     iorb2, isp2, iorb1, isp1 ] = 0.5 * U_prime - 0.5 * JH
            num_elem += 2

    for iorb1 in xrange(n_orb):
        iorb2 = iorb1
        for isp1 in xrange(n_spins):
            isp2 = 1 - isp1
            U_tensor[iorb1, isp1, iorb2, isp2,
                     iorb1, isp1, iorb2, isp2] = -0.5 * U
            U_tensor[iorb1, isp1, iorb2, isp2,
                     iorb2, isp2, iorb1, isp1 ] = 0.5 * U
            num_elem += 2
    print "found ", num_elem, "non zero coefficients in local Hamiltonian"
    return U_tensor, num_elem



n_sites = 2
Uval = 4.7
Jval = 0.25 * Uval
Jval = 1.0
#Jprime = 0.002

V_mat = np.identity(2 * n_sites, dtype=complex)

U_tensor, num_elem = generate_dd_tensor(n_sites, Uval, Jval)

corresp = {}
corresp[(0, 0)] = 0
corresp[(0, 1)] = 2
corresp[(1, 0)] = 3
corresp[(1, 1)] = 1

f = open("Uijkl.txt", "w")
print >>f, num_elem
line = 0
for iorb1 in xrange(n_sites):
    for iorb2 in xrange(n_sites):
        for iorb3 in xrange(n_sites):
            for iorb4 in xrange(n_sites):
                for isp1 in xrange(2):
                    for isp2 in xrange(2):
                        for isp3 in xrange(2):
                            for isp4 in xrange(2):
                                if U_tensor[iorb1,isp1,iorb2,isp2,iorb3,isp3,iorb4,isp4] != 0.0:
                                    print >>f, line, "   ", corresp[(iorb1,isp1)], corresp[(iorb2,isp2)], corresp[(iorb3,isp3)], corresp[(iorb4,isp4)], U_tensor[iorb1,isp1,iorb2,isp2,iorb3,isp3,iorb4,isp4].real, U_tensor[iorb1,isp1,iorb2,isp2,iorb3,isp3,iorb4,isp4].imag
                                    line += 1

f.close()
