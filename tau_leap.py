# -*- encoding: utf-8 -*-

# -------------------------------------------------------------------------------
# Copyright (c) 2014 Vincent Gauthier Telecom SudParis.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# -------------------------------------------------------------------------------

__author__ = """\n""".join(['Vincent Gauthier <vgauthier@luxbulb.org>'])

import numpy as np
from numpy.random import poisson

def division_safe(a,b):
    c = np.zeros((a.size,))
    for i in xrange(c.size):
        if a[i] != 0 and b[i] != 0.0:
            c[i] = a[i]/b[i]
        else:
            c[i] = 0.0
    return c

def population_at_equilibrum(i, sigma, nu, rho, Ni):
    '''
    Definition:
    - i: cellid
    - sigma[i] as the per capita living rate of the region i
    - nu[i,j] is the conditional probability of leaving region i for region j
    - rho[i,j] is rate at which persons in region i are returning to region j
    - Ni is the total population living in the cell i
    '''
    dim, = sigma.shape
    N = np.zeros(dim)
    for j in range(dim):
        if i != j:
            N[j] = Ni*(sigma[i] * nu[i,j]) / (rho[i,j] * (1.0 + sigma[i] * division_safe(nu[i,:],rho[i,:]).sum()))
        else:
            N[j] = Ni/(1.0 + sigma[i] * division_safe(nu[i,:],rho[i,:]).sum())
    return N

def stoc_eqs(Y, tau, beta, gamma, sigma, nu, rho, dim):
    '''
    Definition:
    - beta: infection rate
    - gamma: Recovery rate
    - sigma[i] as the per capita living rate of the region i
    - nu[i,j] the conditional probability of leaving region i for region j
    - rho[i,j] rate at which persons from region i who are visiting region j return to region i
    - pop: total population living in the cell i
    '''
    Y = Y.reshape((3,dim*dim))
    S = Y[0].reshape((dim,dim))
    I = Y[1].reshape((dim,dim))
    R = Y[2].reshape((dim,dim))
    Ni = np.sum(S, axis=0) + np.sum(I, axis=0) + np.sum(R, axis=0)
    Sy = S.copy()
    Iy = I.copy()
    Ry = R.copy()
    #
    # Compute the mobilty transitions
    #

    # for each community
    Rate = np.zeros((dim,dim,8))
    for i in xrange(dim):
        for j in range(dim):
            if i != j:
                # Suceptible departing from i to j
                #print i,j,(rho[i,j]*Sy[i,j])*tau
                Rate[i, j, 0] = min(poisson((rho[i,j]*Sy[i, j])*tau), Sy[i,j])
                # Infected departing from i to j
                Rate[i,j,2] = min(poisson((rho[i,j]*Iy[i,j])*tau), Iy[i,j])
                # Recovered departing from i to j
                Rate[i,j,4] = min(poisson((rho[i,j]*Ry[i,j])*tau), Ry[i,j])
                Sy[i, j] -= Rate[i,j,0]
                Sy[i,i] += Rate[i,j, 0]
                Iy[i,j] -= Rate[i,j,2]
                Iy[i,i] += Rate[i,j,2]
                Ry[i,j] -= Rate[i,j,4]
                Ry[i,i] += Rate[i,j,4]
            else:
                # All possible outgoing
                for jj in xrange(dim):
                    Rate[i,jj,5] = min(poisson((sigma[i]*nu[i,jj]*Sy[i,i])*tau), Sy[i,i])
                    Rate[i,jj,6] = min(poisson((sigma[i]*nu[i,jj]*Iy[i,i])*tau), Iy[i,i])
                    Rate[i,jj,7] = min(poisson((sigma[i]*nu[i,jj]*Ry[i,i])*tau), Ry[i,i])
                    Sy[i,i] -= Rate[i,jj,5]
                    Sy[i,jj] += Rate[i,jj,5]
                    Iy[i,i] -= Rate[i,jj,6]
                    Iy[i,jj] += Rate[i,jj,6]
                    Ry[i,i] -= Rate[i,jj,7]
                    Ry[i,jj] += Rate[i,jj,7]
    #
    # Compute the infection dynamic
    #
    for i in xrange(dim):
        for j in xrange(dim):
            # Suceptible becaming infected
            Rate[i,j,1] = min(poisson(((beta[i]/Ni[j]) * (Sy[i,j]*Iy[:,j]).sum())*tau), Sy[i,j])
            # Infected that recover
            Rate[i,j,3] = min(poisson(gamma*Iy[i,j]*tau), Iy[i,j])
            Sy[i,j] -= Rate[i,j,1]
            Iy[i,j] += Rate[i,j,1]
            Iy[i,j] -= Rate[i,j,3]
            Ry[i,j] += Rate[i,j,3]
    Yy = Sy.reshape(dim*dim).tolist()
    Yy = np.append(Yy, Iy.reshape(dim*dim).tolist())
    Yy = np.append(Yy, Ry.reshape(dim*dim).tolist())
    return Yy
