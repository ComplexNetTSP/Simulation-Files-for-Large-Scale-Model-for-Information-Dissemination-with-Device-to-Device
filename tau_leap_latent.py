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

def stoc_eqs(Y, tau, beta, gamma, sigma, nu, rho, dim,alphaS,alphaI,alphaR,muS,muI,muR,deltaEI):
    '''
    Definition:
    - beta: infection rate
    - gamma: Recovery rate
    - sigma[i] as the per capita living rate of the region i
    - nu[i,j] the conditional probability of leaving region i for region j
    - rho[i,j] rate at which persons from region i who are visiting region j return to region i
    - pop: total population living in the cell i
    '''
    Y = Y.reshape((6,dim*dim))
    S = Y[0].reshape((dim,dim))
    I = Y[1].reshape((dim,dim))
    R = Y[2].reshape((dim,dim))
    ES = Y[3].reshape((dim,dim))
    EI = Y[4].reshape((dim,dim))
    ER = Y[5].reshape((dim,dim))
    Ni = np.sum(S, axis=0) + np.sum(I, axis=0) + np.sum(R, axis=0) + np.sum(ES, axis=0) + np.sum(EI, axis=0) + np.sum(ER, axis=0)
    Sy = S.copy()
    Iy = I.copy()
    Ry = R.copy()
    ESy = ES.copy()
    EIy = EI.copy()
    ERy = ER.copy()
    #
    # Compute the mobilty transitions
    #

    # for each community
    Rate = np.zeros((dim,dim,21))
    for i in xrange(dim):
      for j in range(dim):
        if i != j:
          # Suceptible departing from i to j
          #print i,j,(rho[i,j]*Sy[i,j])*tau
          Rate[i,j,0] = min(poisson((rho[i,j]*Sy[i,j])*tau), Sy[i,j])
          Sy[i,j] -= Rate[i,j,0]
          Sy[i,i] += Rate[i,j,0]

          # Infected departing from i to j
          Rate[i,j,2] = min(poisson((rho[i,j]*Iy[i,j])*tau), Iy[i,j])
          Iy[i,j] -= Rate[i,j,2]
          Iy[i,i] += Rate[i,j,2]

          # Recovered departing from i to j
          Rate[i,j,4] = min(poisson((rho[i,j]*Ry[i,j])*tau), Ry[i,j])
          Ry[i,j] -= Rate[i,j,4]
          Ry[i,i] += Rate[i,j,4]

          # E_Suceptible departing from i to j
          Rate[i,j,8] = min(poisson((rho[i,j]*ESy[i,j])*tau), ESy[i,j])
          ESy[i,j] -= Rate[i,j,8]
          ESy[i,i] += Rate[i,j,8]

          # EInfected departing from i to j
          Rate[i,j,9] = min(poisson((rho[i,j]*EIy[i,j])*tau), EIy[i,j])
          EIy[i,j] -= Rate[i,j,9]
          EIy[i,i] += Rate[i,j,9]

          # ERecovered departing from i to j
          Rate[i,j,10] = min(poisson((rho[i,j]*ERy[i,j])*tau), ERy[i,j])
          ERy[i,j] -= Rate[i,j,10]
          ERy[i,i] += Rate[i,j,10]
        else:
          # All possible outgoing
          for jj in xrange(dim):
            Rate[i,jj,5] = min(poisson((sigma[i]*nu[i,jj]*Sy[i,i])*tau), Sy[i,i])
            Sy[i,i] -= Rate[i,jj,5]
            Sy[i,jj] += Rate[i,jj,5]

            Rate[i,jj,6] = min(poisson((sigma[i]*nu[i,jj]*Iy[i,i])*tau), Iy[i,i])
            Iy[i,i] -= Rate[i,jj,6]
            Iy[i,jj] += Rate[i,jj,6]

            Rate[i,jj,7] = min(poisson((sigma[i]*nu[i,jj]*Ry[i,i])*tau), Ry[i,i])
            Ry[i,i] -= Rate[i,jj,7]
            Ry[i,jj] += Rate[i,jj,7]

            Rate[i,jj,11] = min(poisson((sigma[i]*nu[i,jj]*ESy[i,i])*tau), ESy[i,i])
            ESy[i,i] -= Rate[i,jj,11]
            ESy[i,jj] += Rate[i,jj,11]

            Rate[i,jj,12] = min(poisson((sigma[i]*nu[i,jj]*EIy[i,i])*tau), EIy[i,i])
            EIy[i,i] -= Rate[i,jj,12]
            EIy[i,jj] += Rate[i,jj,12]

            Rate[i,jj,13] = min(poisson((sigma[i]*nu[i,jj]*ERy[i,i])*tau), ERy[i,i])
            ERy[i,i] -= Rate[i,jj,13]
            ERy[i,jj] += Rate[i,jj,13]

    #
    # Compute the infection dynamic
    #
    for i in xrange(dim):
      for j in xrange(dim):
        boundaryCondition_Sy = 0
        boundaryCondition_Iy = 0
        boundaryCondition_Ry = 0
        boundaryCondition_EIy = 0
        boundaryCondition_ESy = 0
        boundaryCondition_ERy = 0

        # Suceptible that become infected
        Rate[i,j,1] = min(poisson(((beta[j]/Ni[j]) * (Sy[i,j]*Iy[:,j]).sum())*tau), Sy[i,j])
        boundaryCondition_Sy -= Rate[i,j,1]
        boundaryCondition_Iy += Rate[i,j,1]

        # Suceptible that become LatentS
        boundaryCondition_Sy = max(Sy[i,j]+boundaryCondition_Sy, 0)
        Rate[i,j,14] = min(poisson(muS*Sy[i,j]*tau), boundaryCondition_Sy)
        boundaryCondition_ESy += Rate[i,j,14]

        # LatentSuceptible that become S
        boundaryCondition_ESy = max(ESy[i,j]+boundaryCondition_ESy,0)
        Rate[i,j,15] = min(poisson(alphaS*ESy[i,j]*tau), boundaryCondition_ESy)

        # Infected that recover
        boundaryCondition_Iy = max(Iy[i,j]+boundaryCondition_Iy,0)
        Rate[i,j,3] = min(poisson((gamma*Iy[i,j])*tau), boundaryCondition_Iy)
        boundaryCondition_Iy -= Rate[i,j,3]
        boundaryCondition_Ry += Rate[i,j,3]

        # Infected that become LatentI
        boundaryCondition_Iy = max(Iy[i,j]+boundaryCondition_Iy, 0)
        Rate[i,j,16] = min(poisson(muI*Iy[i,j]*tau), boundaryCondition_Iy)
        boundaryCondition_EIy += Rate[i,j,16]

        # LatentI that become I
        boundaryCondition_EIy = max(EIy[i,j]+boundaryCondition_EIy, 0)
        Rate[i,j,17] = min(poisson(alphaI*EIy[i,j]*tau), boundaryCondition_EIy)
        boundaryCondition_EIy -= Rate[i,j,17]

        # Recovered that become LatentR
        boundaryCondition_Ry = max(Ry[i,j]+boundaryCondition_Ry, 0)
        Rate[i,j,18] = min(poisson(muR*Ry[i,j]*tau), boundaryCondition_Ry)
        boundaryCondition_ERy += Rate[i,j,18]

        # LatentR that become R
        boundaryCondition_ERy = max(ERy[i,j]+boundaryCondition_ERy, 0)
        Rate[i,j,19] = min(poisson(alphaR*ERy[i,j]*tau), boundaryCondition_ERy)

        # LatentI that become R
        boundaryCondition_EIy = max(EIy[i,j]+boundaryCondition_EIy, 0)
        Rate[i,j,20] = min(poisson(deltaEI*EIy[i,j]*tau), boundaryCondition_EIy)

        # Suceptible that become infected
        if Sy[i,j] - Rate[i,j,1] < 0:
          flow = Sy[i,j]
        else:
          flow = Rate[i,j,1]
        Sy[i,j] -= flow
        Iy[i,j] += flow
        # Suceptible that become LatentS
        if Sy[i,j] - Rate[i,j,14] < 0:
          flow = Sy[i,j]
        else:
          flow = Rate[i,j,14]
        Sy[i,j] -= flow
        ESy[i,j] += flow
        # Suceptible that become LatentS
        if ESy[i,j] - Rate[i,j,15] < 0:
          flow = ESy[i,j]
        else:
          flow = Rate[i,j,15]
        ESy[i,j] -= flow
        Sy[i,j] += flow
        # Infected that recover
        if Iy[i,j] - Rate[i,j,3] < 0:
          flow = Iy[i,j]
        else:
          flow = Rate[i,j,3]
        Iy[i,j] -= flow
        Ry[i,j] += flow
        # Infected that become LatentI
        if Iy[i,j] - Rate[i,j,16] < 0:
          flow = Iy[i,j]
        else:
          flow = Rate[i,j,16]
        Iy[i,j] -= flow
        EIy[i,j] += flow
        # LatentI that become I
        if EIy[i,j] - Rate[i,j,17] < 0:
          flow = EIy[i,j]
        else:
          flow = Rate[i,j,17]
        EIy[i,j] -= flow
        Iy[i,j] += flow
        # Recovered that become LatentR
        if Ry[i,j] - Rate[i,j,18] < 0:
          flow = Ry[i,j]
        else:
          flow = Rate[i,j,18]
        Ry[i,j] -= flow
        ERy[i,j] += flow
        # LatentR that become R
        if ERy[i,j] - Rate[i,j,19] < 0:
          flow = ERy[i,j]
        else:
          flow = Rate[i,j,19]
        ERy[i,j] -= flow
        Ry[i,j] += flow
        # LatentI that become R
        if EIy[i,j] - Rate[i,j,20] < 0:
          flow = EIy[i,j]
        else:
          flow = Rate[i,j,20]
        EIy[i,j] -= flow
        Ry[i,j] += flow

    Yy = Sy.reshape(dim*dim).tolist()
    Yy = np.append(Yy, Iy.reshape(dim*dim).tolist())
    Yy = np.append(Yy, Ry.reshape(dim*dim).tolist())
    Yy = np.append(Yy, ESy.reshape(dim*dim).tolist())
    Yy = np.append(Yy, EIy.reshape(dim*dim).tolist())
    Yy = np.append(Yy, ERy.reshape(dim*dim).tolist())

    return Yy
