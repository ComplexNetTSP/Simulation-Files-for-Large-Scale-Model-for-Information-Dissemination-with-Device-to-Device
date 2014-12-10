

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

#
# Imports
#
import os
import pickle as p
import numpy as np
import argparse

from tau_leap_latent import population_at_equilibrum, stoc_eqs
from progressbar import ProgressBar, Percentage, RotatingMarker, ETA, Bar

import properties

# def save_results(S, I, R, ES, EI, ER, A, directory='Results'):
#     if not os.path.exists(directory):
#         os.makedirs(directory)
#     with open(directory+'/S.p', 'wb') as fp:
#         p.dump(S, fp)
#     with open(directory+'/I.p', 'wb') as fp:
#         p.dump(I, fp)
#     with open(directory+'/R.p', 'wb') as fp:
#         p.dump(R, fp)
#     with open(directory+'/ES.p', 'wb') as fp:
#         p.dump(ES, fp)
#     with open(directory+'/EI.p', 'wb') as fp:
#         p.dump(EI, fp)
#     with open(directory+'/ER.p', 'wb') as fp:
#         p.dump(ER, fp)
#     with open(directory+'/A.p', 'wb') as fp:
#         p.dump(A, fp)

# def get_beta(densitySubPrefecture_filename,
#              polygonPointsSubPrefecture_filename,
#              subPrefectureNumbering_filename,
#              r,
#              c):
#
#     with open(densitySubPrefecture_filename, "rb") as pickleFile:
#         RhoPolygons = p.load(pickleFile)
#     with open(polygonPointsSubPrefecture_filename, "rb") as pickleFile:
#         PolygonPoints = p.load(pickleFile)
#         community = len(PolygonPoints.keys())
#         listing = PolygonPoints.keys()
#     with open(subPrefectureNumbering_filename, "rb") as pickleFile:
#         ConnectionNumber = p.load(pickleFile)
#
#     beta = np.zeros(community)
#     for i in listing:
#         k = RhoPolygons[i]*(np.pi)*r**2
#         beta[ConnectionNumber[i]-1] = -k*np.log(1-c[ConnectionNumber[i]-1])
#
#     return beta


# def initial_population(areaSubPrefecture_filename,
#                        densitySubPrefecture_filename,
#                        polygonPointsSubPrefecture_filename,
#                        subPrefectureNumbering_filename,
#                        Totalpopulation):
#
#     with open(areaSubPrefecture_filename, "rb") as pickleFile:
#         AreaPolygons = p.load(pickleFile)
#     with open(densitySubPrefecture_filename, "rb") as pickleFile:
#         RhoPolygons = p.load(pickleFile)
#     with open(polygonPointsSubPrefecture_filename, "rb") as pickleFile:
#         PolygonPoints = p.load(pickleFile)
#         community = len(PolygonPoints.keys())
#         listing = PolygonPoints.keys()
#     with open(subPrefectureNumbering_filename, "rb") as pickleFile:
#         ConnectionNumber = p.load(pickleFile)
#         #print ConnectionNumber
#     N0 = np.zeros(community)
#     for i in listing:
#         N0[ConnectionNumber[i]-1] = AreaPolygons[i]*RhoPolygons[i]/float(Totalpopulation)
#     return np.identity(community)*N0


# def get_transition_probability(filename):
#     with open(filename, "rb") as pickleFile:
#         Tlist = p.load(pickleFile)
#     #
#     # Transition Probability
#     #
#     Tarray = np.array(Tlist, dtype=np.float)
#     O = np.ones((255, 255)) - np.identity(255)
#     Tarray = Tarray * O
#
#     with np.errstate(invalid='ignore'):
#         res1 = Tarray*(1/Tarray.sum(axis=1))
#
#     for i in xrange(255):
#         for j in xrange(255):
#             if np.isnan(res1[i, j]):
#                 res1[i, j] = 0.0
#     #
#     # Per Capita Leaving Rate
#     #
#     res2 = Tarray.sum(axis=1)/Tarray.sum()
#     return res1, res2

def rate_of_return(dim, rate):
    rho = np.zeros((dim, dim))
    for i in xrange(dim):
        for j in xrange(dim):
            if i != j:
                rho[i, j] = rate
    return rho


# def compute_population_at_equilibrium(N0, dim, sigma, nu, rho, total_population):
#     N = np.zeros((dim, dim))
#     for i in range(dim):
#         Ni = np.sum(N0.reshape(dim, dim), axis=1)
#         N[i, :] = np.floor(population_at_equilibrum(i, sigma, nu, rho, Ni[i])*total_population)
#     return N

#
# MAIN FUNCTION THAT RUN THE SIMULATION
#
def run_simumation(N0, dim, tau, beta, sigma, nu, rho, gamma, total_population, simulation_end_time,alphaS,alphaI,alphaR,muS,muI,muR,deltaEI,initialInfectedCommunity):
    # Steps
    steps = int(simulation_end_time*(1.0/tau))
    # Compute the initial population distribution
    N = compute_population_at_equilibrium(N0, dim, sigma, nu, rho, total_population)
    print 'average population per cellid: ', np.sum(N, axis=0)

    #
    # init the progress bar
    #
    widgets = ['Simulation: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
           ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval=steps).start()

    #
    #
    #
    # Inititial Population in each States
    S = N.copy()
    I = np.zeros((dim, dim))
    R = np.zeros((dim, dim))
    ES = np.zeros((dim,dim))
    EI = np.zeros((dim,dim))
    ER = np.zeros((dim,dim))

    # Infect some nodes
    initital_infection = 100.0
    S[initialInfectedCommunity,initialInfectedCommunity] = S[initialInfectedCommunity, initialInfectedCommunity]-initital_infection
    I[initialInfectedCommunity, initialInfectedCommunity] = initital_infection
    # Stack the differents S.I.R. variables in one vector
    Y = S.reshape(dim*dim).tolist()
    Y = np.append(Y, I.reshape(dim*dim).tolist())
    Y = np.append(Y, R.reshape(dim*dim).tolist())
    Y = np.append(Y, ES.reshape(dim*dim).tolist())
    Y = np.append(Y, EI.reshape(dim*dim).tolist())
    Y = np.append(Y, ER.reshape(dim*dim).tolist())

    Sr = []
    Ir = []
    Rr = []
    ESr = []
    EIr = []
    ERr = []

    InfectionMatrix = np.zeros((steps, 255))
    for step in xrange(steps):
        Ytemp = stoc_eqs(Y, tau, beta, gamma, sigma, nu, rho, dim,alphaS,alphaI,alphaR,muS,muI,muR,deltaEI)
        Ytemp = Ytemp.reshape((6, dim*dim))
        Stemp = Ytemp[0].reshape((dim, dim))
        Itemp = Ytemp[1].reshape((dim, dim))
        Rtemp = Ytemp[2].reshape((dim, dim))
        Sr.append(Stemp.sum())
        Ir.append(Itemp.sum())
        Rr.append(Rtemp.sum())
        EStemp = Ytemp[3].reshape((dim,dim))
        EItemp = Ytemp[4].reshape((dim,dim))
        ERtemp = Ytemp[5].reshape((dim,dim))
        ESr.append(EStemp.sum())
        EIr.append(EItemp.sum())
        ERr.append(ERtemp.sum())
        InfectionMatrix[step, :] = Itemp.sum(axis=0)
        Y = Ytemp
        pbar.update(step)
    pbar.finish()
    return Sr, Ir, Rr, ESr, EIr, ERr, InfectionMatrix

if __name__ == '__main__':
  #
  # Parse argument
  #

  parser = argparse.ArgumentParser(description='Process SIR simulation with latent states.')
  parser.add_argument('--output', help='output directory', required=True)
  parser.add_argument('--duration', type=int, help='simulation duration in days', required=True)
  parser.add_argument('--tau', type=float, help='simulation step (fraction of day)', default=1.0/5)
  parser.add_argument('--mu', type=float, help='simulation mu for latent state (fraction of the population)', default=1.0/10)
  parser.add_argument('--sim-id', type=int, help='simulation step (fraction of day)', default=1.0/5)
  parser.add_argument('--cell-id', type=int, help='initial cellID', default=0)
  parser.add_argument('--gamma', type=float, help='recovery rate', default=1.0/3.0)

  args = parser.parse_args()
  # Simualtion parameters
  simulation_end_time = float(args.duration)
  # Simulation Step
  tau = float(args.tau)
  muS = float(args.mu)
  muI = float(args.mu)
  muR = float(args.mu)
  gamma = float(args.gamma)

  # EI to R
  deltaEI = gamma

  simulation_id=int(args.sim_id)
  cell_id = args.cell_id


  argsdict = vars(args)

  if( args.output and
      args.mu and
      args.tau and
      args.duration and
      args.mu and
      args.sim_id
    ):

    output_dir = argsdict['output']

    if output_dir.endswith('\\'):
      output_dir = output_dir[:-1]

    # if output dire doesn' extist create it
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    #
    # Start Simulation
    #
    beta = get_beta(
        properties.densitySubPrefectureCensusData,
        properties.polygonPointsSubPrefectureCensusData,
        properties.subPrefectureNumbering,
        properties.r,
        properties.c)

    with np.errstate(divide='ignore'):
        (nu, sigma) = get_transition_probability(properties.transitionProbability)

    rho = rate_of_return(properties.dim, properties.return_rate)

    N0 = initial_population(
        properties.areaSubPrefectureCensusData,
        properties.densitySubPrefectureCensusData,
        properties.polygonPointsSubPrefectureCensusData,
        properties.subPrefectureNumbering,
        properties.total_population)

    #
    # Simulation community=0
    #
    S,I,R,ES,EI,ER,InfectionMatrix = run_simumation(N0,
                                                    properties.dim,
                                                    tau,
                                                    beta,
                                                    sigma,
                                                    nu,
                                                    rho,
                                                    gamma,
                                                    properties.total_population,
                                                    simulation_end_time,
                                                    properties.alphaS,
                                                    properties.alphaI,
                                                    properties.alphaR,
                                                    muS,
                                                    muI,
                                                    muR,
                                                    deltaEI,
                                                    cell_id)
    A = InfectionMatrix.T
    save_results(S, I, R, ES, EI, ER, A, output_dir + '/' + str(simulation_id))

    #####################
    #
    # end Simulation
    #
  else:
    parser.print_help()
