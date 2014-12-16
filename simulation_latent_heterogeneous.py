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
import argparse

import pickle as p
import numpy as np

from progressbar import ProgressBar, Percentage, RotatingMarker, ETA, Bar

# import global variables
import properties

# local import
from utils import *
from tau_leap_latent import stoc_eqs

###############################################################################
#
# Function rate_of_return
#
###############################################################################


def rate_of_return(dim, rate, Khi, degree_filename):
    with open(degree_filename, 'rb') as pickleFile:
        k = p.load(pickleFile)

    kmean = np.mean(np.array(k, dtype=np.float))

    rho = np.zeros((dim, dim))
    for i in xrange(dim):
        for j in xrange(dim):
            if i != j:
                if k[j] == 0:
                    rho[i, j] = rate
                else:
                    rho[i, j] = 1.0 / \
                        ((1.0 / rate) * (k[j] ** Khi) / (kmean ** Khi))
    return rho

###############################################################################
#
# MAIN FUNCTION THAT RUN THE SIMULATION
#
###############################################################################


def run_simulation(
        N0,
        dim,
        tau,
        beta,
        sigma,
        nu,
        rho,
        total_population,
        simulation_end_time,
        alphaS,
        alphaI,
        alphaR,
        muS,
        muI,
        muR,
        deltaEI,
        initialInfectedCommunity):
    # Steps
    steps = int(simulation_end_time * (1.0 / tau))
    # Compute the initial population distribution
    N = compute_population_at_equilibrium(
        N0, dim, sigma, nu, rho, total_population)
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
    ES = np.zeros((dim, dim))
    EI = np.zeros((dim, dim))
    ER = np.zeros((dim, dim))

    # Infect some nodes
    initital_infection = 100.0
    S[initialInfectedCommunity, initialInfectedCommunity] = S[
        initialInfectedCommunity, initialInfectedCommunity] - initital_infection
    I[initialInfectedCommunity, initialInfectedCommunity] = initital_infection
    # Stack the differents S.I.R. variables in one vector
    Y = S.reshape(dim * dim).tolist()
    Y = np.append(Y, I.reshape(dim * dim).tolist())
    Y = np.append(Y, R.reshape(dim * dim).tolist())
    Y = np.append(Y, ES.reshape(dim * dim).tolist())
    Y = np.append(Y, EI.reshape(dim * dim).tolist())
    Y = np.append(Y, ER.reshape(dim * dim).tolist())

    Sr = []
    Ir = []
    Rr = []
    ESr = []
    EIr = []
    ERr = []

    InfectionMatrix = np.zeros((steps, 255))
    for step in xrange(steps):
        Ytemp = stoc_eqs(Y, tau, beta, gamma, sigma, nu, rho,
                         dim, alphaS, alphaI, alphaR, muS, muI, muR, deltaEI)
        Ytemp = Ytemp.reshape((6, dim * dim))
        Stemp = Ytemp[0].reshape((dim, dim))
        Itemp = Ytemp[1].reshape((dim, dim))
        Rtemp = Ytemp[2].reshape((dim, dim))
        Sr.append(Stemp.sum())
        Ir.append(Itemp.sum())
        Rr.append(Rtemp.sum())
        EStemp = Ytemp[3].reshape((dim, dim))
        EItemp = Ytemp[4].reshape((dim, dim))
        ERtemp = Ytemp[5].reshape((dim, dim))
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

    parser = argparse.ArgumentParser(
        description='Process SIR simulation with latent states and heterogeneous return probability')
    parser.add_argument('--output', help='output directory', required=True)
    parser.add_argument(
        '--duration', type=int, help='simulation duration in days', required=True)
    parser.add_argument(
        '--tau', type=float, help='simulation step (fraction of day)', default=1.0 / 5)
    parser.add_argument(
        '--mu', type=float, help='simulation mu for latent state (fraction of the population)', default=1.0 / 10)
    parser.add_argument(
        '--sim-id', type=int, help='simulation step (fraction of day)', default=1.0 / 5)
    parser.add_argument(
        '--cell-id', type=int, help='initial cellID', default=0)
    parser.add_argument(
        '--gamma', type=float, help='recovery rate', default=1.0 / 3.0)
    parser.add_argument(
        '--khi', type=float, help='khi recovery rate', default=-0.5)

    args = parser.parse_args()
    # Simualtion parameters
    simulation_end_time = float(args.duration)

    # Simulation Step
    tau = float(args.tau)
    muS = float(args.mu)
    muI = float(args.mu)
    muR = float(args.mu)
    gamma = float(args.gamma)
    khi = float(args.khi)
    simulation_id = int(args.sim_id)
    cell_id = args.cell_id

    # EI to R
    deltaEI = gamma

    argsdict = vars(args)

    if (args.output and
            args.mu and
            args.tau and
            args.duration and
            args.mu and
            args.sim_id):

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
            (nu, sigma) = get_transition_probability(
                properties.transitionProbability)

        rho = rate_of_return(
            properties.dim,
            properties.return_rate,
            khi,
            properties.graphDegree)

        N0 = initial_population(
            properties.areaSubPrefectureCensusData,
            properties.densitySubPrefectureCensusData,
            properties.polygonPointsSubPrefectureCensusData,
            properties.subPrefectureNumbering,
            properties.total_population)

        #
        # Simulation community=0
        #
        S, I, R, ES, EI, ER, InfectionMatrix = run_simulation(N0,
                                                              properties.dim,
                                                              tau,
                                                              beta,
                                                              sigma,
                                                              nu,
                                                              rho,
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
        save_results(
            S, I, R, ES, EI, ER, A, output_dir + '/' + str(simulation_id))

        #####################
        #
        # end Simulation
        #
    else:
        parser.print_help()
