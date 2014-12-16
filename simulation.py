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
import numpy as np
import argparse
from progressbar import ProgressBar, Percentage, RotatingMarker, ETA, Bar

#
# Custom imports
#
from tau_leap import stoc_eqs
from utils import *
import properties

###############################################################################
#
# Function rate_of_return
#
###############################################################################


def rate_of_return(dim, rate):
    rho = np.zeros((dim, dim))
    for i in xrange(dim):
        for j in xrange(dim):
            if i != j:
                rho[i, j] = rate
    return rho

###############################################################################
#
# MAIN FUNCTION THAT RUN THE SIMULATION
#
###############################################################################


def run_simumation(N0, dim, tau, beta, sigma, nu, rho, total_population, simulation_end_time, initialInfectedCommunity):
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

    # Infect some nodes
    initital_infection = 100.0
    S[initialInfectedCommunity, initialInfectedCommunity] = S[
        initialInfectedCommunity, initialInfectedCommunity] - initital_infection
    I[initialInfectedCommunity, initialInfectedCommunity] = initital_infection
    # Stack the differents S.I.R. variables in one vector
    Y = S.reshape(dim * dim).tolist()
    Y = np.append(Y, I.reshape(dim * dim).tolist())
    Y = np.append(Y, R.reshape(dim * dim).tolist())
    Sr = []
    Ir = []
    Rr = []
    InfectionMatrix = np.zeros((steps, 255))
    for step in xrange(steps):
        Ytemp = stoc_eqs(Y, tau, beta, gamma, sigma, nu, rho, dim)
        Ytemp = Ytemp.reshape((3, dim * dim))
        Stemp = Ytemp[0].reshape((dim, dim))
        Itemp = Ytemp[1].reshape((dim, dim))
        Rtemp = Ytemp[2].reshape((dim, dim))
        Sr.append(Stemp.sum())
        Ir.append(Itemp.sum())
        Rr.append(Rtemp.sum())
        InfectionMatrix[step, :] = Itemp.sum(axis=0)
        Y = Ytemp
        pbar.update(step)
    pbar.finish()
    return Sr, Ir, Rr, InfectionMatrix

if __name__ == '__main__':
    #
    # Parse argument
    #
    parser = argparse.ArgumentParser(
        description='Process SIR simulation with nolatent states.')
    parser.add_argument('--output', help='output directory', required=True)
    parser.add_argument(
        '--duration', type=int, help='simulation duration in days', required=True)
    parser.add_argument(
        '--tau', type=float, help='simulation step (fraction of day)', default=1.0 / 5)
    parser.add_argument(
        '--sim-id', type=int, help='simulation step (fraction of day)', default=1.0 / 5)
    parser.add_argument(
        '--cell-id', type=int, help='initial cellID', default=0)
    parser.add_argument(
        '--gamma', type=float, help='recovery rate', default=1.0 / 3.0)

    args = parser.parse_args()

    # Simualtion parameters
    simulation_end_time = float(args.duration)
    tau = float(args.tau)
    simulation_id = int(args.sim_id)
    cell_id = int(args.cell_id)
    gamma = float(args.gamma)
    cell_id = args.cell_id

    argsdict = vars(args)

    conditions_mets = (
        args.output and
        args.tau and
        args.duration and
        args.sim_id)

    if conditions_mets:

        output_dir = argsdict['output']
        if output_dir.endswith('\\'):
            output_dir = output_dir[:-1]

        # if output dire doesn' extist create it
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        #
        # Start Simulation
        #
        beta = get_beta(properties.densitySubPrefectureCensusData,
                        properties.polygonPointsSubPrefectureCensusData,
                        properties.subPrefectureNumbering,
                        properties.r,
                        properties.c)

        with np.errstate(divide='ignore'):
            (nu, sigma) = get_transition_probability(
                properties.transitionProbability)

        rho = rate_of_return(properties.dim, properties.return_rate)

        N0 = initial_population(properties.areaSubPrefectureCensusData,
                                properties.densitySubPrefectureCensusData,
                                properties.polygonPointsSubPrefectureCensusData,
                                properties.subPrefectureNumbering,
                                properties.total_population)

        #
        # Simulation
        #
        S, I, R, InfectionMatrix = run_simumation(N0,
                                                  properties.dim,
                                                  tau,
                                                  beta,
                                                  sigma,
                                                  nu,
                                                  rho,
                                                  properties.total_population,
                                                  simulation_end_time,
                                                  cell_id
                                                  )
        A = InfectionMatrix.T
        save_results(S, I, R, A, output_dir + '/' + str(simulation_id))

        #####################
        #
        # end Simulation
        #
    else:
        parser.print_help()
