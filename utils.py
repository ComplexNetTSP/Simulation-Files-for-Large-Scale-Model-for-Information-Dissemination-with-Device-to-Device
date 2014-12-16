import pickle as p
import numpy as np
import os

#
# Custom function import
#
from tau_leap_latent import population_at_equilibrum

###############################################################################
#
# Function Get beta
#
###############################################################################


def get_beta(densitySubPrefecture_filename,
             polygonPointsSubPrefecture_filename,
             subPrefectureNumbering_filename,
             r,
             c):
    '''
    Compute the beta parameter
    '''
    with open(densitySubPrefecture_filename, "rb") as pickleFile:
        RhoPolygons = p.load(pickleFile)
    with open(polygonPointsSubPrefecture_filename, "rb") as pickleFile:
        PolygonPoints = p.load(pickleFile)
        community = len(PolygonPoints.keys())
        listing = PolygonPoints.keys()
    with open(subPrefectureNumbering_filename, "rb") as pickleFile:
        ConnectionNumber = p.load(pickleFile)

    beta = np.zeros(community)
    for i in listing:
        k = RhoPolygons[i]*(np.pi)*r**2
        beta[ConnectionNumber[i]-1] = -k*np.log(1-c[ConnectionNumber[i]-1])

    return beta

###############################################################################
#
# Function save_results
#
###############################################################################


def save_results(S, I, R, A, ES=None, EI=None, ER=None, directory='Results'):
    """
    Save the result of the number of people in each different statesin in separate files
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(directory + '/S.p', 'wb') as fp:
        p.dump(S, fp)
    with open(directory + '/I.p', 'wb') as fp:
        p.dump(I, fp)
    with open(directory + '/R.p', 'wb') as fp:
        p.dump(R, fp)
    if ES:
        with open(directory + '/ES.p', 'wb') as fp:
            p.dump(ES, fp)
    if EI:
        with open(directory + '/EI.p', 'wb') as fp:
            p.dump(EI, fp)
    if ER:
        with open(directory + '/ER.p', 'wb') as fp:
            p.dump(ER, fp)
    with open(directory + '/A.p', 'wb') as fp:
        p.dump(A, fp)


##############################################################################
#
# Function get_transition_probability
#
##############################################################################


def get_transition_probability(filename):
    with open(filename, "rb") as pickleFile:
        Tlist = p.load(pickleFile)
    #
    # Transition Probability
    #
    Tarray = np.array(Tlist, dtype=np.float)
    O = np.ones((255, 255)) - np.identity(255)
    Tarray = Tarray * O

    with np.errstate(invalid='ignore'):
        res1 = Tarray * (1 / Tarray.sum(axis=1))

    for i in xrange(255):
        for j in xrange(255):
            if np.isnan(res1[i, j]):
                res1[i, j] = 0.0
    #
    # Per Capita Leaving Rate
    #
    res2 = Tarray.sum(axis=1) / Tarray.sum()
    return res1, res2


###############################################################################
#
# Function compute_population_at_equilibrium
#
###############################################################################


def compute_population_at_equilibrium(
        N0,
        dim,
        sigma,
        nu,
        rho,
        total_population):

    N = np.zeros((dim, dim))
    for i in range(dim):
        Ni = np.sum(N0.reshape(dim, dim), axis=1)
        N[i, :] = np.floor(
            population_at_equilibrum(i, sigma, nu, rho, Ni[i]) * total_population)
    return N
