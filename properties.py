################################################################################
#
# Global Definition
#
################################################################################

import numpy as np

# Transition Probability file name
filenameT = 'Transitions/transitionprob.p'
degree_filename = 'Transitions/degree.p'

# Initial Census Data
filenameC = 'PopulationCensus/populationCensusData.p'
areaSubPrefecture_filename = 'PopulationCensus/areaSubPrefectureCensusData.p'
densitySubPrefectureCensusData = 'PopulationCensus/densitySubPrefectureCensusData.p'
polygonPointsSubPrefecture_filename = 'PopulationCensus/polygonPointsSubPrefectureCensusData.p'
subPrefectureNumbering_filename = 'PopulationCensus/subPrefectureNumberingCensusData.p'

# Number of region
dim = 255

# Total Population
total_population = 15686986

# Probability of message transimision per contact
c = np.zeros(dim)
for i in range(dim):
    c[i] = 0.8

# Radius of transmision in km
r = 100.0/1000.0

# Return Rate
return_rate = 1.0/0.5

# alphaS, alphaI, alphaR
alphaS = 1./0.5
alphaI = 1./0.5
alphaR = 1./0.5

################################################################################
#
# End of Global Definition
#
################################################################################
