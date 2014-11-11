import numpy as np
#Transition Probability file name
transitionProbability = '/transitionprob.p'
degree = '/degree.p'
#Initial Census Data
populationCensusData = '/populationCensusData'
areaSubPrefectureCensusData = '/areaSubPrefectureCensusData.p'
densitySubPrefectureCensusData = '/densitySubPrefectureCensusData.p'
polygonPointsSubPrefectureCensusData = '/polygonPointsSubPrefectureCensusData.p'
subPrefectureNumbering = '/subPrefectureNumberingCensusData.p'

dim = 255# Number of region

total_population = 15686986# Total Population

# Probability of message transimision per contact
c=np.zeros(dim)
for i in range(dim):
    c[i]=0.8
    
r = 100.0/1000.0# Radius of transmision in km

gamma = 1.0/3.0# Recovery rate

return_rate = 1.0/0.5# Return Rate

# alphaS, alphaI, alphaR
alphaS = 1./0.5 
alphaI = 1./0.5 
alphaR = 1./0.5 

# muS, muI, muR
muS = 1.0/10.0 
muI = 1.0/10.0 
muR = 1.0/10.0 

deltaEI =  1.0/3.0 #EI to R

Khi = -0.5