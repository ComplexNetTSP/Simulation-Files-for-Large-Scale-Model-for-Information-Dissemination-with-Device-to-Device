import cPickle as p
import shapefile
areaSubPrefecture_filename = 'subPrefectureNumberingCensusData'

with open(areaSubPrefecture_filename,'rU') as pf_s:
  with open(areaSubPrefecture_filename + '.p','wb') as pf_d:
    A = p.load(pf_s)
    p.dump(A, pf_d)
