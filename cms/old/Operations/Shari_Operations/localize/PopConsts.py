
"""Symbolic names for the constants used throughout the code.
"""

import copy

pn_European, pn_EastAsian, pn_WestAfrican, pn_LWK, pn_MKK = 1, 4, 5, 6, 7
popName = { pn_European:'European',pn_EastAsian:'EastAsian',pn_WestAfrican:'WestAfrican' }
pop2name = copy.deepcopy( popName )
popName.update( ( str(k), v ) for k, v in list(popName.items()) )

pop2hm2name = { pn_European: 'CEU', pn_EastAsian: 'JPT+CHB', pn_WestAfrican: 'YRI' }


HM2Ages = ( 10, )
HM2Pops = ( pn_European, pn_EastAsian, pn_WestAfrican )
HM2Freqs = ( 20, 40, 60, 80, 100 )

AllAges = HM2Ages
AllPops = HM2Pops
AllFreqs = HM2Freqs

CAUSAL_POS = 500000

