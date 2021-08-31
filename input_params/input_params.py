import numpy as np
import glob


eV = 1.0
keV = 1e3
GeV = 1e9
au = 27.211386245988*eV
gram = 5.60958884493318e23*GeV
cm = 5.067730214314311e13/GeV
meter = 100*cm
a0  =  5.29177e-11*meter #Bohr radius
mElectron =  511*keV
aEM  =  1.0/137.035999139
#=======================================================================
#=======================================================================
'''

Basic data structure.

'''
class element:
    def __init__(self, name ,C,Z,n_list,E_B,Z_eff,semi_full):
        self.name = name
        self.C = C
        self.Z = Z
        self.n_list = n_list
        self.E_B = E_B
        self.Z_eff = Z_eff
        self.semi_full = semi_full

    def call(self):
        return self.name, self.C , self.Z , self.n_list , self.E_B , self.Z_eff , self.semi_full





