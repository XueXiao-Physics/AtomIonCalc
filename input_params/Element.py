import numpy as np
import glob
au = 27.211386245988
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





