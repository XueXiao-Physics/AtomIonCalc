from input_params import *
from methods import *
import multiprocess 
import sys
import os
import glob
import h5py

'''

Preparations (tools)

'''
# load the file
class ARdata_class():
    def __init__(self):
        self.q_grid=None
        self.kPrime_grid=None
        self.Atomic_Response_W1=None
        self.Atomic_Response_K=None
        self.QN_nllL=None   


T2k = lambda T,m : np.sqrt( (T+m)**2 - m**2 )

k2T = lambda k,m : np.sqrt( k**2 + m**2) - m

p2v = lambda p,m : p/np.sqrt(m**2 + p**2)


def ARinterp(ER,AR,ER_new):       
    AR_new = np.ndarray((len(ER_new),AR.shape[1]))    
    for it in range( AR.shape[1] ):    
        ar = AR[:,it]
        ar_new = np.interp( ER_new , ER , ar ,left=0.,right=0.)
        AR_new[:,it] = ar_new       
    return AR_new


'''
Interface 
'''

print(' > Welcome. Please run the code after you have got the data.')

# Choose element
elements = [Xe,O,Al,Ca,Fe,Mg,Ni,Si,S]
element_dict = {e.name:e for e in elements}
arg = input('\n >>>>> Please type in the element you choose... \n > element : ')

if arg not in element_dict.keys():
    print(' \n > The element you ask for is not available.')
    sys.exit()
else:
    name,C,Z,n_list,E_B,Z_eff,semi_full = element_dict[arg].call()
    datadir = 'OUTPUT/'+name
    files = glob.glob(datadir+'/*')


# Get the Quantum Numbers
ARdatasets = []
for i in range(len(files)):
    ARdata = ARdata_class()
    f = h5py.File(files[i],'r')
    ARdata.q_grid = np.asarray(f['q_grid'])
    ARdata.kPrime_grid = np.asarray(f['kPrime_grid'])
    ARdata.Atomic_Response_W1 = np.asarray(f['Atomic_Response_W1'])
    ARdata.Atomic_Response_K = np.asarray(f['Atomic_Response_K'])
    ARdata.QN_nllL = np.asarray(f['QN_nllL'])
    ARdatasets.append(ARdata)
    f.close()

QNs = [tuple(np.unique( ar.QN_nllL[:,[0,1]] ,axis=0)[0]) for ar in ARdatasets] 
print('\n > The quantum numbers available are \n > ' + str(QNs) )
input('\n >>>>> Press Enter to continue... <<<<<')



# Setting the interpolation options
N_ER = 2**8
ER_new = np.logspace(1, 5, N_ER)
print('\n > ER is fixed between ',[ER_new.min(),ER_new.max()],' N = ',N_ER,' log spaced.')
input('\n >>>>> Press Enter to run... <<<<<')


K = 0
for i in range(len(ARdatasets)):
    ar = ARdatasets[i]
    n,l = QNs[i]
    ER = k2T(ar.kPrime_grid,mElectron) + abs(E_B[n-1][l])
    K += ARinterp(ER,ar.Atomic_Response_K[0],ER_new=ER_new)

try:
    os.mkdir('OUTPUT2/')
except:
    pass
    
filename = 'OUTPUT2/'+name+'_Ktot.hdf5'

try:
    f = h5py.File(filename,'w-')
except OSError:
    print('\n > File name exists. Fail to save.')
    input('\n >>>>> Press Enter to exit... <<<<<')
else:
    f.create_dataset('Ktot',data=K)
    f.create_dataset('ER',data=ER_new)
    f.create_dataset('q',data=ar.q_grid)
    f.close()

    print('\n > The result is saved in',filename)
    input('\n >>>>> Done. Press Enter to exit... <<<<<')

