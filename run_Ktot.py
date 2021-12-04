from methods import *
import multiprocess 
import sys
sys.path.append('input_params')
import os
import glob
import h5py
import scipy.interpolate as si

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


# m^2 + k^2 = (T + m)^2
T2k = lambda T,m : np.sqrt( (T+m)**2 - m**2 ) 

k2T = lambda k,m : np.sqrt( k**2 + m**2) - m

p2v = lambda p,m : p/np.sqrt(m**2 + p**2)

'''
def ARinterp(ER,AR,ER_new):       
    AR_new = np.ndarray((len(ER_new),AR.shape[1]))    
    for it in range( AR.shape[1] ):    
        ar = AR[:,it]
        ar_new = np.interp( ER_new , ER , ar ,left=0.,right=0.)
        AR_new[:,it] = ar_new       
    return AR_new
'''

'''
Interface 
'''

print(' > Welcome. Run the code so you can get the atomic response function K(ER,q).')

# Choose element
element_param_files = glob.glob('input_params/Param_*.py')
element_names = [elem_file.split('/')[-1][6:-3] for elem_file in element_param_files]
arg = input('\n >>>>> Please type in the element you choose... \n > element : ')


if arg not in element_names:
    print(' \n > The element you ask for is not available.')
    sys.exit()
else:
    elem = __import__('Param_'+arg)
    name,C,Z,n_list,E_B,Z_eff,semi_full = elem.elem.call()
    datadir = 'OUTPUT/'+name
    files = glob.glob(datadir+'/*')
    files = [f for f in files if not os.path.isdir(f)] 

    
    
arg2 = input('\n >>>>> Please type in the orbital information ... eg. 4s or 3 or press enter directly to represent all (available) oribial...\n > orbital info: ')
files =  [f for f in files if arg2 in f.split('.')[0]]

if len(files)==0:
    print(' \n > No available orbital.')
    sys.exit()
              


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
N_ER = 2**10
ER_new = np.logspace(1, 5, N_ER)
N_q = 2**9
q_new = np.logspace(3 , 6 , N_q)
print('\n > ER is fixed between ',[ER_new.min(),ER_new.max()],' N = ',N_ER,' log spaced.')
input('\n >>>>> Press Enter to run... <<<<<')


K = 0
for i in range(len(ARdatasets)):
    ar = ARdatasets[i]
    n,l = QNs[i]
    ER = k2T(ar.kPrime_grid,mElectron) + abs(E_B[n-1][l])
    ifun = si.RectBivariateSpline(np.log10(ER),np.log10(ar.q_grid),\
                                  ar.Atomic_Response_K[0])
    mask1 = (ER_new>ER.min())*(ER_new<ER.max())
    mask2 = (q_new>ar.q_grid.min())*(q_new<ar.q_grid.max())
    K += ifun(np.log10(ER_new),np.log10(q_new)) * mask1[:,None]*mask2[None,:]

try:
    os.mkdir('OUTPUT2/')
except:
    pass
    
filename = 'OUTPUT2/'+name+arg2+'_Ktot.hdf5'

f = h5py.File(filename,'w')
f.create_dataset('Ktot',data=K)
f.create_dataset('ER',data=ER_new)
f.create_dataset('q',data=q_new)
f.close()

print('\n > The result is saved in',filename)
input('\n >>>>> Done. Press Enter to exit... <<<<<')

