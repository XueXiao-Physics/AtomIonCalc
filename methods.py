import numpy as np
import tqdm
import matplotlib.pyplot as plt
import scipy.special as ssp
import mpmath as mp
import scipy.integrate
import h5py
import sys
from sympy.physics.wigner import gaunt, wigner_3j

keV = 1e3
GeV = 1e9
mElectron =  511*keV
aEM  =  1.0/137.035999139
a0  =  1/aEM/mElectron #Bohr radius

'''
July 2020
xuexiao@mail.itp.ac.cn
xxueitp@gmail.com

For simplicity, this program can only process one pair of (n,l) at a time.
'''


rv2real = lambda x: np.vectorize(np.real)(x)
rv2imag = lambda x: np.vectorize(np.imag)(x)
rv2float = lambda x : np.vectorize(float)(x)
rv2log = lambda x : np.vectorize(mp.log)(x)
rv2exp = lambda x : np.vectorize(mp.exp)(x)

def Save(file_name,data_name,data,silence=False):

    f = h5py.File(file_name,'a')

    try:
        f.create_dataset(data_name , data = data)
        if silence == False:
            print('>>\''+file_name+'\':', '\''+data_name+'\'' , 'saved')
    except ValueError:
        del f[data_name]
        f.create_dataset(data_name , data = data)
        #data0 = f[data_name]
        #data0[...] = data #renew the data
        if silence == False:
            print('>>\''+file_name+'\':', '\''+data_name+'\'' , 'rewritten.')

    f.close()
    
    

'''

main

'''

class pipeline:

    def __init__(self,C,Z,n_list,E_B,Z_eff,fac,file_name):
        # Get the Whole picture
        self.global_C = C
        self.global_Z = Z
        self.global_n_list = n_list
        self.global_E_B = E_B
        self.global_Z_eff = Z_eff
        self.global_fac = fac # impose a general normalization factor, K -> fac*K
        self.create_global_quantum_number_list()
        self.file_name = file_name
        
        

        self.QN_nllL = None
        self.r_grid = None
        self.kPrime_grid = None
        self.q_grid = None
        self.R_nlr_Table = None
        self.R_final_nllkr_Table = None
        self.Spherical_Jn_Lqr_Table = None
        self.I1 = None
        self.I2 = None
        self.I3 = None
        self.Atomic_Response_W1 = None
        self.Atomic_Response_K = None


   

    def save_which(self,start=0,stop=12):

        name_list = ['QN_nllL','r_grid','kPrime_grid','q_grid',
                'R_nlr_Table','R_final_nllkr_Table','Spherical_Jn_Lqr_Table',
                'I1','I2','I3','Atomic_Response_W1','Atomic_Response_K']

        save_list = [self.QN_nllL,self.r_grid,self.kPrime_grid,self.q_grid,
                self.R_nlr_Table,self.R_final_nllkr_Table,self.Spherical_Jn_Lqr_Table,
                self.I1,self.I2,self.I3,self.Atomic_Response_W1,self.Atomic_Response_K]

        for i in range(start,stop):
            if type(save_list[i])==type(None):
                print(name_list[i], 'is not calculated')
            else:
                Save( self.file_name , name_list[i] , save_list[i] )




    # global_QN = global quantum numbers for the element.
    def create_global_quantum_number_list(self,Lmax=100):

        global_QN = np.empty((0,4))

        for n in range(1,len(self.global_C)+1):    
            for l in range(len(self.global_C[n-1])):        
                for lPrime in range(Lmax+1):                            
                    for L in range(abs(l-lPrime), l+lPrime+1):                
                        vec = np.array([n,l,lPrime,L])
                        global_QN = np.vstack([global_QN,vec])
                                        
        self.global_QN = global_QN.astype(int)




    def set_r_grid(self,rmin,rmax,Nr):

        x_grid,weight = ssp.p_roots(Nr)
        r_grid = 0.5*(rmax-rmin)*x_grid+0.5*(rmax+rmin)

        self.Nr = int(Nr)
        self.r_grid = r_grid
        self.r_weight = weight
        self.rmin = r_grid[0]
        self.rmax = r_grid[-1]
        
        #print('rmin','%.3e'%self.r_grid[0],'rmax','%.3e'%self.r_grid[-1],'1/eV')
        



    def set_kPrime_grid(self,kPrime_grid):


        self.Nk = len(kPrime_grid)
        self.kPrime_grid = kPrime_grid
        
        #print('kmin','%.3e'%self.kPrime_grid[0],'kmax','%.3e'%self.kPrime_grid[-1],'eV')
    



    def set_q_grid(self,q_grid):

        self.Nq = len(q_grid)
        self.q_grid = q_grid

        #print('qmin','%.3e'%self.q_grid[0],'qmax','%.3e'%self.q_grid[-1],'eV')





    def set_demanded_n_l(self,n,l):
    
        Save(self.file_name,'n',n)
        Save(self.file_name,'l',l)
        Save(self.file_name,'E_B',self.global_E_B[n-1][l])
        

        #f1 = np.tile(self.global_QN[:,:],(np.size(n),1,1))
        #f2 = np.tile(n,(len(self.global_QN),1)).T
        #f3 = np.tile(l,(len(self.global_QN),1)).T

        #find = np.where(  (f1[:,:,0]==f2)*(f1[:,:,1]==f3) )
        self.QN_nllL = self.global_QN[np.where( (self.global_QN[:,[0,1]]==[n,l]).all(axis=1)  )]

        
        self.QN_nl = np.unique( self.QN_nllL[:,[0,1]],axis=0)
        self.QN_nll = np.unique( self.QN_nllL[:,[0,1,2]],axis=0)
        self.QN_L = np.unique( self.QN_nllL[:,[3]],axis=0)



    def _R_fun(self,n,l,r):

        r = np.tile(r,(1,1))

        C_vec = np.array(self.global_C[n-1][l])[:,None]
        Z_vec = np.array(self.global_Z[l])[:,None]
        n_vec = np.array(self.global_n_list[l])[:,None]
        factorize_vec =  ssp.factorial(2*n_vec)
        factor_vec = np.power(2*Z_vec, n_vec + 0.5) / np.sqrt(factorize_vec)
        
        rw =  np.sum( np.power(a0,-3./2.) * C_vec * factor_vec * np.power(r/a0, n_vec-1.)\
            * np.exp(-Z_vec * r/a0) ,axis=0)
        return rw



    def calculate_initial_rwf(self):

        size = len(self.QN_nl)
        R_nlr_Table = np.ndarray(( size , self.Nr ))
        for it in range(size): 
            n,l = self.QN_nl[it]
            R_nlr_Table[it] = self._R_fun(n,l,self.r_grid)

        self.R_nlr_Table = R_nlr_Table
        

    
    def _R_final_fun(self,n,l,lPrime,kPrime,r):

        Z_eff = self.global_Z_eff
        # notice it is a complex function
        hyp1f1 = np.vectorize(mp.hyp1f1)
        fac = 1 / np.math.factorial(2*lPrime+1)
        if fac == 0:
            result = 0.
            flag = 0
        else:
            hyp1f1_table = hyp1f1(lPrime+1+1j*Z_eff[n-1][l]/kPrime/a0,(2*lPrime+2),2j*kPrime*r ) 
            part1 = lPrime*np.log(2*kPrime*r)            
            part2 = rv2log(hyp1f1_table)
            part3 = -1J*kPrime*r
            part = rv2real(rv2exp(part1+part2+part3)) # It is to avoid divergence
            
            '''
            result = 4 * np.pi * (2*kPrime*r)**lPrime\
                    * np.exp(np.pi * Z_eff[n-1][l] /2/kPrime/a0  \
                   + np.real(ssp.loggamma(lPrime+1-1j*Z_eff[n-1][l]/kPrime/a0))  )*fac \
                   * np.exp(-1j*kPrime*r) \
                   * hyp1f1(lPrime+1+1j*Z_eff[n-1][l]/kPrime/a0,(2*lPrime+2),2j*kPrime*r )
            '''
            result = 4 * np.pi * part\
                    * np.exp(np.pi * Z_eff[n-1][l] /2/kPrime/a0  \
                   + np.real(ssp.loggamma(lPrime+1-1j*Z_eff[n-1][l]/kPrime/a0))  )*fac
            flag = 1
            
            
        return result , flag



    def calculate_final_rwf(self):
               
        
        kPrime_mesh, r_mesh = np.meshgrid(self.kPrime_grid, self.r_grid, indexing='ij')

        size = len(self.QN_nll)
        R_final_nllkr_Table = np.ndarray((size,self.Nk,self.Nr))  
        for it in range(size) :           
            n,l,lPrime = self.QN_nll[it]            
            val , flag  = self._R_final_fun(n,l,lPrime,kPrime_mesh,r_mesh)
            R_final_nllkr_Table[it] = val
            self.R_final_nllkr_Table = R_final_nllkr_Table # constantly renew the data
            Save(self.file_name,'R_final_nllkr_Table',self.R_final_nllkr_Table,silence=True)
            print('>> Calculating R_final',it+1,'/',size,'n,l,l\'=',(n,l,lPrime),'\r',end='')
            
            if flag == 0:
                print('\n>> Calculation is finished for orbit','n,l=',(n,l,lPrime),'at l\'=',lPrime)
                break     
                
        if flag == 1:
            print('\n>> Calculation is finished for orbit','n,l=',(n,l,lPrime),'at l\'=',lPrime)         


                    

    def calculate_spherical_Jn_Lqr_Table(self):


        q_mesh,r_mesh = np.meshgrid(self.q_grid,self.r_grid,indexing='ij')

        size = len(self.QN_L)
        Spherical_Jn_Lqr_Table = np.ndarray((size,self.Nq,self.Nr))
        for it in range(size):       
            L = self.QN_L[it]              
            Spherical_Jn_Lqr_Table[it] = ssp.spherical_jn(L ,q_mesh*r_mesh)
        
        self.Spherical_Jn_Lqr_Table = Spherical_Jn_Lqr_Table
        
         



    # To calculate I1(q=0,kPrime) with (n,l,lPrime=l,L=0), it is to be subtracted from I1.
    def get_I1q0(self):
        size = len(self.QN_nl)
        I1_q0 = np.ndarray((size,self.Nk))
        
        L0_index = np.where(self.QN_L ==0)[0] # find L=0 in L list

        for it in range(size):
            n,l = self.QN_nl[it]
            index_nl = np.where( (self.QN_nl==np.array([n,l])[None,:]).all(axis=1) )[0]
            index_nll = np.where( (self.QN_nll==np.array([n,l,l])[None,:]).all(axis=1) )[0] # find (n,l,lPrime=l)
            # (List,k,r)
            Integrand = (self.r_grid**2)[None,None,:] * self.R_nlr_Table[index_nl,None,:]\
                    * self.R_final_nllkr_Table[index_nll,:,:] * 1. #ssp.spherical_jn(0 ,0*r_mesh) = 1
            # (List,k)        
            I1_q0[it] = 0.5*(self.rmax-self.rmin)* np.sum(self.r_weight[None,None,:]*Integrand,axis=2)
            
        self.I1_q0 = I1_q0
            
            


    '''
    Calculate Radial Integrals
    '''

    # I1 quantum numbers : (n,l,lPrime,L)
    # [List,k,q,r]
    def get_I1(self):
 
        size = len(self.QN_nllL)
        I1 = np.ndarray((size,self.Nk,self.Nq))

        for it in range(size):
            n,l,lPrime,L = self.QN_nllL[it]
            
            index_nl = np.where( (self.QN_nl[:,:]==[n,l]).all(axis=1) )[0]
            index_nll = np.where( (self.QN_nll[:,:]==[n,l,lPrime]).all(axis=1) )[0]
            index_L = np.where( (self.QN_L[:,:]==[L]).all(axis=1) )[0]            

            Integrand = (self.r_grid**2)[None,None,None,:] * self.R_nlr_Table[index_nl,None,None,:]\
                    * self.R_final_nllkr_Table[index_nll,:,None,:] * self.Spherical_Jn_Lqr_Table[index_L,None,:,:]
                         
            I1[it] = 0.5*(self.rmax-self.rmin)* np.sum(self.r_weight[None,None,None,:]*Integrand,axis=3)       

        self.I1 = I1
        
        
    def get_I2(self):
 
        size = len(self.QN_nllL)
        I2 = np.ndarray((size,self.Nk,self.Nq))

        for it in range(size):
            n,l,lPrime,L = self.QN_nllL[it]
            
            index_nl = np.where( (self.QN_nl[:,:]==[n,l]).all(axis=1) )[0]
            index_nll = np.where( (self.QN_nll[:,:]==[n,l,lPrime]).all(axis=1) )[0]
            index_L = np.where( (self.QN_L[:,:]==[L]).all(axis=1) )[0]            
            
            dRdr = np.diff( self.Spherical_Jn_Lqr_Table[:,:,:] , axis=2 )
            dRdr = np.dstack([dRdr,dRdr[:,:,-1][:,:,None]])
            Integrand = (self.r_grid**2)[None,None,None,:] * self.R_nlr_Table[index_nl,None,None,:]\
                    * self.R_final_nllkr_Table[index_nll,:,None,:] * dRdr[index_L,None,:,:]
                         
            I2[it] = 0.5*(self.rmax-self.rmin)* np.sum(self.r_weight[None,None,None,:]*Integrand,axis=3)

        self.I2 = I2
        
        
    def get_I3(self):
 
        size = len(self.QN_nllL)
        I3 = np.ndarray((size,self.Nk,self.Nq))

        for it in range(size):
            n,l,lPrime,L = self.QN_nllL[it]
            
            index_nl = np.where( (self.QN_nl[:,:]==[n,l]).all(axis=1) )[0]
            index_nll = np.where( (self.QN_nll[:,:]==[n,l,lPrime]).all(axis=1) )[0]
            index_L = np.where( (self.QN_L[:,:]==[L]).all(axis=1) )[0]            
            
            Integrand = self.r_grid[None,None,None,:] * self.R_nlr_Table[index_nl,None,None,:]\
                    * self.R_final_nllkr_Table[index_nll,:,None,:] * self.Spherical_Jn_Lqr_Table[index_L,None,:,:]
                         
            I3[it] = 0.5*(self.rmax-self.rmin)* np.sum(self.r_weight[None,None,None,:]*Integrand,axis=3)
            

        self.I3 = I3
        
        
    '''
    Atomic Response
    '''


    def get_W1_atomic_response(self):
        
        
        '''
        size = len(self.QN_nllL)
        w1_raw = np.zeros([size,self.Nk,self.Nq])
        for i in range(size):
            n,l,lPrime,L = self.QN_nllL[i]
            coeff = (2*l+1)*(2*lPrime+1)*(2*L+1)*wigner_3j(l,lPrime,L,0,0,0)**2
            w1_raw[i] = Integral[i] * coeff * 4 * self.kPrime_grid[None,:,None]**3 /( 2*np.pi )**2
        '''
            
        # sum up
        size0 = len(self.QN_nl)
        W1 = np.zeros([size0,self.Nk,self.Nq])
        for i in range(size0):
            n,l = self.QN_nl[i]
            index_nllL = np.where( ( self.QN_nllL[:,[0,1]]==[n,l]  ).all(axis=1) )[0]
            for j in index_nllL:
                n,l,lPrime,L = self.QN_nllL[j]
                coeff = (2*l+1)*(2*lPrime+1)*(2*L+1)*wigner_3j(l,lPrime,L,0,0,0)**2
                W1[i] = W1[i] + self.I1[j]**2 * coeff * 4 * self.kPrime_grid[:,None]**3 /( 2*np.pi )**3
                if L==0 and l==lPrime:
                    correction = 4*self.kPrime_grid[:,None]**3/( 2*np.pi )**3*\
                                (2*l+1)*(self.I1_q0[i][:,None]**2 - 2*self.I1_q0[i][:,None]*self.I1[j]) 
                    correction0 = np.where(correction<0 , correction , 0.)
                    W1[i] += correction0
            
            
        self.Atomic_Response_W1 = self.global_fac * W1   
        
    
        '''

        # make a general list with m and m'
        QNmm = np.empty((0,6))
    
        for n,l,lPrime,L in self.QN_nllL:
            for m in range(-l,l+1):
                for mPrime in range(-lPrime,lPrime+1):
                    add = [n,l,lPrime,L,m,mPrime]
                    QNmm = np.vstack([QNmm,add])

        QNmm = QNmm.astype(int)
        QNmm_Aux = np.unique(QNmm[:,:4],axis=0,return_inverse=True)[1]


        # calculating integral's coefficients 
        coeff_f12 = np.empty(0)
        for n,l,lPrime,L,m,mPrime in QNmm:            
            n,l,lPrime,L,m,mPrime = int(n),int(l),int(lPrime),int(L),int(m),int(mPrime)            
            add = np.sqrt(4*np.pi) * 1j**L * (-1.)**mPrime * np.sqrt(2*L+1) * \
                        float( gaunt(l,lPrime,L,m,-mPrime,0) )            
            coeff_f12 = np.append(coeff_f12,add)
        raw_f12 = coeff_f12[:,None,None] * self.Integral[QNmm_Aux,:,:]
        del coeff_f12 



        # make a shorter list in order to sum up all the L
        QN_nllmm, QN_nllmm_Aux = np.unique( QNmm[:,[0,1,2,4,5]] , axis=0 , return_inverse=True )
        QN_nl , QN_nl_Aux = np.unique( QN_nllmm[:,[0,1]] , axis=0 , return_inverse=True )
        f12_gen = lambda index : np.sum(raw_f12[np.where(QN_nllmm_Aux==index)[0],:,:],axis=0)
        f12 = np.zeros((len(QN_nllmm) ,  self.Nk   , self.Nq  ), dtype=complex)

        for it in range(len(QN_nllmm)):
            f12[it] = f12_gen(it)
        del raw_f12


        # calculating w1 without summation
        raw_w1 = np.abs(f12)**2
        del f12
        

        # calculating W1
        size = len(QN_nl)
        Atomic_Response_W1 = np.zeros([size, self.Nk, self.Nq])
        for it in range(size):
            Atomic_Response_W1[it] = np.sum( raw_w1[ np.where(QN_nl_Aux==it)[0],:,:] ,axis=0)\
             * 4 * (self.kPrime_grid**3)[None,:,None] / (2*np.pi)**3

        self.Atomic_Response_W1 = self.global_fac * Atomic_Response_W1
        '''
        
        
        







    def get_K_atomic_response(self):
        
        W1 = self.Atomic_Response_W1
        Ee_grid = np.sqrt( self.kPrime_grid**2 + mElectron**2 ) - mElectron
        
        counts = len(self.QN_nl)
        K = np.zeros(( counts ,self.Nk,self.Nq))
 
        for it in range(counts):
            n,l = self.QN_nl[it]
            K[it] = W1[it]*(aEM*mElectron)**2/4./mElectron/Ee_grid[:,None]

        self.Atomic_Response_K = K
        






    def run_all_calculations(self):
        self.save_which(0,4)
        
        self.calculate_initial_rwf()
        self.calculate_final_rwf()
        self.calculate_spherical_Jn_Lqr_Table()
        
        #self.save_which(4,6)
        
        self.get_I1q0()
        self.get_I1()
        #self.get_I2()
        #self.get_I3()
        
        #self.save_which(6,10)
        
        self.get_W1_atomic_response()
        self.get_K_atomic_response()
        
        self.save_which(10,12)




