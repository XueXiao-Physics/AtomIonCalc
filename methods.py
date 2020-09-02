#!/usr/bin/env python
# coding: utf-8

# In[1]:

from input_params import *
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import scipy.special as ssp
import mpmath as mp
import scipy.integrate
import h5py
import sys
from sympy.physics.wigner import gaunt


'''

July 2020
xuexiao@mail.itp.ac.cn
xxueitp@gmail.com

'''


rv2real = lambda x: np.vectorize(np.real)(x)
rv2float = lambda x : np.vectorize(float)(x)


def Save(file_name,data_name,data):

    f = h5py.File(file_name,'a')

    try:
        f.create_dataset(data_name , data = data)
        print('**\''+file_name+'\'' ,'info :'  , '\''+data_name+'\'' , 'data saved')
    except RuntimeError:
        print('**\''+file_name+'\'' ,'error :'  , '\''+data_name+'\'' , 'data already exists!')

    f.close()



'''

main

'''

class pipeline:

    def __init__(self,C,Z,n_list,E_B,Z_eff):

        self.C = C
        self.Z = Z
        self.n_list = n_list
        self.E_B = E_B
        self.Z_eff = Z_eff
        self._create_entire_quantum_number_list()
        self.set_demanded_n_l(None,None)
        
        


        self.QNL_nllL_select = None
        self.r_grid = None
        self.kPrime_grid = None
        self.q_grid = None
        self.R_nlr_list = None
        self.R_final_nllkr_List = None
        self.Spherical_Jn_Lqr_List = None
        self.Integral = None
        self.Atomic_Response_W1 = None
        self.Atomic_Response_K = None
        save_list = [self.QNL_nllL_select,self.r_grid]



    def set_file_name(self, file_name):

        self.file_name = file_name



    def save_which(self,start=0,stop=10):

        name_list = ['QNL_nllL_select','r_grid','kPrime_grid','q_grid',\
                'R_nlr_list','R_final_nllkr_List','Spherical_Jn_Lqr_List',
                'Integral','Atomic_Response_W1','Atomic_Response_K']

        save_list = [self.QNL_nllL_select,self.r_grid,self.kPrime_grid,self.q_grid,\
                self.R_nlr_list,self.R_final_nllkr_List,self.Spherical_Jn_Lqr_List,
                self.Integral,self.Atomic_Response_W1,self.Atomic_Response_K]

        for i in range(start,stop):
            if type(save_list[i])==type(None):
                print(name_list[i], 'is not calculated')
            else:
                Save( self.file_name , name_list[i] , save_list[i] )





    def _create_entire_quantum_number_list(self,Lmax=7):

        QNL_nllL = np.empty((0,4))

        for n in range(1,len(self.C)+1):    
            for l in range(len(self.C[n-1])):        
                for lPrime in range(Lmax+1):                            
                    for L in range(abs(l-lPrime), l+lPrime+1):                
                        vec = np.array([n,l,lPrime,L])
                        QNL_nllL = np.vstack([QNL_nllL,vec])
                                        
        self.QNL_nllL = QNL_nllL.astype(int)









    def set_r_grid(self,rmin,rmax,length):

        x_grid,weight = ssp.p_roots(length)
        r_grid = 0.5*(rmax-rmin)*x_grid+0.5*(rmax+rmin)

        self._N1 = int(length)
        self.r_grid = r_grid
        self.r_weight = weight
        self.rmin = r_grid[0]
        self.rmax = r_grid[-1]
        
        #print('rmin','%.3e'%self.r_grid[0],'rmax','%.3e'%self.r_grid[-1],'1/eV')
        






    def set_kPrime_grid(self,kPrime_grid):


        self._N2 = len(kPrime_grid)
        self.kPrime_grid = kPrime_grid
        
        #print('kmin','%.3e'%self.kPrime_grid[0],'kmax','%.3e'%self.kPrime_grid[-1],'eV')
    








    def set_q_grid(self,q_grid):

        self._N3 = len(q_grid)
        self.q_grid = q_grid

        #print('qmin','%.3e'%self.q_grid[0],'qmax','%.3e'%self.q_grid[-1],'eV')








    def set_demanded_n_l(self,n,l):
        
        if (n==None)*(l==None):
            QNL_nllL_select = self.QNL_nllL

        else:
            f1 = np.tile(self.QNL_nllL[:,:],(np.size(n),1,1))
            f2 = np.tile(n,(len(self.QNL_nllL),1)).T
            f3 = np.tile(l,(len(self.QNL_nllL),1)).T

            find = np.where(  (f1[:,:,0]==f2)*(f1[:,:,1]==f3) )
            QNL_nllL_select = f1[find[0],find[1],:]


        QNL_nl, QNL_nl_inverse = np.unique( QNL_nllL_select[:,[0,1]],axis=0, return_inverse=True)
        QNL_nll, QNL_nll_inverse = np.unique( QNL_nllL_select[:,[0,1,2]],axis=0, return_inverse=True)
        QNL_L, QNL_L_inverse = np.unique( QNL_nllL_select[:,[3]],axis=0, return_inverse=True)

        self.QNL_nllL_select = QNL_nllL_select

        self._QNL_nl = QNL_nl
        self._QNL_nl_inverse = QNL_nl_inverse
        self._QNL_nll = QNL_nll
        self._QNL_nll_inverse = QNL_nll_inverse
        self._QNL_L = QNL_L
        self._QNL_L_inverse = QNL_L_inverse 

        print((n,l))




    def _R_fun(self,n,l,r):

        r = np.tile(r,(1,1))

        C_vec = np.array(self.C[n-1][l])[:,None]
        Z_vec = np.array(self.Z[l])[:,None]
        n_vec = np.array(self.n_list[l])[:,None]
        factorize_vec =  ssp.factorial(2*n_vec)
        factor_vec = np.power(2*Z_vec, n_vec + 0.5) / np.sqrt(factorize_vec)
        
        rw =  np.sum( np.power(a0,-3./2.) * C_vec * factor_vec * np.power(r/a0, n_vec-1.)\
            * np.exp(-Z_vec * r/a0) ,axis=0)
        return rw


    def calculate_initial_rwf(self):

        size = len(self._QNL_nl)
        R_nlr_list = np.ndarray(( size , self._N1 ))
        for it in range(size): 
            n,l = self._QNL_nl[it]
            R_nlr_list[it] = self._R_fun(n,l,self.r_grid)

        self.R_nlr_list = R_nlr_list
        

    
    def _R_final_fun(self,n,l,lPrime,kPrime,r):

        Z_eff = self.Z_eff
        # notice it is a complex function
        hyp1f1 = np.vectorize(mp.hyp1f1)

        result = 4 * np.pi * (2*kPrime*r)**lPrime\
                * np.exp(np.pi * Z_eff[n-1][l] /2/kPrime/a0 + np.real(ssp.loggamma(lPrime+1-1j*Z_eff[n-1][l]/kPrime/a0)) )/ np.math.factorial(2*lPrime+1)\
                * np.exp(-1j*kPrime*r) \
                * hyp1f1(lPrime+1+1j*Z_eff[n-1][l]/kPrime/a0,(2*lPrime+2),2j*kPrime*r )  

        return result



    def calculate_final_rwf(self):
               
        
        kPrime_mesh, r_mesh = np.meshgrid(self.kPrime_grid, self.r_grid, indexing='ij')

        size = len(self._QNL_nll)
        R_final_nllkr_List = np.ndarray((size,self._N2,self._N1))  
        for it in tqdm.tqdm( range(size) ):           
            n,l,lPrime = self._QNL_nll[it]            
            R_final_nllkr_List[it]  = rv2float(rv2real(self._R_final_fun(n,l,lPrime,kPrime_mesh,r_mesh)))
         
        
        self.R_final_nllkr_List = R_final_nllkr_List
                  




                    

    def calculate_spherical_Jn_Lqr_List(self):


        q_mesh,r_mesh = np.meshgrid(self.q_grid,self.r_grid,indexing='ij')

        size = len(self._QNL_L)
        Spherical_Jn_Lqr_List = np.ndarray((size,self._N3,self._N1))
        for it in range(size):       
            L = self._QNL_L[it]              
            Spherical_Jn_Lqr_List[it] = ssp.spherical_jn(L ,q_mesh*r_mesh)   
  
        self.Spherical_Jn_Lqr_List = Spherical_Jn_Lqr_List 
         



       



    def get_I1(self):
 
        size = len(self.QNL_nllL_select)
        Integral = np.ndarray((size,self._N2,self._N3))

        for it in tqdm.tqdm( range(size) ) :
            it_nl = self._QNL_nl_inverse[it]
            it_nll = self._QNL_nll_inverse[it]
            it_L = self._QNL_L_inverse[it]

            Integrand = (self.r_grid**2)[None,None,None,:] * self.R_nlr_list[it_nl,None,None,:]\
                    * self.R_final_nllkr_List[it_nll,:,None,:] * self.Spherical_Jn_Lqr_List[it_L,None,:,:]
                     
            Integral[it] = 0.5*(self.rmax-self.rmin)* np.sum(self.r_weight[None,None,None,:]*Integrand,axis=3)  

        self.Integral = Integral
        
    




    def get_W1_atomic_response(self):

        # make a general list with m and m'
        QNL_nllLmm = np.empty((0,6))
    
        for n,l,lPrime,L in self.QNL_nllL_select:
            for m in range(-l,l+1):
                for mPrime in range(-lPrime,lPrime+1):
                    add = [n,l,lPrime,L,m,mPrime]
                    QNL_nllLmm = np.vstack([QNL_nllLmm,add])

        QNL_nllLmm = QNL_nllLmm.astype(int)
        QNL_nllLmm_Aux = np.unique(QNL_nllLmm[:,:4],axis=0,return_inverse=True)[1]


        # calculating integral's coefficients 
        coeff_f12 = np.empty(0)
        for n,l,lPrime,L,m,mPrime in QNL_nllLmm:            
            n,l,lPrime,L,m,mPrime = int(n),int(l),int(lPrime),int(L),int(m),int(mPrime)            
            add = np.sqrt(4*np.pi) * 1j**L * (-1.)**mPrime * np.sqrt(2*L+1) * \
                        float( gaunt(l,lPrime,L,m,-mPrime,0) )            
            coeff_f12 = np.append(coeff_f12,add)
        raw_f12 = coeff_f12[:,None,None] * self.Integral[QNL_nllLmm_Aux,:,:]
        del coeff_f12 



        # make a shorter list in order to sum up all the L
        QNL_nllmm, QNL_nllmm_Aux = np.unique( QNL_nllLmm[:,[0,1,2,4,5]] , axis=0 , return_inverse=True )
        QNL_nl , QNL_nl_Aux = np.unique( QNL_nllmm[:,[0,1]] , axis=0 , return_inverse=True )
        f12_gen = lambda index : np.sum(raw_f12[np.where(QNL_nllmm_Aux==index)[0],:,:],axis=0)
        f12 = np.zeros((len(QNL_nllmm) ,  self._N2   , self._N3  ), dtype=complex)

        for it in range(len(QNL_nllmm)):
            f12[it] = f12_gen(it)
        del raw_f12


        # calculating w1 without summation
        raw_w1 = np.abs(f12)**2
        del f12


        # calculating W1
        size = len(QNL_nl)
        Atomic_Response_W1 = np.zeros([size, self._N2, self._N3])
        for it in range(size):
            Atomic_Response_W1[it] = np.sum( raw_w1[ np.where(QNL_nl_Aux==it)[0],:,:] ,axis=0)\
             * 4 * (self.kPrime_grid**3)[None,:,None] / (2*np.pi)**3

        self.Atomic_Response_W1 = Atomic_Response_W1
        







    def get_K_atomic_response(self):
        
        W1 = self.Atomic_Response_W1
        Ee_grid = np.sqrt( self.kPrime_grid**2 + mElectron**2 ) - mElectron
        
        counts = len(self._QNL_nl)
        K = np.zeros(( counts ,self._N2,self._N3))
 
        for it in range(counts):
            n,l = self._QNL_nl[it]
            K[it] = W1[it]*(aEM*mElectron)**2/4./mElectron/Ee_grid[:,None]

        self.Atomic_Response_K = K
        






    def run_all_calculations(self):
        self.save_which(0,4)

        self.calculate_initial_rwf()
        self.calculate_final_rwf()
        self.calculate_spherical_Jn_Lqr_List()
        self.get_I1()
        self.save_which(4,8)


        self.get_W1_atomic_response()
        self.get_K_atomic_response()
        self.save_which(8,10)




