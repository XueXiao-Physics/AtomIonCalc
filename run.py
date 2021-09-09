from methods import *
import multiprocess 
import sys
sys.path.append('input_params')
import os
import glob

element_param_files = glob.glob('input_params/Param_*.py')
element_names = [elem_file.split('/')[-1][6:-3] for elem_file in element_param_files]
#=======================================================================
#=======================================================================

'''

Interface

'''
# nomenclature
print(' > Welcome. Run the code to calculate the atomic response of each orbit.')
print(' > All elements available:')
print(' > ',end='')
for elem_name in element_names:
    print(elem_name,end='  ')
print('')



arg = input('\n >>>>> Please type in the element you choose: \n > element : ')

#elements = [Xe,O,Al,Ca,Fe,Mg,Ni,Si,S,Na]



if arg not in element_names:
    print(' > The element you ask for is not available.')
    sys.exit()
else:
    elem = __import__('Param_'+arg)
    name,C,Z,n_list,E_B,Z_eff,semi_full = elem.elem.call()
    
input("\n >>>>> Press Enter to create directory... <<<<<")  


#==========================================================================
# directory
try:
    try:
        os.mkdir('OUTPUT/')
    except  FileExistsError:
        pass
        
    os.mkdir('OUTPUT/'+name)
    print(' > \'OUTPUT/'+name+'\'','created')
    
except FileExistsError:
    print(' > Careful, \'OUTPUT/'+name+'\'','already exists.') 
    input("\n >>>>> Press Enter to continue... <<<<<")  

'''
    print(' > The','\'OUTPUT/'+name+'\'','directionary already exists, \n > please delete or rename it first for data safety reasons,','\n > this program will not do it for you.')
    arg = input("\n >>>>> Are you sure you want to continue?(y/n)... <<<<<")
    if arg =='n':
        sys.exit()
    elif arg =='y':
        input("\n >>>>> Press Enter to use the already-exist directory... <<<<<") 
    else:
        print(' > Argument not recognized.')
        sys.exit()
'''
print("\n\n\n\n ============================== INFO ============================== \n\n")



#==========================================================================  
#preparing

file_name = 'OUTPUT/'+name+'/'+name
orbit = ['s','p','d']
combination = []
names = []
for n in range(1,1+len(C)):
    for l in range(len(C[n-1])):
        combination.append((n,l))
        names.append(file_name+str(n)+orbit[l])
print(' > All (n,l) combinations :')
print(' > ',combination)

input("\n >>>>> If alright, press Enter to continue... <<<<<")


#=======================================================================
print("\n\n\n\n ============================== INFO ============================== \n\n")
# customizezd numerical accuracy
Nr,Nk,Nq = 2**12,2**7,2**8
rmin,rmax = 0,30*a0
kPrime_grid = np.logspace(-1,2,Nk)*keV # Define it freely.
q_grid = np.logspace(0,3,Nq)*keV # Define it freely.


print(' > Gridding parameters : ','\n > ',[Nr,Nk,Nq],'\n')
print(' > rmin , rmax (1/eV) (Legendre Polynomials roots)','\n > ',[rmin,rmax],'\n')
print(' > kPrime_min , kPrime_max (eV)','\n > ', [kPrime_grid.min() , kPrime_grid.max()],'\n')
print(' > q_min , q_max (eV)','\n > ',[q_grid.min(),q_grid.max()],'\n')

input("\n >>>>> If alright, press Enter to start calculation ... <<<<<")
# interface end
#0.5*(rmax-rmin)*sum(weight*integrand)
#=======================================================================
#=======================================================================


pipelines = []

for it in range(len(combination)):

    c = combination[it]

    try: 
        fac = semi_full[c][0]/semi_full[c][1] 
        print(c,'%.2f'%fac)
    except KeyError: 
        fac = 1.
            
    pipe = pipeline(C,Z,n_list,E_B,Z_eff,fac,names[it] )

    pipe.set_demanded_n_l( c[0] , c[1] )
    pipe.set_r_grid(rmin,rmax,Nr) 
    pipe.set_kPrime_grid(kPrime_grid)
    pipe.set_q_grid(q_grid) 
    
    pipelines.append(pipe)
    


jobs = []
for it in range(len(pipelines)):
    
    pipe = pipelines[it]
    job = multiprocess.Process(target = pipe.run_all_calculations, args=() )
    job.start()
    jobs.append(job)

for job in jobs:
    job.join()


input('\n >>>>> Done. Press Enter to exit... <<<<<')


