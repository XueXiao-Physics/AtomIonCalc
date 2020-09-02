from input_params import *
from methods import *
import multiprocess 
import sys
import os

#=======================================================================
#=======================================================================

'''

Interface

'''
# nomenclature
print(' > Welcome.')

arg = input('\n >>>>> Please type in the element you choose... \n > element : ')

elements = [Xe,O,Al,Ca,Fe,Ma,Ni,Si,S,Na]
element_dict = {e.name:e for e in elements}


if arg not in element_dict.keys():
    print(' > The element you ask for is not available.')
    sys.exit()
else:
    name,C,Z,n_list,E_B,Z_eff = element_dict[arg].call()
    
input("\n >>>>> Press Enter to create directory... <<<<<")  


#==========================================================================
# directory
try:
    os.mkdir('OUTPUT/'+name)
    print(' > \'OUTPUT/'+name+'\'','created')
    input("\n >>>>> Press Enter to continue... <<<<<")  
except FileExistsError:
    print(' > The','\'OUTPUT/'+name+'\'','direction already exists, \n > please delete or rename it first for data safety reasons,','\n > this program will not do it for you.')
    arg = input("\n >>>>> Are you sure you want to continue?(y/n)... <<<<<")
    if arg =='n':
        sys.exit()
    elif arg =='y':
        input("\n >>>>> Press Enter to use the already-exist directory... <<<<<") 
    else:
        print(' > Argument not recognized.')
        sys.exit()
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
# customizezd numeral accuracy
N1,N2,N3 = 2**12,2**7,2**7
rmin,rmax = 0,30*a0
kPrime_grid = np.logspace(-1,2,N2)*keV
q_grid = np.logspace(0,3,N3)*keV


print(' > Gridding parameters : ','\n > ',[N1,N2,N3],'\n')
print(' > rmin , rmax (1/eV) (Legendre Polynomials roots)','\n > ',[rmin,rmax],'\n')
print(' > kPrime_min , kPrime_max (1/eV)','\n > ', [kPrime_grid.min() , kPrime_grid.max()],'\n')
print(' > q_min , q_max (1/eV)','\n > ',[q_grid.min(),q_grid.max()],'\n')

input("\n >>>>> If alright, press Enter to start calculation ... <<<<<")
# interface end
#0.5*(rmax-rmin)*sum(weight*integrand)
#=======================================================================
#=======================================================================


pipelines = []

for it in range(len(combination)):

    c = combination[it]
    pipe = pipeline(C,Z,n_list,E_B,Z_eff)

    pipe.set_demanded_n_l( [c[0]] , [c[1]] )
    pipe.set_r_grid(rmin,rmax,N1) 
    pipe.set_kPrime_grid(kPrime_grid)
    pipe.set_q_grid(q_grid) 
    pipe.set_file_name( names[it] )


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


