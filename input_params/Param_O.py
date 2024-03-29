from Element import element,au
import numpy as np

C = [[[0.360063,0.466625,-0.000918,0.208441,0.002018,0.000216,0.000133] ] ,
     [[-0.064363,0.433186,-0.000275,-0.072497,-0.369900,-0.512627,-0.227421], [0.005626, 0.126618, 0.328966, 0.395422,0.231788 ]]]

Z = [ [11.2970, 6.5966,20.5019,9.5546,3.2482,2.1608,1.6411] , [9.6471 ,  4.3323, 2.7502, 1.7525, 1.2473] ]


n_list = [ [1, 1, 3 , 2, 2, 2, 2] , [2, 2, 2, 2, 2]]

E_B = [ [-20.668657 * au] ,  [-1.244315 * au, -0.631906 * au] ]

Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {(2,1):(4,6)}

elem = element('O',C,Z,n_list,E_B,Z_eff,semi_full)
