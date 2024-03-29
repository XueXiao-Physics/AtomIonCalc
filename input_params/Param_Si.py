from Element import element,au
import numpy as np

C = [[[0.377006, 0.454461, 0.200676, 0.001490, 0.001201, -0.000454, -0.000507, 0.000103, -0.000053, 0.000013] ] ,
[[0.064222, -0.472631, 0.055333, 0.233799, 0.781919, 0.096627, 0.000257, -0.001832, 0.000879, -0.000033] ,[0.015661, 0.196557, 0.510448, 0.303956, 0.025586, 0.003153, 0.000167, 0.000156]],[[0.023528,-0.136207,0.019663,0.074362,0.122580,0.206180,0.000048,-0.319063,-0.562578,-0.280471], [-0.001966,-0.057175,-0.068127,-0.114298,-0.001976,0.263703,0.522698,0.314467]] ]


Z = [ [19.5017, 11.7539, 16.9664, 6.3693, 4.5748, 3.3712, 36.5764, 2.4996, 1.6627, 1.1812] ,
     [15.7304, 7.2926, 4.6514, 3.3983, 12.0786, 2.0349, 1.3221, 0.9143] ]

n_list = [ [1, 1, 2, 2, 2, 2, 3, 3, 3, 3] ,
     [2, 2, 2, 2, 3, 3, 3, 3] ]


E_B = [ [-68.812456 * au] ,
       [-6.156538 * au,-4.256054 * au] ,
       [-0.539842 * au,-0.297114 * au]]


Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {(3,1):(2,6)} 

elem = element('Si',C,Z,n_list,E_B,Z_eff,semi_full)

