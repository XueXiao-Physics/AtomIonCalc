from input_params.input_params import element,au
import numpy as np

C = [[[1.176575, -0.119414, 0.000730, -0.086509, 0.000185, -0.000097, 0.000145, -0.000317, 0.000021, -0.000015, 0.000005] ] ,
     [[-0.341543, -0.019904, 1.150570, -0.033142, -0.000762, 0.001467, -0.068959, 0.016219, -0.000184, 0.000124, -0.000044] ,
      [0.001400, 0.348898, 0.406540, 0.003420, 0.300142, -0.000861, 0.000083] ] ,
     [[-0.117803, -0.010168, 0.468206, -0.016715, -0.640251, -0.441778, 0.043103, -0.064309, -0.001947, 0.000314, -0.000058], [-0.000342, -0.076856, -0.433592, 0.877080, -0.015449, 0.312083, 0.020164]],
     [[0.028125,0.002534,-0.113880,0.004109,0.179379,0.115506,-0.014022,0.020813,-0.331747,-0.522462,-0.257798]]]

Z = [ [17.6670, 24.5295, 7.2525, 21.3493, 3.6185, 2.6126, 8.7808, 6.0113, 1.5770, 1.0566, 0.7584] ,
     [28.8909, 11.9319, 6.0967, 2.5605, 9.4272, 2.1782, 1.3295] ]

n_list = [ [1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4] ,
     [2, 2, 2, 2, 3, 3, 3] ]

E_B = [ [-149.363724 * au] , 
       [-16.822743 * au,-13.629270 * au] , 
       [-2.245375 * au,-1.340706 * au] ,
       [-0.195529 * au] ]

Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full={}

elem = element('Ca',C,Z,n_list,E_B,Z_eff,semi_full)
