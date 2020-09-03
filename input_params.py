import numpy as np


eV = 1.0
keV = 1e3
GeV = 1e9
au = 27.211386245988*eV
gram = 5.60958884493318e23*GeV
cm = 5.067730214314311e13/GeV
meter = 100*cm
a0  =  5.29177e-11*meter #Bohr radius
mElectron =  511*keV
aEM  =  1.0/137.035999139

#=======================================================================
#=======================================================================
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

'''

Xe (Xenon)

'''  
C = [[[-0.965401, -0.040350, 0.001890, -0.003868, -0.000263,
         0.000547, -0.000791, 0.000014, -0.000013, -0.000286,
         0.000005, -0.000003, 0.000001] ] ,
     [[0.313912, 0.236118, -0.985333, 0.000229, -0.346825, 
        0.345786, -0.120941, -0.005057, 0.001528, -0.151508, 
        -0.000281, 0.000134, -0.000040] ,
      [0.051242, 0.781070, 0.114910, -0.000731, 0.000458,
       0.083993, -0.000265, 0.000034, 0.009061, -0.000014,
       0.000006, -0.000002] ] ,
     [[-0.140382, -0.125401, 0.528161, -0.000435, 0.494492, 
       -1.855445, 0.128637, -0.017980, 0.000792, 0.333907,
       -0.000228, 0.000191, -0.000037],
      [0.000264, 0.622357, -0.009861, -0.952677, -0.337900,
       -0.026340, -0.000384, -0.001665, 0.087491, 0.000240,
       -0.000083, 0.000026],
      [0.220185, 0.603140, 0.194682, -0.014369, 0.049865,
       -0.000300, 0.000418, -0.000133]] ,
     [[0.064020, 0.059550, -0.251138, 0.000152, -0.252274,
       1.063559, -0.071737, -0.563072, -0.697466, -0.058009,
       -0.018353, 0.00292, -0.000834],
      [0.013769, -0.426955, 0.045088, 0.748434, 0.132850,
       0.059406, -0.679569, -0.503653, -0.149635, -0.014193,
       0.000528, -0.000221],
      [-0.013758, -0.804573, 0.260624, 0.007490, 0.244109, 
       0.597018, 0.395554, 0.039786]] ,
     [[-0.022510, -0.021077, 0.088978, -0.000081, 0.095199,
       -0.398492, 0.025623, 0.274471, 0.291110, 0.011171,
       -0.463123, -0.545266, -0.167779],
      [-0.005879, 0.149040, -0.018716, -0.266839, -0.031096, 
       -0.024100, 0.267374, 0.161460, 0.059721, -0.428353,
       -0.542284, -0.201667]] ]

# Z_lj
Z = [ [54.9179, 47.25, 26.0942, 68.1771, 16.8296,
       12.0759, 31.903, 8.0145, 5.8396, 14.7123,
       3.8555, 2.6343, 1.8124] ,
     [58.7712, 22.6065, 48.9702, 13.4997, 9.8328,
      40.2591, 7.1841, 5.1284, 21.533, 3.4469,
      2.2384, 1.4588] , 
     [19.9787, 12.2129, 8.6994, 27.7398, 15.941, 
      6.058, 4.099, 2.5857] ]


n_list = [ [1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5] ,
     [2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5] ,
     [3, 3, 3, 4, 4, 4, 4, 4] ]

# Binding Energy
E_B = [ [-1224.397767 * au] , 
       [-189.230111 * au,-177.782438 * au] , 
       [-40.175652 * au,-35.221651 * au,-26.118859 * au] , 
       [-7.856291 * au,-6.008328 * au,-2.777871 * au] ,
       [-0.944407 * au,-0.457283 * au] ]



# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {}

Xe = element('Xe',C,Z,n_list,E_B,Z_eff,semi_full)


'''

O

'''  

# C_nlj
C = [[[0.360063,0.466625,-0.000918,0.208441,0.002018,0.000216,0.000133] ] ,
     [[-0.064363,0.433186,-0.000275,-0.072497,-0.369900,-0.512627,-0.227421], [0.005626, 0.126618, 0.328966, 0.395422,0.231788 ]]]

# Z_lj
Z = [ [11.2970, 6.5966,20.5019,9.5546,3.2482,2.1608,1.6411] , [9.6471 ,  4.3323, 2.7502, 1.7525, 1.2473] ]


n_list = [ [1, 1, 3 , 2, 2, 2, 2] , [2, 2, 2, 2, 2]]

# Binding Energy
E_B = [ [-20.668657 * au] ,  [-1.244315 * au, -0.631906 * au] ]

# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {(2,1):(4,6)}

O = element('O',C,Z,n_list,E_B,Z_eff,semi_full)


'''

Al

'''
# C_nlj
C = [[[0.373865,0.456146,0.20250,0.001901,0.000823,-0.000267,-0.000560,0.000083,-0.000044,0.000013] ] ,
[[0.061165,-0.460373,0.055062,0.297052,0.750997,0.064079,0.000270,-0.001972,0.000614,-0.000064] ,[0.015480,0.204774,0.474317,0.339646,0.024290,0.003529,-0.000204,0.000199]],
 [[0.020024 , -0.119051 , 0.017451 , 0.079185 , 0.130917 , 0.13913 , 0.000038 , -0.303750 , -0.547941 , -0.285949], [ -0.001690 , -0.048903 , -0.058101 , -0.090680 , -0.001445 , 0.234760 , 0.496072 , 0.359277  ]]  ]

# Z_lj
Z = [ [18.1792,10.8835,15.7593,5.7600,4.0085,2.8676,33.5797,2.1106,1.3998,1.0003] ,
     [14.4976, 6.6568, 4.2183, 3.0026, 11.0822, 1.6784, 1.0788, 0.7494] ]

n_list = [ [1, 1, 2, 2, 2, 2, 3, 3, 3, 3] ,[2, 2, 2, 2, 3, 3, 3, 3] ]

# Binding Energy
E_B = [ [-58.501026 * au] ,[-4.910672 * au,-3.218303 * au] , [-0.393420 * au,-0.209951 * au] ]

# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {(3,1):(1,6)}

Al = element('Al',C,Z,n_list,E_B,Z_eff,semi_full)


'''

Ca

'''

# C_nlj
C = [[[1.176575, -0.119414, 0.000730, -0.086509, 0.000185, -0.000097, 0.000145, -0.000317, 0.000021, -0.000015, 0.000005] ] ,
     [[-0.341543, -0.019904, 1.150570, -0.033142, -0.000762, 0.001467, -0.068959, 0.016219, -0.000184, 0.000124, -0.000044] ,
      [0.001400, 0.348898, 0.406540, 0.003420, 0.300142, -0.000861, 0.000083] ] ,
     [[-0.117803, -0.010168, 0.468206, -0.016715, -0.640251, -0.441778, 0.043103, -0.064309, -0.001947, 0.000314, -0.000058], [-0.000342, -0.076856, -0.433592, 0.877080, -0.015449, 0.312083, 0.020164]],
     [[0.028125,0.002534,-0.113880,0.004109,0.179379,0.115506,-0.014022,0.020813,-0.331747,-0.522462,-0.257798]]]

# Z_lj
Z = [ [17.6670, 24.5295, 7.2525, 21.3493, 3.6185, 2.6126, 8.7808, 6.0113, 1.5770, 1.0566, 0.7584] ,
     [28.8909, 11.9319, 6.0967, 2.5605, 9.4272, 2.1782, 1.3295] ]

n_list = [ [1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4] ,
     [2, 2, 2, 2, 3, 3, 3] ]

# Binding Energy
E_B = [ [-149.363724 * au] , 
       [-16.822743 * au,-13.629270 * au] , 
       [-2.245375 * au,-1.340706 * au] ,
       [-0.195529 * au] ]

# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full={}

Ca = element('Ca',C,Z,n_list,E_B,Z_eff,semi_full)


'''

Fe

'''
# C_nlj
C = [[[0.942524, 0.063210, 0.000162, 0.008635, -0.000114, -0.000058, 0.000005, 0.000297, -0.000005, 0.000003, -0.000001] ] ,
     [[-0.285853, -0.189028, 1.008545, -0.000649, 0.093973, -0.007693, 0.003112, 0.087617, -0.000360, 0.000209, -0.000070] ,
      [0.000583, 0.323715, 0.434152, 0.006165, 0.288160, -0.000634, 0.000215] ] ,
     [[0.105891, 0.081539, -0.442968, 0.000883, -0.101287, 0.635102, 0.577033, -0.103885, 0.007684, -0.001785, 0.000454], [0.000148, 0.068651, 0.516640, -0.821240, 0.017006, -0.407245, -0.041825], [0.017645,0.216184,0.414852,0.383883,0.129854]],[[0.022650,0.017743,-0.096526,0.000194,-0.028327,0.171358,0.129599,-0.024523,-0.275750,-0.539740,-0.312861]]]

# Z_lj
Z = [ [26.7103, 22.7394, 11.1579, 32.8592, 8.2265, 5.7011, 3.8393, 13.2036, 2.2676, 1.4378, 0.9462] ,
     [42.2924, 16.3351, 8.7201, 4.0978, 13.2514, 3.3690, 2.1919], 
     [12.7458,6.6281,4.0177,2.4057,1.4774]]

n_list = [ [1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4] ,
     [2, 2, 2, 2, 3, 3, 3] ,
    [3, 3, 3, 3, 3]]

# Binding Energy
E_B = [ [-261.373415 * au] , 
       [-31.935515 * au,-27.413707 * au] , 
       [-4.169434 * au,-2.742193 * au, -0.646883 * au],
       [-0.258178 * au]]

# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {(3,2):(6,10)}

Fe = element('Fe',C,Z,n_list,E_B,Z_eff,semi_full)

'''

Mg

'''

# C_nlj
C = [[[0.352464,0.481225,0.198592,0.002259,0.000556,-0.000136,-0.000669,0.000056,-0.000033,0.000011] ],  [[0.059265,-0.447481,0.055907,0.355163,0.696633,0.058440,0.000283,-0.001173,0.000277,-0.000059] ,[0.004178,0.175692,0.420054,0.456246,0.012155]],
[[0.016053 , -0.096426 , 0.014785 , 0.077390 , 0.110979 , 0.082870 , 0.000010 , -0.232777 , -0.494745 , -0.378869]] ]

# Z_lj
Z = [ [17.0241,10.0727,14.6751,5.1514,3.4870,2.5249,29.9018,1.7568,1.1659,0.8244] ,
     [14.9021,6.8076,4.1426,2.7152,1.4623] ]

n_list = [ [1, 1, 2, 2, 2, 2, 3, 3, 3, 3] ,
     [2, 2, 2, 2, 2] ]

# Binding Energy
E_B = [ [-49.031735 * au] ,
       [-3.767721 * au,-2.282225 * au],
      [-0.253052 * au] ]

# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {}

Mg = element('Mg',C,Z,n_list,E_B,Z_eff,semi_full)

'''

Ni

'''

C = [[[0.945523, 0.059970, 0.000155, 0.008101, -0.000099, -0.000046, 0.000008, 0.000226, -0.000005, 0.000003, -0.000001]],
[[0.289565, 0.187069, -1.073281, 0.000139, -0.063716, 0.009081, -0.003291, -0.043470, 0.000399, -0.000214, 0.000071],
[-0.000470, -0.319845, -0.438165, -0.006828, -0.285978, 0.000547, -0.000255]],
[[-0.108484, -0.083030, 0.486257, -0.000950, 0.00334, -0.685620, -0.475026, 0.106095, -0.004158, 0.000568, -0.000119],
[-0.000119, -0.067683, -0.527884, 0.810679, -0.020283, 0.421372, 0.049664],
[0.017517, 0.225100, 0.419212, 0.376659, 0.130041]],
[[0.022468, 0.017547, -0.102978, 0.000219, -0.004844, 0.182371, 0.094382, -0.024948, -0.283769, -0.538297, -0.308953]]]

Z = [ [28.7266, 24.5237, 11.7009, 35.2606, 8.3190, 5.8931, 4.0420, 14.1323, 2.4332, 1.5225, 0.9912] ,
     [46.7514, 17.7675, 9.5726, 4.5914, 14.5116, 3.7535, 2.4658], 
     [13.9008, 7.3628, 4.4602, 2.6365, 1.5875]]

n_list = [ [1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4] ,
     [2, 2, 2, 2, 3, 3, 3] ,
    [3, 3, 3, 3, 3]]

E_B = [ [-305.619027 * au] , 
       [-37.917819 * au,-32.941726 * au] , 
       [-4.887827 * au, -3.277669 * au, -0.706924 * au],[-0.276245 * au]]

Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {(3,2):(8,10)}

Ni = element('Ni',C,Z,n_list,E_B,Z_eff,semi_full)


'''

Si

'''

# C_nlj
C = [[[0.377006, 0.454461, 0.200676, 0.001490, 0.001201, -0.000454, -0.000507, 0.000103, -0.000053, 0.000013] ] ,
[[0.064222, -0.472631, 0.055333, 0.233799, 0.781919, 0.096627, 0.000257, -0.001832, 0.000879, -0.000033] ,[0.015661, 0.196557, 0.510448, 0.303956, 0.025586, 0.003153, 0.000167, 0.000156]],[[0.023528,-0.136207,0.019663,0.074362,0.122580,0.206180,0.000048,-0.319063,-0.562578,-0.280471], [-0.001966,-0.057175,-0.068127,-0.114298,-0.001976,0.263703,0.522698,0.314467]] ]

# Z_lj
Z = [ [19.5017, 11.7539, 16.9664, 6.3693, 4.5748, 3.3712, 36.5764, 2.4996, 1.6627, 1.1812] ,
     [15.7304, 7.2926, 4.6514, 3.3983, 12.0786, 2.0349, 1.3221, 0.9143] ]

n_list = [ [1, 1, 2, 2, 2, 2, 3, 3, 3, 3] ,
     [2, 2, 2, 2, 3, 3, 3, 3] ]

# Binding Energy
E_B = [ [-68.812456 * au] ,
       [-6.156538 * au,-4.256054 * au] ,
       [-0.539842 * au,-0.297114 * au]]

# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {(3,1):(2,6)} 

Si = element('Si',C,Z,n_list,E_B,Z_eff,semi_full)


'''

S

'''

# C_nlj
C = [[[0.367468,0.471254,0.192030,0.000539,0.002074,-0.000864,-0.000442,0.000163,-0.000073,0.000014] ] ,
[[0.070509,-0.492518,0.056472,0.114779,0.846899,0.150553,0.000233,-0.002483,0.001580,0.000045] ,
[-0.000466,0.141231,0.501894,0.403324,-0.006509,0.004375,0.000225,0.000315] ] ,
[[0.028833,-0.160596,0.022476,0.056607,0.107272,0.326955,0.000067,-0.337632,-0.602710,-0.272978], 
[0.001263,-0.047039,-0.090305,-0.159888,0.005017,0.341230,0.519259,0.262504]] ]

# Z_lj
Z = [ [22.2949, 13.5666, 19.4969, 7.5145, 5.7222, 4.2264, 42.1787, 3.3088, 2.1707, 1.5140] ,
     [22.6414, 10.4197, 6.1160, 4.4156, 17.3448, 2.6496, 1.6975, 1.1477] ]

n_list = [ [1, 1, 2, 2, 2, 2, 3, 3, 3, 3] ,
     [2, 2, 2, 2, 3, 3, 3, 3] ]

# Binding Energy
E_B = [ [-92.004449 * au] , 
       [-9.004288 * au,-6.682508 * au] , 
       [-0.879527 * au,-0.437368 * au] ]

# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

semi_full = {(3,1):(4,6)}
S = element('S',C,Z,n_list,E_B,Z_eff,semi_full)


'''

Na



# C_nlj
C = [[[0.387167,0.434278,0.213027,0.002205,0.000627,-0.000044,-0.000649,0.000026,-0.000023,0.000008] ] ,
[[0.053722,-0.430794,0.053654,0.347971,0.608890,0.157462,0.000280,-0.000492,0.000457,0.000016],[0.004308,0.157824,0.388545,0.489339,0.039759]], 
[[0.011568,-0.072430,0.011164,0.057679,0.089837,0.042114,-0.000001,-0.182627,-0.471631,-0.408817]]]

# Z_lj
Z = [ [15.3319,9.0902,13.2013,4.7444,3.1516,2.4047,28.4273,1.3179,0.8911,0.6679] ,
     [13.6175,6.2193,3.8380,2.3633,1.5319] ]

n_list = [ [1, 1, 2, 2, 2, 2, 3, 3, 3, 3] ,[2, 2, 2, 2, 2] ]

# Binding Energy
E_B = [ [-40.478500 * au] , 
       [-2.797026 * au,-1.518140 * au] ,
      [-0.182102 * au]  ]

# Effective charge
Z_eff=[]
for n in range(len(E_B)):
    Z_eff.append([])
    for l in range(len(E_B[n])):
        Z_eff[n].append((n+1) * np.sqrt(-2*E_B[n][l] / au))

Na = element('Na',C,Z,n_list,E_B,Z_eff)
'''
