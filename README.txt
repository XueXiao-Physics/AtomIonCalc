# Please use python 3 instead of 2.
# Use command below to calculate atomic response for each orbital
    python run.py
# then type in the name of the element (eg. "Xe", no ") and follow the guide.
# After finishing the calculation, run
    python run_Ktot.py
# then type in the same element name (eg. "Xe", no ") and follow the guide.
# The result would be (using "Xe" as an example)
    OUTPUT2/"Xe_Ktot.hdf5" 
# In "Xe_Ktot.hdf5" one will see these data
    "Ktot : K(q,ER)
    "q" : the grdding of q (eV)
    "ER" : the gridding of ER (eV)
###################################################################  
# Notice: please do install h5py first
    pip install h5py
    
# Contact:
    xxueitp@gmail.com

