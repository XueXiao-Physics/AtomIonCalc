## Manual
1. Use command below to calculate atomic response for each orbital
    python run.py
    
2. then type in the name of the element, eg.
    Xe
	
3. Then follow the guidance. Using "Xe" as an example, the results will be in the following directionary
    "OUTPUT/Xe/"

4. You can run the following code to calculate the total atomic response
    python run_Ktot.py
    
5. Type in the same element name and follow the guidance
    Xe

6. The result would be (using "Xe" as an example)
    "OUTPUT2/Xe_Ktot.hdf5" 
    
7. In "Xe_Ktot.hdf5" one will see these data
    "Ktot" : K(ER,q) (the first index represents ER)
    "q" : the grid of q (eV)
    "ER" : the grid of ER (recoil energy eV)
    
###################################################################  
> Notice: Only support python 3.
> Notice: Install h5py first.
    pip install h5py
> Notice: We suggest you install a higher mpmath version, for example
    pip install 'mpmath>=1.2.1'
> Contact:
    xxueitp@gmail.com
################################################################### 
## Citation:
If you use our code in your publications, you shall cite following papers
> Yifan Chen, Bartosz Fornal, Pearl Sandick, Jing Shu, Xiao Xue, Yue Zhao, and Junchao Zong. (2021). **Earth Shielding and Daily Modulation from Electrophilic Boosted Dark Matter**, [[arXiv:2110.09685]](https://arxiv.org/abs/2110.09685).
> Catena, R., Emken, T. , Spaldin, N., and Tarantino, W., **Atomic responses to general dark matter-electron interactions**, [[arXiv:1912.08204]](https://arxiv.org/abs/1912.08204).

and share the following link in your footnote.
3. https://github.com/XueXiao-Physics/AtomIonCalc 


