## Manual
> Use command below to calculate atomic response for each orbital
> 
    python run.py
> type in the name of the element, eg. Xe, then follow the guidance. The results will be in "OUTPUT/Xe/".
> You can also run the following code to calculate the total atomic response
> 
    python run_Ktot.py
    
> type in the same element name (eg. Xe) and follow the guidance.
> The result would be in "OUTPUT2/Xe_Ktot.hdf5". 
> In "Xe_Ktot.hdf5" one will see these data
> 
    "Ktot" : K(ER,q) (the first index represents ER)
    "q" : the grid of q (eV)
    "ER" : the grid of ER (recoil energy eV)
    
> Note that the program only supports python 3. To properly run the program, you need to install h5py and multiprocess first.
> 
    pip install h5py
    pip install multiprocess

## Citation:
If you use our code in your publications, you shall cite following papers
> Yifan Chen, Bartosz Fornal, Pearl Sandick, Jing Shu, Xiao Xue, Yue Zhao, and Junchao Zong. (2021). **Earth Shielding and Daily Modulation from Electrophilic Boosted Dark Matter**, [[arXiv:2110.09685]](https://arxiv.org/abs/2110.09685).


> Catena, R., Emken, T. , Spaldin, N., and Tarantino, W., **Atomic responses to general dark matter-electron interactions**, [[arXiv:1912.08204]](https://arxiv.org/abs/1912.08204).

and share the following link in your footnote.
> https://github.com/XueXiao-Physics/AtomIonCalc 


