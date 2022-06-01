
AF 15/02/2022

##############################################
Q-WEIGHTING TOOLS CONSISTS OF THE FOLLOWING :
##############################################

Scripts :
---------
xtal_analysis_functions.py -> all functions required for scripts compiled into one file and installed as a package


Jupyter notebooks :
-------------------
make_qweighted.ipynb -> jupyter notebook to generate a qweighted map for a specific time point


Files :
-------
inputs/FC_dark.mtz -> Calculated structure factors
inputs/*-400nm_fobs_unique1.mtz (for every time point available) -> observed structure factor amplitudes with respective errors
outputs/ -> folder where output maps are stored

#######################
INSTRUCTIONS FOR USE :
#######################

Simply run each cell on the jupyter notebook. Change parameters at the top of the notebook as needed and any of the parameters specified at the top of xtal_analysis_functions.py. 

This applies q-weighting the same way as Marius Schmidt's FORTRAN scripts, though the scaling is applied in python. For better maps, use the files created here for noise filtering via PCA (in ../apply_PCA_denoise folder)

Currently changes in xtal_analysis_functions.py are applied because installed as "python setup.py develop" rather than "python setup.py install".
Might look into making a parameter file once development version finished. 