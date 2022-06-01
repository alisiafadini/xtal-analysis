
AF 15/02/2022

######################################

OBTAIN A Nbgmax VALUE FOR EACH DATA
POINT STARTING FROM OBSERVED STRUCTURE
FACTORS AND SAVE EACH 
(LIGHT - Nbgmax*DARK) MAP TO EXAMINE

######################################
PIPELINE CONSISTS OF THE FOLLOWING :
######################################

Scripts :
---------
xtal_analysis_functions.py -> all functions required for scripts compiled into one file and installed as a package


Jupyter notebooks :
-------------------
qweight_correlate_pipeline.ipynb -> jupyter notebook guiding identification of Nbg_max for a specific time point


Files :
-------
dark.pdb -> refined dark model 
input_maps/FC_dark.mtz -> Calculated structure factors
input_maps/*-400nm_fobs_unique1.mtz (for every time point available) -> observed structure factor amplitudes with respective errors


#######################
INSTRUCTIONS FOR USE :
#######################

Simply run each cell on the jupyter notebook. Change parameters at the top of the notebook as needed and any of the parameters specified at the top of xtal_analysis_functions.py. 

Currently changes in xtal_analysis_functions.py are applied because installed as "python setup.py develop" rather than "python setup.py install".
Might look into making a parameter file once development version finished. 

