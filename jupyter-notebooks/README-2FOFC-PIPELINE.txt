
AF 15/02/2022

######################################

OBTAIN A Nbgmax VALUE FOR EACH DATA
POINT STARTING FROM 2mFo-DFc MAPS 
AND SAVE EACH (LIGHT - Nbgmax*DARK) 
MAP TO EXAMINE

######################################
PIPELINE CONSISTS OF THE FOLLOWING :
######################################

Scripts :
---------
xtal_analysis_functions.py -> all functions required for scripts compiled into one file and installed as a package


Jupyter notebooks :
-------------------
2Fo-Fc_correlate_pipeline.ipynb -> jupyter notebook guiding identification of Nbg_max for a specific time point


Files :
-------
dark.pdb -> refined dark model 
input_maps/*_nochrom_2mFo-DFc_map.ccp4 (for every time point available) -> 2mFo-DFc map for every time point generated as specified in input_maps/NOTES.txt 


#######################
INSTRUCTIONS FOR USE :
#######################

Simply run each cell on the jupyter notebook. Change parameters at the top of the notebook as needed and any of the parameters specified at the top of xtal_analysis_functions.py. 

Currently changes in xtal_analysis_functions.py are applied because installed as "python setup.py develop" rather than "python setup.py install".
Might look into making a parameter file once development version finished. 

