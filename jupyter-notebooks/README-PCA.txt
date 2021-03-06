  
AF 15/02/2022

######################################

OBTAIN DIFFERENCE QWEIGHTED MAPS
RECONSTRUCTED WITHOUT NOISE THROUGH
PCA ANALYSIS 

######################################
FOLDER CONSISTS OF THE FOLLOWING :
######################################

Scripts :
---------
xtal_analysis_functions.py -> all functions required for scripts compiled into one file and installed as a package


Jupyter notebooks :
-------------------
make_chromophore_maps.ipynb -> writes maps that are masked around a specific region, given an input CCP4 map. Choice of spherical or rectangular masks, xyz limits specified by user. 

apply_PCA.ipynb -> loads in N maps, creates a data matrix, and performs PCA. Allows visualization of PCA results and selection of number of components. Maps can be reconstructed and written to output_maps folder. 


Files :
-------
dark.pdb -> refined dark model 

input_maps/*_diff.ccp4 -> qweighted maps (entire protein) generated from ../apply_qweighting methods. One for each time point collected

output_maps/*_masked.ccp4 -> qweighted maps where a spherical mask around chromophore has been applied. These can be used as inputs for the PCA analysis instead of the full protein maps

output_maps/*_rect_masked.ccp4 -> qweighted maps where a rectangular mask around chromophore has been applied. These can be used as inputs for the PCA analysis instead of the full protein maps

output_maps/pca_reconstructed_t*_{mask?}.map -> reconstructud CCP4 difference maps, one for each time point. Mask status (entire protein or chromophore) specified



#######################
INSTRUCTIONS FOR USE :
#######################

Simply run each cell on the jupyter notebook. Change parameters at the top of the notebook as needed and any of the parameters specified at the top of xtal_analysis_functions.py. 

Decide whether to use full protein maps or one of the masked options generated by make_chromophore_maps.ipynb

Currently changes in xtal_analysis_functions.py are applied because installed as "python setup.py develop" rather than "python setup.py install".
Might look into making a parameter file once development version finished. 

