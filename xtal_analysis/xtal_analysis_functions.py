#!/usr/bin/env python


import matplotlib.pyplot as plt
import numpy as np
import reciprocalspaceship as rs
import scipy.optimize as opt
import gemmi as gm
import os
from   xtal_analysis.params import *
from   tqdm   import tqdm
from skimage.restoration import denoise_tv_chambolle
import diptest
from scipy.stats import norm, kurtosis, skew, differential_entropy


"""
This file contains the functions required to run the following jupyter notebooks :

	(1) qweight_correlate_pipeline.ipynb
	(2) 2fo-fc_correlate_pipeline.ipynb
	(3) pca_maps.ipynb

Other requirements :
    
    edit params.py file for all dataset-specific/personal choices --> RIGHT NOW THIS WORKS IF LIBRARY IS SETUP AS 'DEVELOP'
	make_ccp4map.sh script (appropriately filled in) in the same path as specified by the screen_qweighted function


AF 15/02/2022

"""

    
def map2mtz(map, mtz_name, high_res):

    """
    Writes an mtz from a map file

    Input :

    1. map path – string
    2. output mtz path – string
    3. high resolution cutoff – float

    """

    m = gm.read_ccp4_map(map, setup=True)
    sf = gm.transform_map_to_f_phi(m.grid, half_l=False)
    data = sf.prepare_asu_data(dmin=high_res-0.03, with_sys_abs=True)
    mtz = gm.Mtz(with_base=True)
    mtz.spacegroup = sf.spacegroup
    mtz.set_cell_for_all(sf.unit_cell)
    mtz.add_dataset('unknown')
    mtz.add_column('FWT', 'F')
    mtz.add_column('PHWT', 'P')
    mtz.set_data(data)
    mtz.switch_to_asu_hkl()
    mtz.write_to_file(mtz_name)
    

def make_test_set(df, percent, Fs, out_name):

    """
    Writes an mtz where data from the original mtz has been divided in a "fit" set and a "test" set.
    Also saves test set indices as numpy object.

    Input :

    1. dataset mtz to split – in reciprocalspaceship format
    2. fraction of reflections to keep for test set – float e.g. 0.03
    3. labels for structure factors to split – string
    4. output file specification – string

    """

    choose_test = np.random.binomial(1, percent, df[Fs].shape[0]).astype(bool)
    test_set = df[Fs][choose_test] #e.g. 3%
    fit_set  = df[Fs][np.invert(choose_test)] #97%
    
    df["fit-set"]   = fit_set
    df["fit-set"]   = df["fit-set"].astype("SFAmplitude")
    df["test-set"]  = test_set
    df["test-set"]  = df["test-set"].astype("SFAmplitude")
    
    df.write_mtz("split-{}.mtz".format(out_name))
    np.save("test_flags-{}.npy".format(out_name), choose_test)
    
    return test_set, fit_set

def TV_filter(map, l, name):
    
    """
    Applies TV filtering to a map and saves it in directory.
    Computes negentropy and diptest statistic for denoised map.

    Input :

    1. map element (as loaded by GEMMI)
    2. weighting for filtering – float
    3. name specification for output map – string

    Returns :

    1. filtered map – 3D numpy array
    2. negentropy value – float
    3. diptest value – float
    
    """

    map_TV = denoise_tv_chambolle(np.array(map.grid), weight=l)
    entropy    = negentropy(map_TV.flatten())
    dip        = diptest.dipstat(map_TV.flatten())
    save_map('','{name}_map_TV_{l:.5f}'.format(name=name, l=l), map_TV)

    return map_TV, entropy, dip

def positive_Fs(df, phases, Fs, phases_new, Fs_new):

    """
    Converts between mtz format where deltaFs are saved as both positive and negative, to format where they are only positive.

    Input :

    1. mtz (ideally cut at resolution already) –reciprocalspaceship format
    2. label for phases in original mtz – string
    3. label for Fs in original mtz – string
    4. label for new phases in output mtz – string
    5. label for new Fs in output mtz – string

    Returns :

    1. mtz dataset with new columns added – reciprocalspaceship format
    
    """
    
    new_phis = df[phases].copy(deep=True)
    new_Fs   = df[Fs].copy(deep=True)
    
    negs = np.where(df[Fs]<0)
    
    for i in negs:
        new_phis.iloc[i]  = df[phases].iloc[i]+180
        new_Fs.iloc[i]    = np.abs(new_Fs.iloc[i])

    df_new = df.copy(deep=True)
    df_new[Fs_new]  = new_Fs
    df_new[Fs_new]  = df_new[Fs_new].astype("SFAmplitude")
    df_new[phases_new]  = new_phis
    df_new[phases_new]  = df_new["new_og-Phis"].astype("Phase")
    
    return df_new


def negentropy(x):
    
    # negetropy is the difference between the entropy of samples x
    # and a Gaussian with same variance
    # http://gregorygundersen.com/blog/2020/09/01/gaussian-entropy/
    
    std = np.std(x)
    neg_e = np.log(std*np.sqrt(2*np.pi*np.exp(1))) - differential_entropy(x)
    #neg_e = 0.5 * np.log(2.0 * np.pi * std ** 2) + 0.5 - differential_entropy(x)
    #assert neg_e >= 0.0
    
    return neg_e



def scale_iso(data1, data2, ds):

    """
    Isotropic resolution-dependent scaling of data2 to data1.
    (minimize [dataset1 - c*exp(-B*sintheta**2/lambda**2)*dataset2]

    Input :

    1. dataset1 in form of 1D numpy array
    2. dataset2 in form of 1D numpy array
    3. dHKLs for the datasets in form of 1D numpy array

    Returns :

    1. entire results from least squares fitting
    2. c (as float)
    3. B (as float)
    2. scaled dataset2 in the form of a 1D numpy array

    """
        
    def scale_func(p, x1, x2, qs):
        return x1 - (p[0]*np.exp(-p[1]*(qs**2)))*x2
    
    p0 = np.array([1.0, -20])
    qs = 1/(2*ds)
    matrix = opt.least_squares(scale_func, p0, args=(data1, data2, qs))
    
    return matrix, matrix.x[0], matrix.x[1], (matrix.x[0]*np.exp(-matrix.x[1]*(qs**2)))*data2

def map_from_Fs(path, Fs, phis, map_res):
    
    mtz  = gm.read_mtz_file('{}'.format(path))
    ccp4 = gm.Ccp4Map()
    ccp4.grid = mtz.transform_f_phi_to_map('{}'.format(Fs), '{}'.format(phis), sample_rate=map_res)
    ccp4.update_ccp4_header(2, True)
    
    return ccp4

def load_ccp4(ccp4):
    
    """
    Loads ccp4 map with GEMMI

    Input :
	
	1. electron density in CCP4 map format

    Returns :

	1. entire map as a GEMMI object
	2. GEMMI grid of map
	3. map data points as a 3D numpy array  
	4. flattened map array (to be used e.g. for constructing PCA data matrix)

    """
 
    dataset = gm.read_ccp4_map(ccp4, setup=True)
    #dataset.setup()
    grid    = dataset.grid
    array   = np.array(dataset.grid, copy=True)
    array_f = np.array(dataset.grid, copy=True).flatten()
    
    return dataset, grid, array, array_f
    
def load_mtz(mtz):
    """
    Loads mtz file (observed structure factors or intensities)

    Input :

	1. mtz file

    Returns :

	1. reciprocalspaceship object (effectively a pandas dataframe)

    """
    dataset = rs.read_mtz(mtz)
    dataset.compute_dHKL(inplace=True)
    
    return dataset
    
def scale(data1, data2) :

    """
    Scales the second dataset to the first. 
    Uses non-linear least squares to fit y = mx function to data (scipy.optimize.curve_fit gave best results).

    Input :

	1. dataset1 in form of 1D numpy array
	2. dataset2 in form of 1D numpy array

    Returns :

	1. the m coefficient from the fit
	2. scaled dataset2 in the form of a 1D numpy array 
    

    """
    
    def func(x, m): return m*x
    m = opt.curve_fit(func, data1, data2)[0][0]
    
    scaled = data2/m

    return m, scaled
    
def compute_weights(df, sigdf, alpha):
    
    """
    Compute weights for each structure factor based on deltaF and its uncertainty

    Input :

	1. 1D numpy array of difference structure factors
	2. 1D numpy array of errors for the difference structure factors
	3. value of alpha to be used in the qweighting (as a float)

    Returns :

	1. 1D numpy array of 1/weights (to be multiplied with the array of deltaFs for weighting)

    """
    w = (1 + (sigdf**2 / (sigdf**2).mean()) + alpha*(df**2 / (df**2).mean()))
    return w**-1
    
def solvent_mask(pdb, sampling, map_array) :

    """
    Creates a solvent mask to an electron density map

    Input :

	1. path to pdb structure to be used for solvent mask (as a string)
	2. spacing to be used to make the new grid (this should match the grid dimensions of the map you want to eventually mask)
    3. 3D map array to be masked

    Returns :

	1. flattened numpy mask array

    Parameters rprobe, rshrink, atomic radius should be optimized for purpose

    """
    
    st = gm.read_structure(pdb)
    solventmask = gm.FloatGrid()
    solventmask.setup_from(st, spacing=sampling)
    solventmask.set_unit_cell(gm.UnitCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]))
    masker.put_mask_on_float_grid(solventmask, st[0])
    nosolvent = np.where(np.array(solventmask)==0, map_array, 0)
    
    return nosolvent.flatten()
    
def res_cutoff(df, h_res, l_res) :

    """
    Given a dataset of intensities (reciprocal spaceship dataframe), applies specified low and high resolution cutoffs

    """

    df = df.loc[(df['dHKL'] >= h_res) & (df['dHKL'] <= l_res)]
    
    return df
    
def get_mapmask(grid, position, r) :
    
    """
    Returns spherical mask of a map (grid element) : radius 'r' (float) with center 'position' (vector of xyz coordinates)

    """

    grid
    grid.fill(0)
    grid.set_points_around(gm.Position(position[0], position[1], position[2]), radius=r, value=1)
    grid.symmetrize_max()
    
    return np.array(grid, copy=True)

def save_map(path, name, data) :

    """
    Writes and saves a CCP4 map using GEMMI functions.
    
    Input :
    
	1. path where output map should be saved (as string)
	2. name to identify map (as string)
	3. 3D numpy array of map (grid spacing should match grid specified in function)

    Returns :
   
	(1) CCP4 map as a GEMMI object

    """

    og = gm.Ccp4Map()
    
    og.grid = gm.FloatGrid(np.zeros((grid_size[0], grid_size[1], grid_size[2]), dtype=np.float32))
    og.grid.set_unit_cell(gm.UnitCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]))
    
    for point in tqdm(og.grid):
        u = og.grid.get_point(point.u, point.v, point.w).u
        v = og.grid.get_point(point.u, point.v, point.w).v
        w = og.grid.get_point(point.u, point.v, point.w).w
        og.grid.set_value(point.u, point.v, point.w, data[u,v,w])

    og.grid.spacegroup = gm.find_spacegroup_by_name(space_group)
    og.grid.symmetrize_max()
    og.update_ccp4_header()
    og.write_ccp4_map('{path}{name}.ccp4'.format(path=path, name=name))
    
    return og
    
def make_Nbg_map(query_map_data, background_map_data, Nbg_value):

    """Calculate query_map_data - Nbg*ref_map_data (from Panddas)"""

    Nbg_map_data = (
        (query_map_data) -
        (background_map_data * Nbg_value)
        )
    
    return Nbg_map_data
    
def save_Nbg_map(map_data, Nbg, path, name) :

    """ Saves background subtracted map calculated by make_Nbg_map function (params specific for SACLA Cl-rsEGFP2 2021 data)"""
    
    og = gm.Ccp4Map()
    og.grid = gm.FloatGrid(np.zeros((grid_size[0], grid_size[1], grid_size[2]), dtype=np.float32))
    og.grid.set_unit_cell(gm.UnitCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]))
    map_data = map_data.reshape(grid_size[0], grid_size[1], grid_size[2])
    
    for point in tqdm(og.grid):
        u = og.grid.get_point(point.u, point.v, point.w).u
        v = og.grid.get_point(point.u, point.v, point.w).v
        w = og.grid.get_point(point.u, point.v, point.w).w

        og.grid.set_value(point.u, point.v, point.w, map_data[u,v,w])

    og.grid.spacegroup = gm.find_spacegroup_by_name(space_group)
    og.grid.symmetrize_max()
    og.update_ccp4_header()
    og.write_ccp4_map('{path}/{name}_2mFo-DFc_Nbg{Nbg}.ccp4'.format(path=path, name=name, Nbg=np.round(Nbg, decimals=3)))
    
def CC_vs_n(on_map, off_map) :

    """ Obtains Pearson correlation coefficient for two maps (input at 1D numpy arrays)"""
        
    CC = np.corrcoef(on_map, off_map)[0,1]
    
    return CC
    
def CC_vs_n_2FoFc(on_map, off_map, fractions) :
        
    CCs = [np.corrcoef(make_Nbg_map(on_map, off_map, n), off_map)[0,1] for n in fractions]
    
    return CCs

def screen_qweighted(on_s, off_s, sig_on, sig_off, calc_df, Nbg, name, alpha, h_res, l_res, path, chrom_center, chrom_radius, sampling) :

    """

    For a specific Nbg value:
  	FIRST calculates the background subtracted, qweighted map given dark and light structure factors
 	THEN computes the difference between local and global correlation coefficients between that map and the dark map (by calling the get_corrdiff function)

    Input :

	1.     light structure factors (as 1D numpy array)
	2.     dark structure factors (as 1D numpy array)
	3.     light structure factor errors (as 1D numpy array)
	4.     dark structure factor errors (as 1D numpy array)
	5.     structure factors calculated from dark model (loaded as reciprocal spaceship dataframe e.g. with load_mtz function)
	6.     Nbg value to be used (as float)
	7.     name for this map iteration (as a string)
	8-9.   high and low resolution cutoffs (as floats)
	10.    path (as string) where maps should be saved. The bash script make_ccp4map.sh should be copied here.
        11-13. chromophore center (xyz coordinate vector) and radius (float) to be used for the mask, as well as grid sampling (float)

    Returns :

	1.   diff = difference in local and global correlation coefficients (float)
	2.   Nbg used
	3/4. correlation coefficients computed for the local and global regions respectively (floats) 

    """

    diffs = on_s - Nbg * off_s
    sig_diffs = np.sqrt(sig_on**2 + sig_off**2)
    
    ws = compute_weights(diffs, sig_diffs, alpha=alpha)
    diffs_w = ws * diffs
    
    calc_df["DF"]	= diffs

    calc_df["DF"]	= calc_df["DF"].astype("SFAmplitude")

    calc_df["WDF"]	= diffs_w

    calc_df["WDF"]	= calc_df["WDF"].astype("SFAmplitude")

    calc_df["SIGDF"]	= sig_diffs
    calc_df["SIGDF"]	= calc_df["SIGDF"].astype("Stddev")

    calc_df.write_mtz("{path}/{name}_diffmap_{alpha}_{Nbg}.mtz".format(path=path, name=name, alpha=alpha, Nbg=np.round(Nbg, decimals=3)))
    Nbg_map  = map_from_Fs("{path}/{name}_diffmap_{alpha}_{Nbg}.mtz".format(path=path, name=name, alpha=alpha, Nbg=np.round(Nbg, decimals=3)), "WDF", "PHI_D", map_res)
    dark_map = map_from_Fs("{path}/{name}_diffmap_{alpha}_{Nbg}.mtz".format(path=path, name=name, alpha=alpha, Nbg=np.round(Nbg, decimals=3)), "FC_D", "PHI_D", map_res)

    diff, CC_loc, CC_glob = get_corrdiff(path, Nbg_map, dark_map, alpha, Nbg, chrom_center, chrom_radius, sampling)
    
    return diff, Nbg, CC_loc, CC_glob
    
def get_corrdiff(loop_path, Nbg_map, dark_map, alpha, Nbg, chrom_center, chrom_radius, sampling) :

    """

    FIRST applies solvent mask to dark and background subtracted maps for a specific Nbg value. 
    THEN applies a chromophore mask around specified region
    RETURNS difference between local and global correlation coefficients (real-space part of analysis)

    """

    off_a             = np.array(dark_map.grid)
    on_a              = np.array(Nbg_map.grid)
    on_nosolvent      = np.nan_to_num(solvent_mask(pdb_path, sampling, on_a))
    off_nosolvent     = np.nan_to_num(solvent_mask(pdb_path, sampling, off_a))
    mask              = get_mapmask(Nbg_map.grid, chrom_center, chrom_radius)
   
    ### optional save this mask as a map ##
    #masked_map               = save_map(loop_path, name, mask)
   
    chrom 	= np.array(mask, copy=True).flatten().astype(bool)
    CC_loc 	= CC_vs_n(on_a.flatten()[chrom], off_a.flatten()[chrom])
    CC_glob	= CC_vs_n(on_nosolvent[np.logical_not(chrom)], off_nosolvent[np.logical_not(chrom)])
    
    diff 	= np.array(CC_glob) -  np.array(CC_loc)
    
    return diff, CC_loc, CC_glob
    

    
    
