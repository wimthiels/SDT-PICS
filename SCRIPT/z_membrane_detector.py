'''
Created on 1 nov 2019 
extracts all cells from a vtp based on parentID
		
@author: wimth

# debug output files : 
# 5_ddzsmooth_peak_overlay.tif : shows you the all the detected peaks for the stack (already filtered on distance to membrane) (control by peak_thresh_rel_slice)
# 6_bool_peak_img.tif : check to see the peaks after the S/N interval cutoff  (control with SN_thresh), can be bypassed via (z_metric_mode )
# 7_extended a_peaks_extend_over_original : check to see the peaks after extension (control by peak_thresh_rel_slice) = final result : check this first !
'''
import numpy as np
import re,os
import pandas as pd
from param_XML import Param_xml
from tifffile import TiffFile,imshow,imsave
from skimage.filters import gaussian,median
from skimage.morphology import reconstruction
from skimage.morphology import disk
from skimage.filters.rank import maximum, minimum, autolevel
from skimage.draw import ellipsoid,circle
# from skimage.metrics import normalized_root_mse, mean_squared_error, peak_signal_noise_ratio  #valid for skimage v17 and up
from skimage.measure import compare_nrmse, compare_psnr, compare_mse

from pathlib import Path
from helper_functions import examine
from scipy.ndimage.morphology import distance_transform_edt
from skimage.feature import peak_local_max

import os,sys


def detect_z_membrane(a_input,a_membrane_mask,a_membrane,a_exterior_mask,
	l_spacing=[1,0.1083474,0.1083474],sigma_ddz_smoothing = 4,min_dist_peak_membrane = 5,
	min_dist_between_peaks=3,peak_thresh_rel=0.4,peak_thresh_rel_slice=0,SN_thresh=2,z_metric_mode="top-bottom",
	path_output_folder=None,verbose=False):
	""" 
	in : 
		- a_input : nd_array of the original input
		- a_membrane_overlay : an that can hide the membrane already detected = the bright membrane.  (> 0 indicates membrane pixel) 
		- sigma_ddz_smoothin : smoothing of the ddz values = d_gaussian2
		- min_dist_peak_membrane : the minimum distance of a peak to a membrane pixels.  Used to give peaks near the center of the membrane

	out: 
		nd-array : boolean: True = location of a z-membrane

	This module will try to detect parallel z-membranes
	- the strong signal that is already detected is removed (=the binary membrane image after thresholding the frangi (XY) output
	- makes use of an extended membrane mask i.e.  the all pixels above and below the membrane is added (this is because these pixels will have an incorrect z-gradient)
	- no z-smoothing is done ! This is because a membrane can occur over 1 Z-pixel, so we do not want to blur this dimension.  We do want to aggregate in XY to share plane info
	- because we only want parallel Z planes we only look at ddz info. (and ) ddx and ddy information turned out not to be informative)
	"""

	d_gaussian_1= {'sigma': [0,1,1],'mode': 'nearest','preserve_range': True,'truncate': 4.0 }
	d_gaussian_2= {'sigma': [0,sigma_ddz_smoothing,sigma_ddz_smoothing],'mode': 'nearest','preserve_range': True,'truncate': 4.0 }
	d_peak_local_max = {'min_distance':min_dist_between_peaks, 'threshold_abs':None, 'threshold_rel':peak_thresh_rel, 'exclude_border':True, 'indices':False, 'num_peaks' :10000,'footprint':None, 'labels':None, 'num_peaks_per_label':10000}
	d_reconstruction = {'method':'dilation'}
	# d_peak_signal_noise_ratio = {'image_true':None,'image_test':None,'data_range':None}
	d_peak_signal_noise_ratio = {'im_true':None,'im_test':None,'data_range':None}
	XY_VS_Z = l_spacing[0]/l_spacing[1]
	if path_output_folder:
		path_output_folder.mkdir(parents=True,exist_ok=True)
	d_log = {};dlog1={};dlog2={};dlog3={}

	#step 1 : remove strong membrane signal-----------------------------------------------------------------
	a_membrane_overlay = np.where(a_membrane_mask>0,0,a_input) #fuse membrane (problem : leaves much bright signal near the membrane)
	# a_membrane_overlay = np.where(a_input_file_frangi==0,a_input,0) #only keep region with 0 frangi ! not smart, because will also kill some vague z-membrane signal

	a_extended_mask = np.where(a_membrane_overlay>0,0,1)
	ix_z_max = a_membrane_overlay.shape[0] - 1 
	for ix_z, a_slice in enumerate(a_membrane_overlay):
		if ix_z>0:
			a_extended_mask[ix_z - 1] = np.where(np.logical_and(a_slice>0,a_extended_mask[ix_z - 1]==1),0,a_extended_mask[ix_z - 1])
		if ix_z<ix_z_max:
			a_extended_mask[ix_z + 1] = np.where(np.logical_and(a_slice>0,a_extended_mask[ix_z + 1]==1),0,a_extended_mask[ix_z + 1])	
	a_original_extended_mask = np.where(a_extended_mask==1,a_input,0) #will not be used_just for vizualization

	if path_output_folder:
		imsave(path_output_folder/"1_membrane_overlay.tif", a_membrane_overlay.astype('int16'))
		imsave(path_output_folder/"1_extended_mask.tif", a_extended_mask.astype('int16'))
		imsave(path_output_folder/"1_original_extended_mask.tif(not used)", a_original_extended_mask.astype('int16'))

	#step 2 : smooth in XY----------------------------------------------------------------------------
	a_smoothedXY = gaussian(a_original_extended_mask, **d_gaussian_1)
	a_smoothedXY_filled = np.where(a_extended_mask==0,a_input,a_smoothedXY)
	if path_output_folder:
		imsave(path_output_folder/"2_smoothedXY.tif", a_smoothedXY.astype('uint16'))
		imsave(path_output_folder/"2_smoothedXY_filled.tif", a_smoothedXY_filled.astype('uint16'))


	#step 3 : gradient dzdz, dydy, dxdx-----------------------------------------------------------
	a_ddz = np.gradient(np.gradient(a_smoothedXY_filled,l_spacing[0],axis=0),l_spacing[0],axis=0)
	a_ddz = np.where(a_ddz>0,0,a_ddz*(-1)) #we only want negative second derivatives ; originale range (the max is  3.3 , the min is  -4.2)
	if path_output_folder:
		a_ddz_viz = np.where(a_extended_mask==0,0,a_ddz) #range is now : the max is  3.1695729534251993 , the min is  0.0
		a_ddz_viz = a_ddz_viz*100 #just to create some range for visualization
		imsave(path_output_folder/"3_dzdz_extended_mask.tif", a_ddz_viz.astype('int16'))

	# #step 4 : smooth ddz over XY-----------------------------------------------------------
	a_ddz_overlay = np.where(a_extended_mask==0,0,a_ddz)
	a_ddz_smoothedXY_overlay = gaussian(a_ddz_overlay, **d_gaussian_2)
	if path_output_folder:
		imsave(path_output_folder/"4_ddz_smoothedXY_overlay.tif", (a_ddz_smoothedXY_overlay*100).astype('int16'))

	# step 5 : find local maxima
	#  first create mask so that peaks are only allowed at some distance from the membrane
	a_dist_transf = distance_transform_edt(a_membrane_overlay==0, sampling = [XY_VS_Z,1,1], return_distances=True, return_indices=False) 
	a_dist_transf_membrane = distance_transform_edt(a_membrane==0, sampling = [XY_VS_Z,1,1], return_distances=True, return_indices=False) #added for scoring
	a_labels = a_dist_transf>min_dist_peak_membrane

	d_peak_local_max['labels']=a_labels

	use_masked_peak_image=True #True is default
	if use_masked_peak_image:
		if verbose:print('using the masked image for peak detection..')
		a_ddz_smoothedXY_overlay_masked = np.where(a_exterior_mask==0,0,a_ddz_smoothedXY_overlay)
		a_bool_peak_img = peak_local_max(a_ddz_smoothedXY_overlay_masked,**d_peak_local_max)
		d_peak_local_max['indices']=True
		a_l_indices = peak_local_max(a_ddz_smoothedXY_overlay_masked,**d_peak_local_max)
		dlog3['peak_threshold'] = np.max(a_ddz_smoothedXY_overlay_masked)*peak_thresh_rel
		a_ddzsmooth_peak_overlay = a_ddz_smoothedXY_overlay_masked.copy()
		if verbose:print("{0} is the max of the preprocessed stack so peak threshold for the stack will be set at {1}".format(np.max(a_ddz_smoothedXY_overlay_masked),np.max(a_ddz_smoothedXY_overlay_masked)*peak_thresh_rel))
	else:
		a_bool_peak_img = peak_local_max(a_ddz_smoothedXY_overlay,**d_peak_local_max)
		d_peak_local_max['indices']=True
		a_l_indices = peak_local_max(a_ddz_smoothedXY_overlay,**d_peak_local_max)
		dlog3['peak_threshold'] = np.max(a_ddz_smoothedXY_overlay)*peak_thresh_rel
		a_ddzsmooth_peak_overlay = a_ddz_smoothedXY_overlay.copy()
		if verbose:print("{0} is the max of the preprocessed stack so peak threshold for the stack will be set at {1}".format(np.max(a_ddz_smoothedXY_overlay),np.max(a_ddz_smoothedXY_overlay)*peak_thresh_rel))

	#  visualize
	if verbose:print("{0} peaks detected".format(a_l_indices.shape[0]))
	a_original_peak_overlay = a_input.copy()
	z,y,x = a_original_peak_overlay.shape
	for peak_idx in a_l_indices:
		rr,cc = circle(peak_idx[1],peak_idx[2],radius=4,shape =(y,x))
		a_ddzsmooth_peak_overlay[peak_idx[0],rr,cc] = 0
		a_original_peak_overlay[peak_idx[0],rr,cc] = 0

	if path_output_folder:
		imsave(path_output_folder/"5_dist_transf.tif", (a_dist_transf).astype('int16'))
		imsave(path_output_folder/"5_labels.tif", (a_labels).astype('int16'))
		imsave(path_output_folder/"5_ddzsmooth_peak_overlay.tif", (a_ddzsmooth_peak_overlay*100).astype('int16'))
		imsave(path_output_folder/"5_original_peak_overlay.tif", a_original_peak_overlay.astype('int16'))

	# step 6 : remove peaks in lower layers where signal deteriorates (=heuristic 1)
	a_z_metric = np.sum(np.where(a_bool_peak_img,a_dist_transf,0),axis=(1,2))
	#a_z_metric = np.sum(np.where(a_bool_peak_img,a_dist_transf_membrane,0),axis=(1,2)) #update 0622 : for scoring we need to original membrane, otherwise some peaks will get filtered where no membrane is present
	l_z_metric = list(a_z_metric)
	dlog2['l_z_metric'] = l_z_metric
	if verbose:print("a_z_metric = {0}".format(a_z_metric))

	#           determine z-metric peaks from both sides
	ix_min = 0
	z_metric_high = 0
	
	for ix_z, z_metric_i in enumerate(l_z_metric):
		if z_metric_i < z_metric_high:
			break
		else:
			z_metric_high=z_metric_i
			ix_min = ix_z


	ix_max = len(l_z_metric) - 1
	z_metric_high = 0
	a_nb_peaks = np.sum(a_bool_peak_img,axis=(-2,-1))

	for ix in range(len(l_z_metric)-1,-1,-1):
		if l_z_metric[ix] < z_metric_high:
			ix_max = ix +1
			if a_nb_peaks[ix_max] < 2 and ix_max +1 < len(a_nb_peaks):
				ix_max += 1
			break
		else:
			z_metric_high=l_z_metric[ix]
 

	# extra metric : track S/N of every slice (use as switch) (=heuristic 2)
	a_membrane_mask = np.where(a_membrane_mask>0,a_input,0)  #=the signal (membrane)  in the original image
	if path_output_folder:
		imsave(path_output_folder/"6_membrane_mask.tif", a_membrane_mask.astype('int16'))
	l_MS_signal = []
	l_MS_noise = []
	for ix_z in range(a_membrane_mask.shape[0]):
		slice_signal = a_membrane_overlay[ix_z]
		slice_noise  = a_membrane_mask[ix_z]
		#l_MS_signal.append(np.square(slice_signal[slice_signal>0]).mean())
		l_MS_signal.append(np.mean(np.square(slice_signal[slice_signal>0])))
		l_MS_noise.append(np.mean(np.square(slice_noise[slice_noise>0])))

	l_SN_ratio = [signal/noise for signal,noise in zip(l_MS_signal,l_MS_noise)]
	dlog1['SN_ratio'] = l_SN_ratio
	if l_SN_ratio[ix_min] < SN_thresh:  # we use the first peak slice as a proxy for overall stack quality.  If S/N ratio of this slice drops below 2, we switch of z-membrane completely
		ix_max = ix_min
		if verbose:print("z_membrane switched OFF for this stack, signal quality stack is to poor:")

	if z_metric_mode=='top-bottom':
		a_bool_peak_img[:ix_min+1,...] = 0
		a_bool_peak_img[ix_max:,...] = 0
		if verbose:print('ix_min={2}, ix_max = {3}. top-bottom mode : Peaks are only kept in index-interval [{0},{1}]'.format(ix_min+1,ix_max-1,ix_min,ix_max))
	elif z_metric_mode =='bottom':
		a_bool_peak_img[ix_max:,...] = 0
		if verbose:print('ix_min={2}, ix_max = {3}. bottom mode :  Peaks are only kept in index-interval [{0},{1}]'.format(ix_min+1,ix_max-1,0,ix_max))
		ix_min = 'NA'
	else:
		if verbose:print('z-metric not used')
		ix_min = 'NA'
		ix_max = 'NA'
	dlog2['l_z_metric'].append(ix_min)
	dlog2['l_z_metric'].append(ix_max)
	dlog2['l_z_metric'].append(a_l_indices.shape[0]) #nb of peaks
		
	if path_output_folder:
		imsave(path_output_folder/"6_bool_peak_img.tif", a_bool_peak_img.astype('int16'))
		

	#step 7 : expand the peaks to a region with a downhill filter
	def filter_reconstruction(img, ix_slice):
		'''
		This is the downhill filter.  This reconstruction is variable ! it will only be done for thresholds that fall below the minimum
		it expects as input the frangi output
		'''
		thresh = a_abs_thresh_slice[ix_slice]
		thresh_seed = a_abs_thresh_slice_seed[ix_slice]

		seed = img.copy()
		seed[seed < thresh_seed] = 0
		img = reconstruction(seed=seed, mask=img, **d_filter_arg['reconstruction'])
		#             print('reconstruction applied to slice {0} for threshold {1} and thresh_seed {2}'.format(ix_slice + 1,thresh,thresh_seed))

		return img

	
	a_peaks_extended = a_bool_peak_img.copy()
	a_seed = np.where(a_bool_peak_img,a_ddz_smoothedXY_overlay,0)
	if use_masked_peak_image:
		a_iter = a_ddz_smoothedXY_overlay_masked
	else:
		a_iter = a_ddz_smoothedXY_overlay
	for ix_z, a_slice_i in enumerate(a_iter):
		# a_slice_i = np.where(a_bool_peak_img[ix_z], np.max(a_slice_i) + 1, a_slice_i) # Intensity of seed image must be less than that of the mask image for reconstruction by dilation.
		if np.sum(a_seed[ix_z])==0:  #no peaks or peaks with zero value
			a_peaks_extended[ix_z] = a_bool_peak_img[ix_z]
			dlog3.setdefault("peak threshold slice", []).append("no peaks or peaks with zero value")
			dlog3.setdefault("max value slice(=ref)", []).append(np.max(a_slice_i))
		else:
			if 	peak_thresh_rel_slice>0:
				if verbose:print("threshold_extension logic active :{0}".format(peak_thresh_rel_slice))
				# threshold_extension = np.max(a_slice_i) *( 2/3) # no exact science, 2/3 seems ok
				threshold_extension = np.max(a_slice_i) *(peak_thresh_rel_slice) #
			else:
				threshold_extension = np.min(a_seed[ix_z][np.nonzero(a_seed[ix_z])]) #Iimportant, not trivial, the threshold will be extended to the lowest peak detected. this way all peaks will be included and the lowest peak will be 1 pixel
				threshold_extension *= 0.98 #give the lowest peak a little spread (no exact science..)
			# print('temp: threshold extension', threshold_extension, 'for z', ix_z)
			dlog3.setdefault("peak threshold slice", []).append(threshold_extension)
			dlog3.setdefault("max value slice(=ref)", []).append(np.max(a_slice_i))
			seed = a_slice_i.copy()
			seed[seed < threshold_extension] = 0
			seed[a_bool_peak_img[ix_z]==0 ]= 0
			a_slice_i[a_slice_i < threshold_extension] =0
			# a_peaks_extended[ix_z] = reconstruction(seed=a_seed[ix_z], mask=a_slice,**d_reconstruction)
			a_peaks_extended[ix_z] = reconstruction(seed=seed, mask=a_slice_i,**d_reconstruction)

	if path_output_folder:
		a_peaks_extend_over_original = np.where(a_peaks_extended > 0, np.max(a_input), a_input)
		imsave(path_output_folder/"7_extended a_peaks_extend_over_original.tif", a_peaks_extend_over_original.astype('int16'))


	dflog={'z_SN_ratio':pd.DataFrame(dlog1),'z_metric':pd.DataFrame(dlog2), 'z_peak_threshold_slice':pd.DataFrame(dlog3,columns = ['peak threshold slice', 'peak_threshold', 'max value slice(=ref)' ])}
	return a_peaks_extended,dflog


if __name__== "__main__":

	#read files
	input_file_original      = "/home/wth/Downloads/testinfra/ISOLATED/z_membrane_detector/TL11_crop_1ch_emb1of1.tif"
	input_file_membrane_mask = "/home/wth/Downloads/testinfra/ISOLATED/z_membrane_detector/threshold(acceptance_level=10).tif"
	input_file_frangi        = "/home/wth/Downloads/testinfra/ISOLATED/z_membrane_detector/frangi.tif"
	output_folder            = Path("/home/wth/Downloads/testinfra/ISOLATED/z_membrane_detector/")
	ix_stack = 0

	with TiffFile(input_file_original) as tif:
		a_input =tif.asarray()[ix_stack,:].astype(float) #NP.GRADIENT ONLY WORKS WITH FLOATS !!!!!

	with TiffFile(input_file_membrane_mask) as tif:
		a_input_file_membrane_mask =tif.asarray()[ix_stack,:].astype(float)

	with TiffFile(input_file_frangi) as tif:
		a_input_file_frangi =tif.asarray()[ix_stack,:].astype(float)

	detect_z_membrane(a_input,a_input_file_membrane_mask,path_output_folder=output_folder,verbose=True)

	print("END")  