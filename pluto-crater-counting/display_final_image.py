#! /usr/bin/env python

"""Reads in the flat map image and the all the crater coordinate
files and plots the crater coordinates over the flat map.

This is to be both sanity check and a final image to be displayed
in the presentation.

Author:

    C.M. Gosmeyer, Nov. 2015
    
Use:
    
    >>> python display_final_image.py
    
Outputs: 

"""

# Need make sure that in the sub-sections the proper number of
# x and y pixels are added from the 0,0, defined as bottom corner.
# (Note that in the uncertainty calculation, 0,0 is defined
# as the center of the body.)
# -- all this info should be in the titles of the text files.

import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import ascii
from skimage import data, color
from skimage.io import imread


path = '/path/to/main/dir/containing/outputs/data/and/scripts/dirs/'
# Read in all the given body's crater count files.

# For a given body, the scaling should probably be the same in all the 
# sub-sections. The task then is to add to each crater's 
# xdata, ydata the sub-section's displacement, and to scale the radius
# (pixels) to the full image.  
# How do pixels in the original image transform to meters?

# For each region,
# Plot the xdata, ydata in color for given region (as sanity check)

# Save figure.

# Note that if matplotlib is going to be a jerk about saving the full-size
# image I'll have to find another package or do this step in IDL. 


def write_to_dict(data, data_dict, region, xlower, ylower):
	"""
	"""
    # First check whether region is already a key.
    # If not, create new key with 6 emtpy list values.
    if region not in dict.keys():
    	dict[region] = [[], [], [], [], [], []]

    for i in range(len(data)):

        # The displaced x and y need be added to xdata,ydata.
	    # xdata
        (data_dict[region])[0].append((data[0])[i] + xlower)

        # ydata
        (data_dict[region])[1].append((data[1])[i] + ylower)

        # radius
        (data_dict[region])[2].append((data[2])[i])

        # xcanvas
        (data_dict[region])[3].append((data[3])[i])

        # ycanvas
        (data_dict[region])[4].append((data[4])[i])

        # scalar
        (data_dict[region])[5].append((data[5])[i])

	return data_dict


def loop_through_files(body):
	"""
	"""
	count_files = glob.glob(path+'outputs/'+body+'*.txt')

    # Initialize dictionaries whose keys will be regions. Use 'north' as 
    # initializer because it exists in both pluto and charon data.
    # THe values will be lists of xdata, ydata, radius, xcanvas, ycanvas, scaler
    data_dict = {'north': [[], [], [], [], [], []]}


	for count_file in count_files:
        # Parse the file name for region.
        region = ((count_file.split('/'))[len(count_file.split('/'))-1]).split('_')[1]

        # And for the xlower
        xlower = (((count_file.split('/'))[len(count_file.split('/'))-1]).split('_')[2]).split('-')[0]

        # And for the ylower
        xlower = (((count_file.split('/'))[len(count_file.split('/'))-1]).split('_')[2]).split('-')[1]

		# Read columns
        data = ascii.read(count_file)
        
        # Append to lists in appropriate region.
        data_dict = write_to_dict(data, data_dict, region, float(xlower), float(ylower))

	return data_dict


def plot_craters_on_body(body, data_dict):
    """
    """

    # Read in the full image.
    if body == 'pluto':
        image = path + 'data/pluto.jpg'
    elif body == 'charon':
    	image = path + 'data/charon.jpg'

    data = imread(image)
    
    # For each region (key in data_dict), plot a circle at the xdata, ydata.
    # Need scale the crater radii using the scaler? 
    
    for region in data_dict.keys():




if __name__=='__main__':
    body = 'charon'
    data_dict = loop_through_files(body)

