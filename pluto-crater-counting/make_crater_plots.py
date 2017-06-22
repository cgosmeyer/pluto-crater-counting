#! /usr/bin/env python

"""Reads in the files and creates plots.

Author:

    C.M. Gosmeyer, Nov. 2015
    
Use:
    
    >>> python make_crater_plots.py
    
Outputs: 

"""

import glob
import numpy as np
import os
from astropy.io import ascii
import matplotlib.pyplot as plt
import pylab
import pyraf
# Probably should replace the iraf linear-fitting if run this again.
from pyraf import iraf
from iraf import tables, ttools, tcreate, tlinear, tdump


path = '/path/to/main/dir/containing/outputs/data/and/scripts/dirs/'

charon_diam_pixels = 1628.
pluto_diam_pixels = 4312.




def transform_meters_to_pixels(body, meters):
    """Transforms meters to pixels using scale of the given body.

    Radii from Stern, S.A. et al. (2015) "The Pluto System:
    Initial results from its exploration by New Horizons." Science 350 (6258): 249-352
    """
    if body == 'charon':
        radius_meters = 606000.  # mean
        radius_pixels = charon_diam_pixels / 2. # +/- 5.
    elif body == 'pluto':
        radius_meters = 1187000.  # mean
        radius_pixels = pluto_diam_pixels / 2.  # +/- 5.

    pixels_per_meter = radius_pixels / radius_meters
    pixels = pixels_per_meter * meters    
 
    return pixels


def transform_pixels_to_meters(body, pixels):
    """Transforms pixels to meters using scale of the given body.

    Note this does NOT take into account the warping at the edges.
    So really must take a weight...

    Radii from Stern, S.A. et al. (2015) "The Pluto System:
    Initial results from its exploration by New Horizons." Science 350 (6258): 249-352
    """
    if body == 'charon':
        radius_meters = 606000.  # mean
        radius_pixels = charon_diam_pixels / 2. # +/- 5.
    elif body == 'pluto':
        radius_meters = 1187000.  # mean
        radius_pixels = pluto_diam_pixels / 2.  # +/- 5.

    meters_per_pixel = radius_meters / radius_pixels
    meters = meters_per_pixel * pixels       
 
    return meters


def apply_weights(body, xdist, ydist, radius):
    """Maybe. Make this easy to adjust. For now leave as 1. 
    Applies weight based on distance from body's 'center'.
    """
    if body == 'charon':
        center_coords = ''  # read in image and get the dimensions?
        body_radius = ''
        weight = 1.  # Do a math thing based on x and y dist from center

    return weight


def read_in_region_files(body, region):
    """Reads in all the files of the given region.
    Returns the coordinates and radii.

    xcoord_data, ycoord_data, radius, xcoord_canvas, ycoord_canvas, scaler
    """

    path_to_files = path + 'outputs/' + body + '/'

    # Initialize the lists.
    xcoord_data_ls = []
    ycoord_data_ls = []
    radius_ls = []
    xcoord_canvas_ls = []
    ycoord_canvas_ls = []
    scaler_ls = []


    region_files = glob.glob(path_to_files + '*' + region + '*txt')
    print region_files

    for region_file in region_files:
        data = ascii.read(region_file)

        for i in range(len(data)):
            xcoord_data_ls.append(float(data['col1'][i]))
            ycoord_data_ls.append(float(data['col2'][i]))
            radius_ls.append(float(data['col3'][i]))
            xcoord_canvas_ls.append(float(data['col4'][i]))
            ycoord_canvas_ls.append(float(data['col5'][i]))
            #scaler_ls.append(float(data['col6'][i]))
    
    return xcoord_data_ls, ycoord_data_ls, radius_ls, \
           xcoord_canvas_ls, ycoord_canvas_ls #, scaler_ls




def bin_and_transform_crater_sizes(body, radii_all, areas_all, plot_type, bin_size=.5):
    """Transforms the crater radii from pixels to meters.
    Converts radii to diameter.
    Bins the craters cumalitively by diameter. 
    Removes the bins with zero counts.


    Parameters:
        body :
        radii :
            The radii of all craters in given surface areas. (pixels)
        areas : 
            The surface area (in km2) where each radius was measured.
        bin_size : float 
            The width of bin in km. Bins by diameter. 
        plot_type:
            'R' bins by intervals.
            'cum' bins cumulatively

    """

    # Bin the cumalative craters by diameter.
    # Initialize dictionary that will have diameters as keys and 
    # counts as values.
    """ diam_dict = {} """
    
    # Find the max of all the radii. Need use the diameter converted to km
    # as bin end for ALL sub-region boxes.

    # Make radii_all and areas_all lists of lists if not.
    if type(radii_all[0]) != list:
        radii_all = [radii_all]
        areas_all = [areas_all]

    # Find the max so will have a constant upper bin limit for region.
    radii_all_flat = [val for sublist in radii_all for val in sublist]
    max_diam_pix = np.max(radii_all_flat) * 2.
    max_diam_km = transform_pixels_to_meters(body, max_diam_pix)/1000.
    print "Max in pix and km: ", max_diam_pix, max_diam_km

    # For sake of getting uncertainty.
    areas_for_uncert = [val for sublist in areas_all for val in sublist]
    max_areas = np.max(areas_for_uncert)

    # Loop through each area. In end will append numbers together?
    counts_all = []
    # Initialize count (not normalized to area) that will use to calc error, N^.2
    counts_nonnorm = 0
    # Initialize list for ALL D values, to be plotted in R plot.
    D_all = []



    for radii, areas in zip(radii_all, areas_all):
        # Convert radii to diams.
        diams_pix = (np.array(radii)) * 2.
        #print diams_pix

        # Convert diams from pixels to kms.
        diams = []
        for diam in diams_pix:
            diams.append(transform_pixels_to_meters(body, diam)/1000.)

        print diams
        # Initialize the cumulative bins.
        # The max end must be largest diam (km) where number of craters = 1.
        bins = np.arange(0, max_diam_km, bin_size)
        print bins

        # Initialize the counting list. Each entry corresponds to the number of 
        # craters in the cumalative bin.
        counts = np.zeros(len(bins)-1)

        # Initialize count of diameters to be used in R plot.
        diam_counts = np.zeros(len(bins)-1)

        # Initialize the D value list to be used in the R plots 
        # R = (D^3)*n / A[bin], D = sum(diams)^-n
        D = []

        """
        bins = np.arange(.5, bin_limit, bin_size)
        # Insert bins as keys and set values (counts) to 0.
        for bin in bins:
           # Convert pixels to meters.
           diam_dict[transform_pixels_to_meters(body, bin)] = 0


        # Sort the dictionary so that bins actually make sense.
        diam_dict = collections.OrderedDict(sorted(diam_dict.items()))
        bins_m = diam_dict.keys()
        """
    
        if plot_type == 'cum':
            # Now bin cumulatively. 
            for i in range(len(bins)-1):
                for diam in diams:
                    if (diam > bins[i]):
                        counts[i] += 1 
                        counts_nonnorm +=1
        elif plot_type == 'R':
            for i in range(len(bins)-2):
                for diam in diams:
                    if ((diam > bins[i]) and (diam <= bins[i+1])):
                        diam_counts[i] += diam
                        counts[i] += 1   
                        counts_nonnorm += 1    
                # Get the D value.
                D.append(np.sum(diam_counts[i]) ** (1./counts[i]))
            D.append(diams[len(diams)-1])



        """
        # Get dict into list form to make this easy.
        nonzeros = np.where(np.array(diam_dict.values()) != 0)[0]

        crater_bins = (np.array(diam_dict.keys()))[nonzeros]
        crater_counts = (np.array(diam_dict.values()))[nonzeros]

        for bin_i in range(len(bins_m)-2):
            for diam in diams:
                if ((diam > bins_m[bin_i]) and (diam <= bins_m[bin_i+1])):
                    diam_dict[bins_m[bin_i]] += 1
        """
        if plot_type == 'cum':
            counts_all.append(np.array(counts)/areas[0])
        elif plot_type == 'R':
            #counts_all.append(( (np.array(D)**3)*np.array(counts) ) / areas[0])
            counts_all.append(counts)
            D_all.append(D)


    if plot_type == 'cum':
        # Now add together the area-normalized counts or R values.
        counts_added = np.zeros(len(bins)-1)
        for counts_ls in counts_all:
            counts_added += np.array(counts_ls)

        # Remove the zeros for plot.
        nonzeros = np.where(np.array(counts_added) != 0)[0]

        # Below is diameters.
        crater_bins = bins[nonzeros]
        # Below is cumulative counts.
        
        crater_counts = np.array(counts_added[nonzeros]) 
    
    if plot_type == 'R':
        # Calculate the D value for each bin over all areas.
        D_added = np.zeros(len(bins)-1)
        for D_ls, counts_ls, areas in zip(D_all, counts_all, areas_all):
            D_added += (np.array(D_ls)**3) * np.array(counts_ls) / areas[0]
        # Remove the zeros for plot.
        nonzeros = np.where(np.array(D_added) != 0)[0]

        # Below is either geom mean of diameters (D).
        crater_bins = bins[nonzeros]
        # Below is either R.
        crater_counts = np.array(D_added[nonzeros]) 


    if plot_type == 'cum':
        print "Total craters in region: ", counts_nonnorm
        print "Error for region's counts in counts^.5: ", np.sqrt(counts_nonnorm)
        print "Uncertainty in the first bin is counts^0.5 / max(area): ", (np.sqrt(counts_nonnorm) / max_areas)
        print "Uncertainty in cumumlative plots in each bin is (counts/area)/N^.5"
        print (np.array(crater_counts) / np.sqrt(crater_bins))
    elif plot_type == 'R':
        print "Total craters in region: ", counts_nonnorm
        print "Uncertainty in the R value is R[i] / N^.5: "
        print crater_counts / np.sqrt(counts_nonnorm)

    return crater_bins, crater_counts


def return_color(body, region):
    """
    """
    if body == 'charon':
        if 'belt' in region:
            color = 'green'
        elif 'harad' in region:
            color = 'orange'
        elif 'mordor' in region:
            color = 'r'
        elif 'north' in region:
            color = 'b'
        elif 'south' in region:
            color = 'violet'

    elif body == 'pluto':
        if 'south' in region:
            color = 'indigo'
        elif 'rbbelt' in region:
            color = 'magenta'
        elif 'bbelt' in region:
            color = 'red'
        elif 'ltbelt' in region:
            color = 'brown'
        elif 'tbelt' in region:
            color = 'pink'
        elif 'lheart' in region:
            color = 'lime'
        elif 'rheart' in region:
            color = 'green'
        elif 'lamp' in region:
            color = 'orange'
        elif 'snakeskin' in region:
            color = 'lightgrey'
        elif 'rnorthterm' in region:
            color = 'skyblue'
        elif 'rnorth' in region:
            color = 'cyan'
        elif 'north' in region:
            color = 'blue'

    return color

# -----------------------------------------------------------------------------
# line-fitting.
# -----------------------------------------------------------------------------


def plot_linear_fit(fitdata_file, plot_on=True, scol=True):
    """Plots the linear fit and prints slope on plot.
    
    Parameters: 
        fitdata_file : text file
            Assumed to be ``<filter>_<amp>_tlinear.dat``, containing
            the linear data found from ``IRAF/TLINEAR``.
        plot_on : {True, False}
            True by default. To turn off line plotting, set to False.
        scol : {True, False}
            True by default. Set False if don't want to weight the fit
            by the deltaflux error.
    
    Returns:
        slope : float
            The slope of the linear fit.
            
    Outputs:
        A linear fit overplotted on a plot containing the points that
        were fit with ``IRAF/TLINEAR``. (Plot must be created before 
        this function is called.) 
    """
    # Read contents of '<filter>_<amp>_tlinear.dat' to lists.
    xfit = []
    yfit = []
    with open(fitdata_file) as fitdata:
        for line in fitdata.readlines()[0:]:
            line = line.split()
            xfit.append(float(line[0]))  #0
            if scol:
                yfit.append(float(line[2]))  #3
            else:
                yfit.append(float(line[1]))
        
        
        # Find slope of line:
        try:
            slope = (yfit[len(yfit)-1] - yfit[0]) / (xfit[len(xfit)-1] - xfit[0])
        except:
            print "The start and end dates match -- no slope calculable."
            slope = 0.0
        
        if plot_on:
            pylab.plot(xfit, yfit, c='black', linewidth=2.5)
#            pylab.figtext(0.91, 0.28 + plotnum, amp + ' slope:', size='x-small')
#            pylab.figtext(0.91, 0.25 + plotnum, str(slope), size='x-small')
            # Prints the slope in the lower left-hand corner, rounded to 6 s.d.
            pylab.annotate("slope: " + str(round(slope, 6)),[0.03,0.08], \
                           xycoords='axes fraction', fontname='monospace', size='large')
                       
    return slope


#-------------------------------------------------------------------------------#

def run_tcreate(y, x, yerror, region, loc):
    """Creates an ``IRAF`` table from delta-flux vs time points.
    
    This intermediate step is necessary in order to run ``IRAF/TLINEAR``
    because, as far as we know, it only takes ``IRAF`` tables.
    
    Parameters: 
        deltaflux : list of floats
            The delta-fluxes calculated for a given filter, detector,
            and amp.
        time : list of floats
            The times corresponding to the delta-fluxes.
        deltafluxerr : list of floats
            The formal error in deltaflux.
        rootname : list of strings
            Rootnames of files. 
        destination : string
            Where you want the files generated.
            
    Returns:
        nothing        
            
    Outputs:
        ``IRAF`` table file named ``tcreate.tab``.
    """
    # First create the input ascii file.
    cols = {'#x':x, 'y':y, 'yerror':yerror}
    print y
    print x
    ascii.write(cols, loc + 'fluxvstime.dat', \
                names=['#x', 'y', 'yerror'], \
                formats={'#x':'%9.4f', 'y':'%0.4f', 'yerror':'%0.9f'})
    
    # Second create the column definition file.
    name = ['x', 'y', 'yerror']
    datatype = ['r', 'r', 'r']
    decimals= ['f9.7', 'f9.7', 'f0.9']
    cols = {'#colname':name, 'datatype':datatype, 'decimals':decimals}
    ascii.write(cols, loc + 'coldef.dat', \
                names=['#colname', 'datatype', 'decimals'], \
                formats={'#colname':'%s', 'datatype':'%s', 'decimals':'%s'})
    
    # Finally run IRAF's tcreate.
    iraf.tcreate(loc+'tcreate.tab', loc+'coldef.dat', loc+'fluxvstime.dat')
    
    # Delete the column definition and input files.
    os.remove(loc+'coldef.dat')
    os.remove(loc+'fluxvstime.dat')


#-------------------------------------------------------------------------------#

def run_tlinear(loc, region, scol=True):
    """Runs ``IRAF/TLINEAR`` (linear fit) on data in file ``tcreate.tab``.
    
    Parameters:
        scol : {True, False}
            True by default. Set False if don't want to weight the fit
            by the deltaflux error.
           
    Returns:
        nothing       
                    
    Outputs:
        text file. The linear fit stored in 
        ``<filtername>_<amp><outfilemark>_tlinear.dat``.
        text file. The parameters of the fit stored in
        ``<filtername>_<amp><outfilemark>_tlinear.hdr.dat``
    """
    if scol:
        iraf.tlinear(loc+'tcreate.tab', loc+'tlinear.tab', \
                     xcol='x', ycol='y', scol='yerror' )
    else:
        iraf.tlinear(loc+'tcreate.tab', loc+'tlinear.tab', \
                     xcol='x', ycol='y')
        
    # Convert IRAF table to ascii file.
    iraf.tdump(loc+'tlinear.tab', \
                datafile = loc + region + '_tlinear.dat', \
                pfile = loc + region + '_tlinear.hdr.dat')

    # Delete the IRAF tables.
    os.remove(loc + 'tcreate.tab')
    os.remove(loc+ 'tlinear.tab')




# -----------------------------------------------------------------------------
# drrr
# -----------------------------------------------------------------------------

def create_R_plots(body, region, body_dict):
    """ yaxis: R 
        xaxis: Crater diameter, D
    

    geom_mean[of diam bin] = (sum(crater_diams))**(1/n_craters)

    R = ((geom_mean)**3) * n_craters) / A

    plot for this bin, log(geom_mean) vs log(R)

    References:
        pg 258-260 of text book

        Arvidson, R.E., et al. (1979) Standard techniques for presentation 
        and analysis of crater size-frequency data, Icarus 37: 467-474
    """
    
    # Have option to plot vs crater diam or vs the geom mean.
        # Get the bins of crater counts for crater diams.
    # Initialize lists.
    radii_all = []
    areas_all = []

    # Append together all lists of same region. 
    # Make lists of lists of same region.
    for key in body_dict.keys():
        if region in key:
            # Get radii from dict.
           radii = (body_dict[key])[2]
           radii_all.append(list(radii)) #+= list(radii)

           # Get the x and y lens in pixels
           x_len_pixels = (body_dict[key])[5]
           y_len_pixels = (body_dict[key])[6]

           # Transform to units of km
           x_len_km = transform_pixels_to_meters(body, x_len_pixels) / 1000.
           y_len_km = transform_pixels_to_meters(body, y_len_pixels) / 1000.

           # Get the area of the sub region in units of km2.
           print "Area in pixels, ", (x_len_pixels * y_len_pixels)
           print "Area in km, ", (x_len_km * y_len_km)
           area_km = (x_len_km * y_len_km)
           #crater_counts = crater_counts / (x_len_km * y_len_km)) #* ((x_len_pixels * y_len_pixels) / (x_len_km * y_len_km))
           areas = np.ones(len(radii)) * area_km
           areas_all.append(list(areas)) #+= list(areas)

    # Make sure that the list of radii is not empty.
    if radii_all != []:

        # Get the bins of crater counts for crater diams.
        crater_bins, crater_counts = bin_and_transform_crater_sizes(body, radii_all, areas_all, 'R') 

        print " "
        #print "Cumulative crater density of region ", region, ": ", crater_counts[1]

        fig=pylab.figure(figsize=(14.5,10.5))

        pylab.loglog(crater_bins, crater_counts, marker='d', linestyle='none') #size=100

        pylab.xlabel('LOG D^-n_craters [km]', fontsize=22, weight='bold')
        pylab.ylabel('LOG R', fontsize=22, weight='bold') 
        pylab.tick_params(axis='both', which='major', labelsize=20)
        pylab.title(body + ", " + region, fontsize=22)
        pylab.xlim([1, 100])
        if body == 'charon':
            pylab.ylim([.0001, 10])
        elif body == 'pluto':
            pylab.ylim([.00001, 10])
        pylab.savefig(path+'outputs/'+body+'_plots/R/'+region+'.png', \
                  bbox_inches='tight')   
        #pylab.show()
        pylab.close()

    else:
        print "No crater counts for region, " + region
        print "Skipping.... no R plot created."



def create_R_oplots(body, regions, body_dict):
    """ yaxis: R 
        xaxis: Crater diameter, D
    

    geom_mean[of diam bin] = (sum(crater_diams))**(1/n_craters)

    R = ((geom_mean)**3) * n_craters) / A

    plot for this bin, log(geom_mean) vs log(R)

    References:
        pg 258-260 of text book

        Arvidson, R.E., et al. (1979) Standard techniques for presentation 
        and analysis of crater size-frequency data, Icarus 37: 467-474
    """
    
    fig=pylab.figure(figsize=(14.5,10.5))
    pylab.xscale('log')
    pylab.yscale('log')

    for region in regions:

        # Initialize lists.
        radii_all = []
        areas_all = []

        # Append together all lists of same region. 
        # Make lists of lists of same region.
        for key in body_dict.keys():
            if region in key:
                # Get radii from dict.
               radii = (body_dict[key])[2]
               radii_all.append(list(radii)) #+= list(radii)

               # Get the x and y lens in pixels
               x_len_pixels = (body_dict[key])[5]
               y_len_pixels = (body_dict[key])[6]

               # Transform to units of km
               x_len_km = transform_pixels_to_meters(body, x_len_pixels) / 1000.
               y_len_km = transform_pixels_to_meters(body, y_len_pixels) / 1000.

               # Get the area of the sub region in units of km2.
               print "Area in pixels, ", (x_len_pixels * y_len_pixels)
               print "Area in km, ", (x_len_km * y_len_km)
               area_km = (x_len_km * y_len_km)
               #crater_counts = crater_counts / (x_len_km * y_len_km)) #* ((x_len_pixels * y_len_pixels) / (x_len_km * y_len_km))
               areas = np.ones(len(radii)) * area_km
               areas_all.append(list(areas)) #+= list(areas)

            print ">>>>>>>"
            print "radius_all:"
            print radii_all
            print "areas_all:"
            print areas_all

        # Make sure that the list of radii is not empty.
        if radii_all != []:
            # Get the bins of crater counts for crater diams.
            crater_bins, crater_counts = bin_and_transform_crater_sizes(body, radii_all, areas_all, 'R') 

            print " "
            print "R value of region ", region, ": ", crater_counts[1]

            print "NOW PLOTTING"
            print crater_bins
            print crater_counts

            color = return_color(body, region)

            pylab.scatter(crater_bins, crater_counts, marker='d', \
                c=color, edgecolors='none') #, linestyle='none', markeredgecolor='none'            
            pylab.plot(crater_bins, crater_counts, c=color)    

        else:
            print "No crater counts for region, " + region
            print "Skipping.... no R plot created."
    

    pylab.xlabel('LOG D^-n_craters [km]', fontsize=22, weight='bold')
    pylab.ylabel('LOG R', fontsize=22, weight='bold') 
    pylab.tick_params(axis='both', which='major', labelsize=20)
    pylab.title(body + ' R Plot', fontsize=22)
    pylab.xlim([1, 100])
    if body == 'charon':
        pylab.ylim([.0001, 10])
    elif body == 'pluto':
        pylab.ylim([.00001, 10])

    # Create the legend. 

    if body == 'charon':
        region_ls = ['south', 'belt', 'north', 'harad', 'mordor']
        for region in region_ls:
            pylab.scatter([], [], marker='d', c=return_color(body, region), label=region, edgecolors='none')
    elif body == 'pluto':
        region_ls = ['south', 'bbelt', 'tbelt', 'rbbelt', 'ltbelt', 'lheart', 'rheart', \
                   'lamp', 'snakeskin', 'rnorth','rnorthterm', 'north'] 
        for region in region_ls:
            pylab.scatter([], [], marker='d', c=return_color(body, region), label=region, edgecolors='none')
    pylab.legend(scatterpoints=1, loc='upper left') #, colors)

    pylab.savefig(path+'outputs/'+body+'_plots/R/all_regions.png', \
                  bbox_inches='tight')   
    pylab.close()




# -----------------------------------------------------------------------------
# non-normalized region plotting functions.
# -----------------------------------------------------------------------------

def create_diam_vs_cum_region_plots(body, region, radii):
    """Plots the crater count vs crater diameter over an entire region.
    """
    # Get the bins of crater counts for crater diams.
    
    crater_bins, crater_counts = bin_and_transform_crater_sizes(body, radii, 'cum')
    print crater_bins
    print crater_counts
    print "Cumulative crater density of region ", region, ": ", crater_bins[0]

    fig=pylab.figure(figsize=(14.5,10.5))
    pylab.scatter(np.log(crater_bins/1000.), np.log(crater_counts), marker='d') #size=100

    pylab.xlabel('LOG Crater diameter D [km]', fontsize=22, weight='bold')
    pylab.ylabel('LOG Cumalitive craters/region with diameter > D', fontsize=22, weight='bold') 
    pylab.tick_params(axis='both', which='major', labelsize=20)
    pylab.title(body + ', ' + region, fontsize=22)
    pylab.xlim([0,5])
    pylab.ylim([-.1,6])
    #axarr[0, 0].set_ylim([1.90,2.45])
    pylab.savefig(path+'outputs/'+body+'_plots/diam_vs_cum_region/'+region+'.png', \
                  bbox_inches='tight')   
    pylab.close()



def create_diam_vs_cum_region_oplots(body, body_dict):
    """OVER-plots the crater count vs crater diameter for all regions.
    
    To compare to Greenstreet et al. predictions for the three models...
    For Charon, see Table 5.
    For Pluto, see Table 4
    """
    
    fig=pylab.figure(figsize=(14.5,10.5))

    for region in body_dict.keys():

        crater_bins, crater_counts = bin_and_transform_crater_sizes(body, (body_dict[region])[2])
        

        color = return_color(body, region)

        pylab.scatter(crater_bins/1000., crater_counts, marker='d', \
            label=region, c=color, edgecolors='none') #size=100
    
    pylab.xlabel('LOG Crater diameter D [km]', fontsize=22, weight='bold')
    pylab.ylabel('LOG Cumalitive craters/region with diameter > D', fontsize=22, weight='bold') 
    pylab.tick_params(axis='both', which='major', labelsize=20)
    pylab.title(body + ' TOTAL (non-normalized) crater counts', fontsize=22)
    #pylab.xlim([0,5])
    #pylab.ylim([-.1,6])
    pylab.legend(scatterpoints=1)
    pylab.savefig(path+'outputs/'+body+'_plots/diam_vs_cum_region/all_regions.png', \
                  bbox_inches='tight')   
    pylab.close()


# -----------------------------------------------------------------------------
# km2 plotting functions.
# -----------------------------------------------------------------------------


def create_cum_km2_plots(body, region, body_dict):
    """Plots the crater count vs crater diameter over a km2 cutout
    of region.

    Parameters:
        region : string
            Region you want all the cutouts for. 
    """
    # Initialize lists.
    radii_all = []
    areas_all = []

    # Append together all lists of same region. 
    # Make lists of lists of same region.
    for key in body_dict.keys():
        if region in key:
            # Get radii from dict.
           radii = (body_dict[key])[2]
           radii_all.append(list(radii)) #+= list(radii)

           # Get the x and y lens in pixels
           x_len_pixels = (body_dict[key])[5]
           y_len_pixels = (body_dict[key])[6]

           # Transform to units of km
           x_len_km = transform_pixels_to_meters(body, x_len_pixels) / 1000.
           y_len_km = transform_pixels_to_meters(body, y_len_pixels) / 1000.

           # Get the area of the sub region in units of km2.
           print "Area in pixels, ", (x_len_pixels * y_len_pixels)
           print "Area in km, ", (x_len_km * y_len_km)
           area_km = (x_len_km * y_len_km)
           #crater_counts = crater_counts / (x_len_km * y_len_km)) #* ((x_len_pixels * y_len_pixels) / (x_len_km * y_len_km))
           areas = np.ones(len(radii)) * area_km
           areas_all.append(list(areas)) #+= list(areas)

    # Make sure that the list of radii is not empty.
    if radii_all != []:

        # Get the bins of crater counts for crater diams.
        crater_bins, crater_counts = bin_and_transform_crater_sizes(body, radii_all, areas_all, 'cum') 

        print " "
        print "Cumulative crater density of region ", region, ": ", crater_counts[1]

        fig=pylab.figure(figsize=(14.5,10.5))

        pylab.loglog(crater_bins, crater_counts, marker='d', linestyle='none') #size=100

        pylab.xlabel('LOG Crater diameter D [km]', fontsize=22, weight='bold')
        pylab.ylabel('LOG Cumalitive craters/km2 with diameter > D', fontsize=22, weight='bold') 
        pylab.tick_params(axis='both', which='major', labelsize=20)
        pylab.title(body + ", " + region, fontsize=22)
        pylab.xlim([1, 1000])
        if body == 'charon':
            pylab.ylim([.00001, .01])
        elif body == 'pluto':
            pylab.ylim([.00001, .01])
        pylab.savefig(path+'outputs/'+body+'_plots/cum_km2/'+region+'.png', \
                  bbox_inches='tight')   
        #pylab.show()
        pylab.close()

    else:
        print "No crater counts for region, " + region
        print "Skipping.... not cumulative plot created."




def create_cum_km2_oplots(body, body_dict, regions):
    """OVER-plots the crater count vs crater diameter for all regions.
    """
    # Get the bins of crater counts for crater diams.
    
    fig=pylab.figure(figsize=(14.5,10.5))
    pylab.xscale('log')
    pylab.yscale('log')

    for region in regions:

        # Initialize lists.
        radii_all = []
        areas_all = []

        # Append together all lists of same region. 
        for key in body_dict.keys():
            if region in key:
                print region
                # Get radii from dict.
                radii = (body_dict[key])[2]
                radii_all.append(list(radii)) #+= list(radii)

                # Get the x and y lens in pixels
                x_len_pixels = (body_dict[key])[5]
                y_len_pixels = (body_dict[key])[6]

                # Transform to units of km
                x_len_km = transform_pixels_to_meters(body, x_len_pixels) / 1000.
                y_len_km = transform_pixels_to_meters(body, y_len_pixels) / 1000.

                # Get the area of the sub region in units of km2.
                print "Area in pixels, ", (x_len_pixels * y_len_pixels)
                print "Area in km, ", (x_len_km * y_len_km)
                area_km = (x_len_km * y_len_km)
 
                areas = np.ones(len(radii)) * area_km
                areas_all.append(list(areas)) #+= list(areas)

        print ">>>>>>>"
        print "radius_all:"
        print radii_all
        print "areas_all:"
        print areas_all
        # Make sure that the list of radii is not empty.
        if radii_all != []:
            # Get the bins of crater counts for crater diams.
            crater_bins, crater_counts = bin_and_transform_crater_sizes(body, radii_all, areas_all, 'cum')  
            print " "
            print "Cumulative crater density of region ", region, ": ", crater_counts[1]

            print "NOW PLOTTING"
            print crater_bins
            print crater_counts

            color = return_color(body, region)

            pylab.scatter(crater_bins, crater_counts, marker='d', \
                c=color, edgecolors='none') #, linestyle='none', markeredgecolor='none'
        else:
            print "No crater counts for region, " + region
            print "Skipping.... not cumulative plot created."
    

    pylab.xlabel('LOG Crater diameter D [km]', fontsize=22, weight='bold')
    pylab.ylabel('LOG Cumalitive craters/km2 with diameter > D', fontsize=22, weight='bold') 
    pylab.tick_params(axis='both', which='major', labelsize=20)
    pylab.title(body + ' Cumulative Plot', fontsize=22)
    pylab.xlim([1, 1000])
    if body == 'charon':
        pylab.ylim([.00001, .01])
    elif body == 'pluto':
        pylab.ylim([.00001, .01])

    # Create the legend. 

    if body == 'charon':
        region_ls = ['south', 'belt', 'north', 'harad', 'mordor']
        for region in region_ls:
            pylab.scatter([], [], marker='d', c=return_color(body, region), label=region, edgecolors='none')
    elif body == 'pluto':
        region_ls = ['south', 'bbelt', 'tbelt', 'rbbelt', 'ltbelt', 'lheart', 'rheart', \
                   'lamp', 'snakeskin', 'rnorth','rnorthterm', 'north'] 
        for region in region_ls:
            pylab.scatter([], [], marker='d', c=return_color(body, region), label=region, edgecolors='none')
    pylab.legend(scatterpoints=1) #, colors)


    pylab.savefig(path+'outputs/'+body+'_plots/cum_km2/all_regions.png', \
                  bbox_inches='tight')   
    pylab.close()


def return_region_files(body):
    """List the region files that I have selected to be further
    cut with km2 boxes. 
    """
    if body == 'charon':
        region_files = ['charon_belt_413-826_413-826', 'charon_belt_413-826_826-1239', 
                        'charon_belt_826-1239_413-826', 'charon_belt_826-1239_826-1239', 
                        'charon_harad_826-1239_826-1239', 'charon_harad_1239-1652_413-826',
                        'charon_harad_1239-1652_826-1239',
                        'charon_mordor_1239-1652_413-826', 'charon_mordor_1239-1652_826-1239',
                        'charon_north_826-1239_413-826', 'charon_north_826-1239_826-1239',
                        'charon_south_0-413_413-826', 'charon_south_413-826_826-1239']
        # If using the full frame, ...
        x_l_corners = [0, 0,
                       160, 0,
                       0, 0,
                       0, 
                       190, 0,
                       0, 200,
                       0, 0]
        x_r_corners = [413, 413,
                       413, 413,
                       115, 413,
                       413,
                       413, 215,
                       413, 413,
                       413, 413]
        y_b_corners = [0, 230,
                       0, 0,
                       320, 0,
                       0,
                       180, 170,
                       290, 200,
                       200, 0]
        y_t_corners = [413, 413,
                       320, 220,
                       413, 200,
                       110,
                       413, 360,
                       413, 413,
                       413, 200]

    elif body == 'pluto':
        region_files = ['pluto_bbelt_0-724_1448-2172', 'pluto_bbelt_724-1448_724-1448',
                        'pluto_bbelt_724-1448_1448-2172', 'pluto_bbelt_1448-2172_724-1448',
                        'pluto_lamp_1448-2172_0-724', 'pluto_lamp_1448-2172_724-1448',
                        'pluto_ltbelt_2896-3620_724-1448',
                        'pluto_north_2896-3620_724-1448', 'pluto_north_2896-3620_1448-2172',
                        'pluto_north_2896-3620_2172-2896',
                        'pluto_rbbelt_724-1448_2896-3620', 
                        'pluto_rheart_724-1448_2896-3620', 'pluto_rheart_1448-2172_2896-3620',
                        'pluto_rnorth_2172-2896_2896-3620', 'pluto_rnorth_2896-3620_2896-3620',
                        'pluto_rnorthterm_2896-3620_2172-2896', 
                        'pluto_south_0-724_1448-2172', 
                        'pluto_tbelt_1448-2172_724-1448', 'pluto_tbelt_2172-2896_724-1448',
                        'pluto_tbelt_2896-3620_724-1448', 'pluto_tbelt_2896-3620_1448-2172']
        x_l_corners = [0, 0,
                       0, 0,
                       320, 0,
                       0,
                       100, 550,
                       0,
                       350,
                       380, 400,
                       200, 0,
                       300, 
                       0,
                       0, 0,
                       350, 0]
        x_r_corners = [724, 724,
                       530, 724,
                       724, 110,
                       220,
                       350, 724,
                       400,
                       620,
                       724, 724,
                       724, 550,
                       724, 
                       724,
                       724, 724,
                       724, 570]
        y_b_corners = [530, 0, 
                       0, 0, 
                       0, 130,
                       0,
                       20, 0,
                       100, 
                       200,
                       350, 0, 
                       350, 0,
                       100,
                       250,
                       150, 0,
                       0, 0]
        y_t_corners = [724, 724, 
                       724, 200,
                       724, 330,
                       724,
                       450, 310,
                       724, 
                       400,
                       724, 724,
                       724, 724,
                       724,
                       550,
                       724, 724,
                       380, 724]


    return region_files, x_l_corners, x_r_corners, y_b_corners, y_t_corners


def cutout_km2_box(body, region_file, x_l, x_r, y_b, y_t):
    """Cuts out a 1000m x 1000m box from a region file.
    Returns all radii in that box.

    ACTUALLY km is too small. That reaches limit of my resolution.
    Perhaps just divide the counts by size of region file in km2??
    And then scale by multiplying by pixel area.

    x_l, y_b define the lower left corner.
    """


    # Read in the region file.
    # Trick this function by giving entire rootname of the region file
    xcoord_data_ls, ycoord_data_ls, radius_ls, xcoord_canvas_ls, \
        ycoord_canvas_ls = read_in_region_files(body, region_file)

    
    # Loop through coordinates. Retain only radii that have coords
    # in that box.
    xcoord_data_box_ls = []
    ycoord_data_box_ls = []
    radius_box_ls = []
    xcoord_canvas_box_ls = []
    ycoord_canvas_box_ls = []

    for i in range(len(xcoord_data_ls)-1):
        x = xcoord_data_ls[i]
        y = ycoord_data_ls[i]
        if (((x > x_l) and (x < x_r)) and ((y > y_b) and (y < y_t))):
            xcoord_data_box_ls.append(xcoord_data_ls[i])
            ycoord_data_box_ls.append(ycoord_data_ls[i])
            radius_box_ls.append(radius_ls[i])
            xcoord_canvas_box_ls.append(xcoord_canvas_ls[i])
            xcoord_canvas_box_ls.append(ycoord_canvas_ls[i])        

    return xcoord_data_box_ls, ycoord_data_box_ls, radius_box_ls, \
           xcoord_canvas_box_ls, ycoord_canvas_box_ls 



# -----------------------------------------------------------------------------
# Wrapping functions.
# -----------------------------------------------------------------------------

def loop_through_km2boxes(body):
    """
    """
     
    region_files, x_l_corners, x_r_corners, y_b_corners, y_t_corners = return_region_files(body)

    body_dict = {}

    for i in range(len(region_files)-1):
        xcoord_data_ls, ycoord_data_ls, radius_ls, xcoord_canvas_ls, \
            ycoord_canvas_ls = cutout_km2_box(body, region_files[i], \
                x_l_corners[i], x_r_corners[i], y_b_corners[i], y_t_corners[i])

        x_len_pixels = x_r_corners[i] - x_l_corners[i]
        y_len_pixels = y_t_corners[i] - y_b_corners[i]

        body_dict[region_files[i]] = [xcoord_data_ls, ycoord_data_ls, radius_ls, \
                                  xcoord_canvas_ls, ycoord_canvas_ls, x_len_pixels, y_len_pixels]
    
    return body_dict


def plot_km2boxes(body, body_dict):
    """
    """
    # The regions lists.
    if body == 'charon':
        regions = ['south', 'belt', 'north', 'harad', 'mordor'] 
    elif body == 'pluto':
        regions = ['south', 'bbelt', 'tbelt', 'rbbelt', 'ltbelt', 'lheart', 'rheart', \
                   'lamp', 'snakeskin', 'rnorth','rnorthterm', 'north'] 

    for region in regions:
        print region
        create_cum_km2_plots(body, region, body_dict)

    create_cum_km2_oplots(body, body_dict, regions)



def plot_regions(body, body_dict):
    """Loops through all regions creating the plots.
    """
    
    for region in body_dict.keys():
        print region
        create_diam_vs_cum_region_plots(body, region, (body_dict[region])[2])
    create_diam_vs_cum_region_oplots(body, body_dict)



def loop_through_regions(body):
    """Loops through all the regions
    of a given body.  Returns all of that body's stats.
    Stores the regions in a dict.
    """

    # The regions lists.
    if body == 'charon':
        regions = ['south', 'belt', 'north', 'harad', 'mordor'] 
    elif body == 'pluto':
        regions = ['south', 'bbelt', 'tbelt', 'rbbelt', 'ltbelt', 'lheart', 'rheart', \
                   'lamp', 'snakeskin', 'rnorth','rnorthterm', 'north'] 

    # Initialize count of ALL craters.
    crater_tot = 0

    # The dictionary is for entire body.
    # The keys will be regions.
    # The values for each key will be the six lists returned from file.
    body_dict = {}
    
    for region in regions:

        xcoord_data_ls, ycoord_data_ls, radius_ls, xcoord_canvas_ls, \
            ycoord_canvas_ls = read_in_region_files(body, region)

        body_dict[region] = [xcoord_data_ls, ycoord_data_ls, radius_ls, \
                             xcoord_canvas_ls, ycoord_canvas_ls]

        region_tot = len(radius_ls)
        print ">>>><<<<"
        print body, " ", region, " number of craters: ", region_tot
        print ">>>><<<<"

        crater_tot += region_tot
    
    print ">>>>>><<<<<<"
    print body, " TOTAL number of craters: ", crater_tot
    print ">>>>>><<<<<<"

    return body_dict


def plot_R(body, body_dict):
    """
    """
    # The regions lists.
    if body == 'charon':
        regions = ['south', 'belt', 'north', 'harad', 'mordor'] 
    elif body == 'pluto':
        regions = ['south', 'bbelt', 'tbelt', 'rbbelt', 'ltbelt', 'lheart', 'rheart', \
                   'lamp', 'snakeskin', 'rnorth','rnorthterm', 'north'] 

    for region in regions:
        print region
        create_R_plots(body, region, body_dict)

    create_R_oplots(body, regions, body_dict)



# -----------------------------------------------------------------------------
# Main functions.
# -----------------------------------------------------------------------------

def print_basins():
    print "Charon basin?"
    print transform_pixels_to_meters('charon', 130.002852871) * 2
    print "Pluto basin?"
    print transform_pixels_to_meters('pluto', 266.867276935) * 2


def R_main():
    """
    """
    charon_dict = loop_through_km2boxes('charon')
    plot_R('charon', charon_dict)

    pluto_dict = loop_through_km2boxes('pluto')
    plot_R('pluto', pluto_dict)

def km2_main():
    """
    """
    charon_dict = loop_through_km2boxes('charon')
    plot_km2boxes('charon', charon_dict)

    pluto_dict = loop_through_km2boxes('pluto')
    plot_km2boxes('pluto', pluto_dict)
    
def region_main():
    """
    """
    # Get all the data.
    charon_dict = loop_through_regions('charon')
    pluto_dict = loop_through_regions('pluto')

    # Plot the data.
    plot_regions('charon', charon_dict)
    plot_regions('pluto', pluto_dict)

def charon_vs_pluto_main():
    """
    """


if __name__=='__main__':

    #region_main()
    #km2_main()
    R_main()
    #print_basins()
    #print "Pluto box length:", transform_pixels_to_meters('pluto', 724)
    #print "Charon error in km (2 pix): ", transform_pixels_to_meters('charon', 2)
    #print "Pluto error in km (2 pix): ", transform_pixels_to_meters('pluto', 2)


