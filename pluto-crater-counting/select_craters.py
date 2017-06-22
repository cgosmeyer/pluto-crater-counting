#! /usr/bin/env python

"""Interactively records a crater diameter given a center pixel the
curser hovers over. Uses matplotlib event handling. 

Author:

    C.M. Gosmeyer, Oct. 2015
    
Use:
    
    >>> python select_craters.py
    
References:

matplotlib.org/users/event_handling.html
stackoverflow.com/questions/33569626/matplotlib-responding-to-click-events
stackoverflow.com/questions/13714454/specifying-and-saving-a-figure-with-exact-size-in-pixels


Still To Do:
    * Check math on radius scaling ... need make sure I can
      report the number of pixels consistently across images.
      Right now, radius, x, and y are all in units of pixels.
      When writing up report, need convert pixels to meters.
      (So what is meters/pixel?)
    * Need record an 'error weighting' of the craters based on the 
      latitude and longitude. 
      Say center of body is 0,0. So increase the uncertainty with distance
      from center.

      Now that I'm set to not use a transformed image, I need come up
      with a reasonable error estimate based on these somewhat related
      terms:
      1. Surface illumination. 
      2. Proximity to terminator (might be same as 1.)
      3. Proximity to limb.

      Item 3. will be the largest one.  It needs depend on lat/long.
      So it will increase with distance from 0,0 at body center.
      I should be able to calculate it beased on xdata,ydata + 
      x,y displacement from 0,0 at corner.

      I can write a function to handle the transoformation from
      0,0 at corner to 0,0 at center.

      So there will be a weight for both 
      1. radius (from limb proximity)
      2. count number (from illumination and pixel resolution)

    
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
#from mpldatacursor import datacursor
import os
from skimage import data, color
from skimage.io import imread

# Because matplotlib's backend is vector-based, it won't size figures in pixels.
# drrr
# See stackoverflow.com/questions/13714454/specifying-and-saving-a-figure-with-exact-size-in-pixels
monitor_dpi = 192
path = '/path/to/main/dir/containing/outputs/data/and/scripts/dirs/'


# -----------------------------------------------------------------------------
# Region-displaying functions
# -----------------------------------------------------------------------------

def return_region_dict(body, width, height):
    """Returns a dictionary where the geological regions are
    matched to image subsections.

    Parameters
    ----------
    body : string
        Either 'pluto' or 'charon'.
    width : float
        Width of the body's image, in pixels.
    height : float
        Height of the body's image, in pixels. 

    Returns
    -------
    region_dict : dictionary
        Keys of region id, values of list of coords of upper
        and lower region corners, in pixels.
    """

    if body == 'pluto':
        divisor = 6.
    elif body == 'charon':
        divisor = 4.

    # In for loop, fill out the x,y coords of corners.
    x_subwidth = int(width / divisor)
    y_subwidth = int(width / divisor)

    # Initialize key and dictionary
    region_dict = {}

    keyname = 0
    for xlower in range(0, x_subwidth*int(divisor), x_subwidth):
        for ylower in range(0, y_subwidth*int(divisor), y_subwidth):
            xupper = xlower + x_subwidth
            yupper = ylower + y_subwidth
 
            # Fill out the dictionary.
            region_dict[keyname] = [xlower, ylower, xupper, yupper]
            print keyname, xlower, ylower, xupper, yupper
            keyname += 1

    return region_dict


def display_histogram(imdata):
    """ Displays histogram of the region's brightness, to aid
    in adjusting the contrast in the final region image to be
    used for counting craters.

    Parameters
    ----------
    imdata : image array
        Array of the region.
    """

    plt.hist(imdata.ravel(), bins=256,  fc='k', ec='k')


def display_image(image, image_id, region=''):
    """ Displays the region image on which will interactively select and
    measure craters.

    Parameters
    ----------
    image : string
        JPG image of the region.
    image_id : int
        ID of the region of full image (generated in 
        :func:'return_region_dict')
    region : string
        My name for the region, to be used in categorizing craters
        by surface type. 
    """
    
    # Find the body depicted in image.
    global body
    body = (image.split('/')[len(image.split('/'))-1]).split('.jpg')[0]
    print body

    # Read in image.
    # Need flip first cuz python
    # http://stackoverflow.com/questions/14589642/python-matplotlib-inverted-image
    imdata = np.flipud(imread(image))
    #print imdata

    # Display the histogram of pixel values. Use this to clip the subimage
    # so get best contrast. 
    display_histogram(imdata)

    # Get the width and height.
    width, height = imdata.shape
    print "Width, Height:", width, height
    
    # Use the dictionary to get desired subimage.
    region_dict = return_region_dict(body, width, height)
    xlower = (region_dict[image_id])[0]
    ylower = (region_dict[image_id])[1]
    xupper = (region_dict[image_id])[2]
    yupper = (region_dict[image_id])[3]
    imdata_sub = imdata[xlower:xupper, ylower:yupper]
    print xupper, xlower, yupper, ylower
    print imdata_sub.shape
    
    # Create filename for text file that will hold crater counts.
    global filename
    filename = body + '_' + region + '_' + str(xlower)+'-'+str(xupper)+'_'+str(ylower)+'-'+str(yupper) + '.txt'
    print filename

    # Scale the image so can fit into matplotlib window (figsize must be inches and cannot 
    # take units of pixels, drrrr). Need keep track of 'scaler' here -- need use it to
    # normalize all crater radii I measure. 
    global scaler
    scaler = monitor_dpi / 7. 

    # Initialize the figure and plot.
    global fig
    global axes
    
    # Resolution might actually be better if don't specify fig size, but
    # rather resize by hand.
    fig, axes = plt.subplots() #figsize=(width_scaled, height_scaled))

    # Display image in greyscale.
    # clim estimated from histogram
    axes.imshow(imdata_sub, cmap=cm.Greys_r, origin='lower', clim=(10, 235))  

    # Begin interactivity. 
    # Click on center of crater and select a radius.
    # A circle will be drawn at that location and the x,y and radius
    # saved to a file.
    # A marked crater may be deleted by selecting it and typing 'd'.

    # Connection ID 1. Select crater radius via a key press event.
    cid_1 = fig.canvas.mpl_connect('key_press_event', onpress)
    # Set off picking event, which can be used to delete picked circle object.
    fig.canvas.mpl_connect('pick_event', onpick)
    
    # Connection ID 2. Define own crater radius via mouse click event.
    cid_2 = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()


def save_crater_records(xcoord_data, ycoord_data, radius, \
                        xcoord_canvas, ycoord_canvas, scaler):
    """Opens a file, or creates if it doesn't already exist, and 
    appends the x and y coords, and radius found from selection.

    Parameters
    ----------
    xcoord_data : string/float
        The x-coordinate of crater center, in pixels.
    ycoord_data : string/float
        The y-coordinate of crater center, in pixels
    radius : string/float
        The radius of the crater, in pixels.
    xcoord_canvas : string/float
        The x-coordinate of crater center on matplotlib canvas.
    ycoord_canvas : string/float
        The y-coordinate of crater center on matplotlib canvas.
    scaler : string/float
        The scaling factor to reduce size so can fit into matplotlib
        window.
    """
    
    with open(path+'outputs/' + body + '/'+filename, "a+") as f:
        f.write(str(xcoord_data) + "\t")
        f.write(str(ycoord_data) + "\t")
        f.write(str(radius) + "\t")
        f.write(str(xcoord_canvas) + "\t")
        f.write(str(ycoord_canvas) + "\t")
        f.write(str(scaler))
        f.write("\n")


def delete_crater_records(xcoord_data, ycoord_data, radius):
    """Open file and remove the line containing the given
    data x and y coords and radius (it will be assumed the 
    combination of these values is unique!)

    Parameters
    ----------
    xcoord_data : string/float
        The x-coordinate of crater center, in pixels.
    ycoord_data : string/float
        The y-coordinate of crater center, in pixels
    radius : string/float
        The radius of the crater, in pixels.
    
    """

    # Open file and read lines.
    file_read = open(path+'outputs/'+body+'/'+filename, 'r')
    lines = file_read.readlines()
    print lines # for debugging
    file_read.close()     
        
    # Now open file again, but for writing.
    # Then read back in lines, it will be overwritten.
    file_write = open(path+'outputs/'+body+'/'+filename, 'w')
    # Find line where match radius and x,y coords.
    # Remove the matched line and write the truncated lines
    # back into file.
    for line in lines:
        print '----'
        print str(xcoord_data), str(ycoord_data)
        if ((str(xcoord_data) not in line) and \
            (str(ycoord_data) not in line)):

            file_write.write(line)

    file_write.close()


def approximate_error(xcoord_data, ycoord_data, xlower, ylower, body_radius):
    """Maybe emulate GRE: the weight will be a percent.
    For example, a crater I find right at the limb edge (89deg) will
    have a weight of ((90-89)/(radius-dist_from_center)*(radius/90) 

    In final count, this point will be multiplied by the weight.
    (So count as '...%' of a crater in statistics of all craters 
        at that radius or whatever)

    A crater right at the center 0,0 will have a weight of 1.
    (So 100 percent in final crater count.)
    """

    # never had time for this! in end, didn't make much difference, since
    # uncertainty in crater-counting is overweighted by its subjectivity 
    # anyway... better error measurement if can compare crater-counts 
    # from multiple counters.

    return error_weight


def return_crater_dict():
    """Returns dict of pre-defined crater radii (pixels).
    """
    crater_dict = {'1':0.1, 
                   '2':0.2, 
                   '3':0.3,
                   '4':0.4,
                   '5':0.5,
                   '6':0.6,
                   '7':0.7,
                   '8':0.8,
                   '9':0.9,
                   '0':1.0} 

    return crater_dict



# -----------------------------------------------------------------------------
# The event-handling functions
# -----------------------------------------------------------------------------


def define_radius(event):
    """ Calculates the radius by measuring the distance
    between user-defined lines.
    """
    xy = plt.ginput(1)
    startx = event.xdata
    starty = event.ydata
    x = [startx, xy[0][0]]
    y = [starty, xy[0][1]]
    print "x of line", x
    print "y of line", y

    # The radius is found from the distance formula.
    radius = np.sqrt( (x[1] - x[0])**2  + (y[1] - y[0])**2 )

    # Finally draw and record the crater.
    draw_and_record_crater(event, radius)


def draw_and_record_crater(event, radius):
    """
    """
    xcoord_data = event.xdata
    ycoord_data = event.ydata
    xcoord_canvas = event.x
    ycoord_canvas = event.y

    print xcoord_data, ycoord_data, radius, xcoord_canvas, ycoord_canvas, scaler

    # Save the crater circle in file. Also record the scaling of the figure (what
    # integer the height and width of figure were divided by).
    save_crater_records(xcoord_data, ycoord_data, radius, xcoord_canvas, ycoord_canvas, scaler)

    circ = plt.Circle((event.xdata, event.ydata), radius=radius, color='r', fill=False, picker=True)
    axes.add_patch(circ)
    axes.figure.canvas.draw()


def delete_crater(event):
    """
    Removes selected object from canvas.
    """
    ax = plt.gca()
    if axes.picked_object:
        print "deleting object", axes.picked_object
        axes.picked_object.remove()
        # Delete from file.
        delete_crater_records(axes.picked_object.center[0], \
                              axes.picked_object.center[1], \
                              axes.picked_object.radius)
        # Delete from canvas.
        axes.picked_object = None
        axes.figure.canvas.draw()




def onclick(event):
    """
    Define own radius by double-clicking center and single-clicking rim.
    A circle will then be drawn from that radius.
    """
    if event.dblclick:
        if event.button == 1:
            define_radius(event)



def onpress(event):
    """Logic for deciding how to handle a command following
    a press of mouse.  
    """
    # Get dictionary of pre-defined crater radii.
    crater_dict = return_crater_dict()
    print event.key
    # Delete a circle from canvas
    if event.key == 'd':
        delete_crater(event)
    # Save current canvas.
    elif event.key == 'q':
        plt.savefig(path+"outputs/" + body + "/"+filename.split('txt')[0]+'.jpg', dpi=monitor_dpi, bbox_inches='tight')
    # Use a pre-defined radius (if crater especially small)
    elif event.key in crater_dict.keys():
        radius = float(event.key)
        draw_and_record_crater(event, radius)
    else:
        pass


def onpick(event):
    """
    Handles the pick event - if an object has been picked, store a 
    reference to it. We do this by adding a reference named 'stored_pick'
    to the axes object. Note that in python we can dynamically add an 
    attribute variable (stored_pick) to an existing object - even one 
    that is produced by a library as in this case.
    """
    this_artist = event.artist #the picked object is avable as event.artist
    print this_artist
    plt.gca().picked_object = this_artist


# -----------------------------------------------------------------------------
# The Main.
# -----------------------------------------------------------------------------

if __name__=='__main__':
    # User needs intialize 
    # 1. The image ['pluto.jpg', 'charon.jpg']
    image = os.path.join(path, 'data/pluto.jpg')

    # 2. The ID on image grid [0-35] for Pluto, [0-15] for Charon 
    image_id = 8
 
    # 3. The region ['south', 'bbelt', 'tbelt', 'rbbelt', 'ltbelt' lheart', 'rheart', 
    #                'lamp', 'snakeskin', 'rnorth','rnorthterm', 'north'] for Pluto
    #               ['south', 'belt', 'north', 'harad', mordor'] for Charon
    # Attempting to list regions from the bottom corner left-right, bottom-up
    region = 'tbelt'

    display_image(image, image_id, region)


