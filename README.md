# pluto-crater-counting
Python scripts to aid user in selecting craters on Charon and Pluto, measure the craters' diameteters and locations, and create crater density plots in order estimate surface ages.

Written as a project for a 2015 Johns Hopkins University course "Planetary Surface Processes," taught by Kevin Lewis. 

Before running any scripts, change the 'path' variable at the top to the path containing this repo, containing 'data', 'outputs', and 'pluto-crater-counting' directories. 

1. Transform full Pluto and Charon images to greyscale and place in 'data'.

2. In main of `select_craters.py`, select the image of body, the id of the image region, and the name of the region to run the crater selector on. When run the script, will be given pop-up matplotlib windows on which you can select the center of a suspected crater and measure the radius. The values of the selected crater coordinates and radius will be saved in 'outputs', as well as the image of the region showing the selected craters. 

3. In main of `make_crater_plots.py`, select the plots/statistic function you wish to run. Crater density plots will appear in 'outputs'.

Can ignore `mapflat_test.pro` and `mapflat.pro`, failed attempts to project Pluto and Charon as flat maps, with pixels equating to equal physical areas. 