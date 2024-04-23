Notes provided by Johannes Mohrmann and added to file by Ehsan Erfani on Feb. 2022:


The code for extracting the GOES data along trajectories was part of the original CSET project. The code is still accessible at UW, but it was never packaged to the UW trajectory since 2015. The most recent version of the code is here: /home/disk/p/jkcm/Code/Lagrangian_CSET/GOES_extractor.py.

1) for this project, the gridded pixel-level data is specified with the goes_folder variable, (goes_folder = '/home/disk/eos4/jkcm/Data/CSET/GOES/VISST_pixel'). The data here is the pixel-level data, not the lat/lon gridded data that you posted. We initially got these pixel images as part of the CSET project, but the data is downloadable here: https://satcorps.larc.nasa.gov/prod/goes-west/visst-pixel-netcdf. I will warn that the site is sometimes very slow to respond.

2) the code block for the project is in part in that Goes_extractor.py script, and partly in the imported hsproject/plots.py and hsproject/util.py scripts. Note that the scripts here will extract variables as well as make plots of them. 

3) the way the code worked is that it reads every trajectory file (generated using HYSPLIT) in a given list (here everything matching '/home/disk/eos4/jkcm/Data/CSET/Trajectories/*/analysis.UW_HYSPLIT_GFS.*.airmass_trajectories_500m_-48.txt'), and for each file extracts GOES data/makes plots.
