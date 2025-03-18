# %%
from osgeo import gdal
import numpy as np

"""
Step 1:
    rasterize waterbodies
        parameters: waterbodies as vector format
    rasterize cutting area
        parameters: cutting area
Step 2:
    Calculate distance to waterbodies
        parameters: waterbodies
Step 3:
    Define initial raster for the riparian area
        parameters: minimum width,cutting area,distance2waterbody
Step 4:
    loop to dtw values
        for i in initial_raster:
            dtw_max = 
Step 5:
    vektorize the riparian zone
        attributes = mean(dtw),max(dtw)

"""
