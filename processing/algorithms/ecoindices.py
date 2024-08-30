# -*- coding: utf-8 -*-
"""
Ecological indices of emaptools

@author: mkesala
"""

import geopandas as gpd
import numpy as np
from geotools2 import valuesByDistance
from geostats import getLstats

def mixedPatchIndex(gdf,height,ndvi,radius):
    gdf[ndvi] = gdf[ndvi]+1
    new_gdf = valuesByDistance(gdf,height,radius)
    new_gdf = valuesByDistance(gdf,ndvi,radius)

    new_gdf['lstat']= new_gdf['list_'+height].apply(lambda x: getLstats(x))
    new_gdf[height+'_lcv'] = new_gdf['lstat'].apply(lambda x: x[0])
    new_gdf[height+'_lskew'] = new_gdf['lstat'].apply(lambda x: x[1])
    
    new_gdf['lstat2']= new_gdf['list_'+ndvi].apply(lambda x: getLstats(x))
    new_gdf[ndvi+'_lcv'] = new_gdf['lstat2'].apply(lambda x: x[0])
    new_gdf[ndvi+'_lskew'] = new_gdf['lstat2'].apply(lambda x: x[1])

    new_gdf['mpi'] = (1 - (abs(new_gdf[height+'_lcv']-0.33))) + (1 - (abs(new_gdf[ndvi+'_lcv']-0.33))) + (1-abs(new_gdf[height+'_lskew'])) + (1-abs(new_gdf[ndvi+'_lskew']))
    new_gdf['mpi'] = np.round(new_gdf['mpi'] / 4,3)  
    

    #new_gdf = new_gdf.drop(['lstat','lstat2','list_'+height,'list_'+ndvi,'distances'],axis=1)
    new_gdf = new_gdf.drop(['lstat','lstat2'],axis=1)
        
    return new_gdf

