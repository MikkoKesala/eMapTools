# -*- coding: utf-8 -*-
"""
Ecological indices of emaptools

@author: mkesala
"""

import geopandas as gpd
import numpy as np
from .geotools2 import valuesByDistance
from .geostats import getLstats
from collections import Counter

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
    

    new_gdf = new_gdf.drop(['lstat','lstat2','list_'+height,'list_'+ndvi,'distances'],axis=1)
    #new_gdf = new_gdf.drop(['lstat','lstat2'],axis=1)
        
    return new_gdf

def diversityEquations(specieslist,method):
    
    # Calculate species count using Counter
    species_counts = Counter(specieslist)

    # Total number of individuals (N)
    N = sum(species_counts.values())
    L = len(species_counts.values())
    if N>0 and L>1:
        proportions = [i/N for i in species_counts.values()]
        
        methods = {'simpsonOld':np.round(1 - np.sum(np.fromiter([(count * (count - 1)) / (N * (N - 1)) for count in species_counts.values()],float)),3),
                    'simpson':np.round(1-np.sum(np.fromiter([p**2 for p in proportions],float)),3),
                    'simpson_evenness':np.round((1-np.sum(np.fromiter([p**2 for p in proportions],float)))/(1-1/L),3),
                    'shannon':np.round(-np.sum(np.fromiter([p * np.log(p) for p in proportions],float)),3),
                    'pielou':np.round(-np.sum(np.fromiter([p * np.log(p) for p in proportions],float) / np.log(N)),3)}
        value = methods[method]
    else:
        value = 0.0
    
    return value

def diversityIndices(gpd,treespecies,distance,method):

    new_df = gpd
        
    # local factors
    new_df = valuesByDistance(gpd,treespecies,distance)
    if type(method) == list:
        for m in method:
            new_df[m] = new_df['list_'+treespecies].apply(lambda x: diversityEquations(x,m))
    else:
        new_df[method] = new_df['list_'+treespecies].apply(lambda x: diversityEquations(x,method))
    
    new_df = new_df.drop(['distances'],axis=1)
    
    return new_df