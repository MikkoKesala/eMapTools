# -*- coding: utf-8 -*-
"""
geostatistic algorithms of emaptools

Created: 20.08.2024

@author: mkesala
"""

import numpy as np
import geopandas as gpd
import pandas as pd
from .geotools2 import valuesByDistance
from collections import Counter

def getLstats(values):
    """
    L-stastistics calculation based on ...

    Parameters
        values: list of values
    Outputs
        lcv: Gini coeffiecient
        lskew: Skewness of L-moments
    """
    n = len(values)
    sample = np.array(values)
    sample = np.sort(sample.reshape(n))[::-1]
    b0 = np.mean(sample)
    b1 = np.array([(n - j - 1) * sample[j] / n / (n - 1)
           for j in range(n)]).sum()
    b2 = np.array([(n - j - 1) * (n - j - 2) * sample[j] / n / (n - 1) / (n - 2)
           for j in range(n - 1)]).sum()
    lmom1 = b0
    lmom2 = 2 * b1 - b0
    lmom3 = 6 * (b2 - b1) + b0
    lcv = np.round(lmom2 / lmom1,3)
    lskew = np.round(lmom3 / lmom2,3)
    return lcv,lskew

def getisord(gpd,value,distance,relationship):
    """
    Getis-ord (Gi*) algorithm based on Getis & Ord (1992) article.

    Parameters
        gpd: geopandas dataframe
        value: fieldname of value to the Getis-Ord z-value calculation (high values means hotspot area of the value and vice versa)
        distance: distance threshold (linear,quadric,gaussian,fixed)
    Outputs:
        new_df: geopandas dataframe with GiZ (Getis-Ord Z-score value) and nneigh (number of neighbors within the distance)
    """
    
    #global factors
    new_df = gpd
    x = np.mean(new_df[value])
    n = len(new_df)
    s = np.std(new_df[value])

    #local factors
    new_df = valuesByDistance(gpd,value,distance)
    new_df['wlist'] = new_df['distances'].apply(lambda x:spatialWeights(x,relationship,distance))
    new_df['w_sum'] = new_df['wlist'].apply(sum)
    new_df['wx_sum'] = new_df.apply(lambda new_df: np.sum([a*b for a,b in zip(new_df['wlist'],new_df['list_'+value])]),axis=1)   
    new_df['nneigh'] = new_df['distances'].apply(lambda x:len([1 for i in x]))
    
    # spatial indicator
    new_df['slocal'] = np.sqrt((new_df['w_sum'] * (n-new_df['w_sum']))/(n-1))
    new_df[value+'_GiZ'] = (new_df['wx_sum'] - (x * new_df['w_sum'])) / (s * new_df['slocal'])
    
    # cleaning the dataframe
    new_df = new_df.drop(['wlist','wx_sum','w_sum','slocal','distances','list_'+value],axis=1)
   
    return new_df
    

def spatialWeights(distances,method,threshold):

    methods = {'linear': [(threshold - x) / threshold for x in distances],
                'quadric':[((threshold - x) / threshold)**2 for x in distances],
                'squared':[np.sqrt((threshold - x) / threshold) for x in distances],
                'gaussian':[np.exp(-(x**2/(2*threshold**2))) for x in distances],
                'fixed':[1 for x in distances]}
    
    return methods[method]

def diversityIndicator(gdf,value,radius,method,abudance_round=0):
    
    #local values
    new_gdf = gdf
    new_gdf['round_'+value] = np.round(gdf[value],abudance_round)
    new_gdf = valuesByDistance(new_gdf,'round_'+value,radius)

    if method == 'raoq':
        new_gdf['localRaoQ'+str(radius)+value] = new_gdf.apply(lambda x: raoQindex(x['round_'+value],x['list_round_'+value],x['distances']),axis=1)

        #global value
        analy_df = new_gdf
        new_gdf['RaoQ'+str(radius)+value] = new_gdf.apply(lambda x:np.sum(analy_df.loc[[i for i in x['nearindices']]]['localRaoQ'+str(radius)+value]),axis=1)
    
    else:
        new_gdf['gini'+str(radius)+value] = new_gdf['list_round_'+value].apply(lambda x:giniCoefficient(x))
    
    new_gdf = new_gdf.drop(['distances','nearindices','round_'+value,'list_round_'+value],axis=1)

    gdf = None
    analy_df = None

    return new_gdf

def raoQindex(value,values,distances):
    count_values = Counter(values)
    N = sum(count_values.values())
    L = len(count_values.values())

    pi = count_values[value] / N
    if L>1:
        raoq = np.sum([pi*(count_values[pj]/N)*distances[c] for c,pj in enumerate(values)])
    else:
        raoq = 0

    return raoq

def giniCoefficient(data):
    """
    Calculate the Gini coefficient for a list of values.
    
    Parameters:
    data (list): List of values (e.g., tree sizes).
    
    Returns:
    float: Gini coefficient.
    """
    
    data = np.array(data)
    if np.amin(data) < 0:
        data -= np.amin(data)  # Values cannot be negative
    data = data.astype(float)  # Ensure data is float to avoid casting issues
    data += 0.0000001  # Values cannot be 0
    data = np.sort(data)  # Sort the data
    index = np.arange(1, len(data) + 1)  # Index array
    n = len(data)
    
    gini = (np.sum((2 * index - n - 1) * data)) / (n * np.sum(data))
    return gini