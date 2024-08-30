# -*- coding: utf-8 -*-
"""
geostatistic algorithms of emaptools

Created: 20.08.2024

@author: mkesala
"""

import numpy as np
import geopandas as gpd
import pandas as pd


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