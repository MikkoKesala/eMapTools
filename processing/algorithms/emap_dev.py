# -*- coding: utf-8 -*-
"""
Develop tools for qgis

@author: mkesala
"""

#%%
import numpy as np
from osgeo import gdal,gdal_array
import geopandas as gpd
from scipy.ndimage import gaussian_filter
import pandas as pd
from qgis.core import QgsVectorLayer,QgsField,QgsFields,QgsFeature,QgsGeometry
from geotools2 import simple_treesegmentation,zonal_stastics,treeTops,calculateNDVI,singleTreeMapping,gpd2qgis
import matplotlib.pyplot as plt
import rasterio
import tempfile
from rasterio.plot import show
from rasterio.enums import Resampling
from shapely.geometry import Point
from PyQt5.QtCore import QVariant
# %%

from qgis.core import (
    QgsRasterLayer,
    QgsProject,
    QgsApplication,
    QgsProcessingFeatureSourceDefinition,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterNumber,
)
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
#import processing

# %%

# Initialize QGIS application (for standalone script)
QgsApplication.setPrefixPath("C:/Program Files/QGIS 3.28.4/apps/qgis-ltr/bin", True)
qgs = QgsApplication([], False)
qgs.initQgis()


# %%
def localConnectivity(input:QgsVectorLayer,field:str):
    print ("calculating local connecitity for the "+str(input.name())+" at field "+str(field))

# %%
def globalMoran(input:QgsVectorLayer):
    print ("calculating global Moran I")
# %%


# %%



# %%


# %%


# %%
testCHM = "testfiles/testCHM.tif"
testResample = "testfiles/watershed_resample.tif"
treeSegments = "testfiles/watershed_result2.tif"
testOrto = "testfiles/testOrto.tif"
testNDVI = "testfiles/testNDVI23.tif"

#testND = calculateNDVI(testOrto,None)

#print (testND)

trees = singleTreeMapping(testCHM,testOrto,2)

print (trees)
#simple_treesegmentation(testCHM,2,treesegments=treeSegments)

"""
with rasterio.open(treeSegments) as src:
    ts = src.read(1)
    transform = src.transform
    crs = src.crs
    print (src.width,src.height)
with rasterio.open(testOrto) as src:
    orto = src.read(1)
    transform = src.transform
    crs = src.crs
    print (src.width,src.height)"""

#resample_segments = resample_raster(treeSegments,testResample,0.5) 
#print (np.max(ts))
#zonals = zonal_stastics(testOrto,treeSegments,1)

#zonals.sort_values("count")
# %%
# Save the watershed result to a new raster file
output_path = 'testfiles/watershed_result.tif'

# %%
# %%
# Plot the original CHM
fig, ax = plt.subplots(figsize=(20, 16))
chm = rasterio.open(testCHM)
show(chm,ax=ax,cmap='viridis')
ax.set_title('CHM with Local Maxima')

# Overlay local maxima coordinates
#ax.scatter(coordinates[:, 0], coordinates[:, 1], color='red', s=1, marker='x', label='Local Maxima')

# Add legend and show plot
ax.legend()
plt.show()
# %%
testpack = "testfiles/testPack.gpkg"
testlayer = "testtreemap"
testfield = "chm_height"
test_qgis = QgsVectorLayer(testpack+"|layername="+testlayer,testlayer,"ogr")
data_df = gpd.read_file(testpack,layer = testlayer)

testtest = gpd2qgis(data_df)
c = testtest.featureCount()
field_names = [field.name() for field in testtest.fields()]
print (field_names,c)

# %%
localConnectivity(test_qgis,testfield)
# %%
gaussian = lambda x, y,sigma,size: (1 / (2 * np.pi * sigma**2)) * np.exp(-((x - (size - 1) / 2)**2 + (y - (size - 1) / 2)**2) / (2 * sigma**2))
print (gaussian(0,0,1,2))


# %%
test = np.array([[1,2,3],
                 [4,5,6],
                 [7,8,9]])

test2 = np.extract(test>4,test)
print(np.mean(test2))
# %%
