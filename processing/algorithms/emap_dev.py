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
from geotools2 import gaussian_kernel,convolve_2d,calcFocal,gpd2qgis,valuesByDistance,singleTreeMapping,resample_raster
from geostats import getLstats
import matplotlib.pyplot as plt
import rasterio
import tempfile
from rasterio.plot import show
from rasterio.enums import Resampling
from shapely.geometry import Point
from concurrent.futures import ProcessPoolExecutor
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
testpack = "testfiles/testPack.gpkg"
testlayer = "chm_points"
testfield = "chm_height"
test_qgis = QgsVectorLayer(testpack+"|layername="+testlayer,testlayer,"ogr")
data_df = gpd.read_file(testpack,layer = testlayer)


new_df = valuesByDistance(data_df,'height',30)


#new_df['lstat']= new_df['list_height'].apply(lambda x: getLstats(x))
new_df['lcv']= new_df['lstat'].apply(lambda x: x[0])
new_df['lskew']= new_df['lstat'].apply(lambda x: x[1])
new_df = new_df.drop('lstat',axis=1)

new_df['list_height']=new_df['list_height'].apply(lambda x: str(x).replace('[','').replace(']','').replace('\'', '').replace(' ',''))
new_df['distances']=new_df['distances'].apply(lambda x: str(x).replace('[','').replace(']','').replace('\'', '').replace(' ',''))
new_df.to_file(testpack,layer='chm_points_lstat2',driver="GPKG")

#print (new_df)

# %%

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



# %%
testpack = "testfiles/testPack.gpkg"
testlayer = "trees"
test_qgis = QgsVectorLayer(testpack+"|layername="+testlayer,testlayer,"ogr")
data_df = gpd.read_file(testpack,layer = testlayer)

new = mixedPatchIndex(data_df,'height','ndvi_mean2',30)

new['list_height']=new['list_height'].apply(lambda x: str(x).replace('[','').replace(']','').replace('\'', '').replace(' ',''))
new['list_ndvi_mean2']=new['list_ndvi_mean2'].apply(lambda x: str(x).replace('[','').replace(']','').replace('\'', '').replace(' ',''))
new['distances']=new['distances'].apply(lambda x: str(x).replace('[','').replace(']','').replace('\'', '').replace(' ',''))

new.to_file(testpack,layer='mpitest4',driver="GPKG")

#new

# %%


# %%
def localConnectivity(input:QgsVectorLayer,field:str):
    print ("calculating local connecitity for the "+str(input.name())+" at field "+str(field))

# %%
def globalMoran(input:QgsVectorLayer):
    print ("calculating global Moran I")
# %%
"""
This python scirpt for Ripley's K function. Performance of the algoritmh needs improvement.

"""
from qgis.core import (
    QgsProject,
    QgsVectorLayer,
    QgsDistanceArea,
    QgsFeature
)
import numpy as np
import matplotlib.pyplot as plt

# Load the point layer
testpack = "testfiles/testPack.gpkg"
testlayer = "testtreemap"
testfield = "chm_height"
layer = QgsVectorLayer(testpack+"|layername="+testlayer,testlayer,"ogr")
#layer = QgsVectorLayer('/path/to/your/point_layer.shp', 'points', 'ogr')
if not layer.isValid():
    print("Layer failed to load!")

# Calculate the average density of points
num_points = layer.featureCount()
extent = layer.extent()
area = extent.width() * extent.height()
density = num_points / area

# Function to calculate distance between two points
def calculate_distance(point1, point2):
    d = QgsDistanceArea()
    return d.measureLine(point1, point2)

# Calculate Ripley's K function
distances = np.linspace(0, 1000, 100)  # Adjust the range and step as needed
K_values = []

for d in distances:
    count = 0
    for feature1 in layer.getFeatures():
        point1 = feature1.geometry().asPoint()
        for feature2 in layer.getFeatures():
            if feature1.id() != feature2.id():
                point2 = feature2.geometry().asPoint()
                if calculate_distance(point1, point2) <= d:
                    count += 1
    K = count / (num_points * density)
    K_values.append(K)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(distances, K_values, label="Ripley's K function")
plt.plot(distances, np.pi * distances**2, 'r--', label="CSR (Complete Spatial Randomness)")
plt.xlabel('Distance')
plt.ylabel("K(d)")
plt.legend()
plt.title("Ripley's K function")
plt.show()

# %%



# %%
def statisticInValue(slicevalues):
    #slicevalues = np.extract(zone==value,ref)
    mean = np.mean(slicevalues)
    std = np.std(slicevalues)
    count = len(slicevalues)
    maxvalue = np.max(slicevalues)
    minvalue = np.min(slicevalues)
        
    return mean,std,count,maxvalue,minvalue
    
# %%

def simple_treesegmentation2(chm,minheight,treesegments):
    """
    This produce treesegments using simple watershed algorithm
    Tree segmentation is based on Kaartinen et al. (2012) study
    Parameters:
    -chm: canopy height model in 2d numpy array format
    -peakarray: tree tops in 2d numpy array format. Tree tops need to have unique values and other are 0. Array need to be same size with chm.

    Outputs
    -labels: treesegments as labeled with peakarray values
    """

    # Define treesegments as output
    if treesegments is None:
        treesegments = tempfile.NamedTemporaryFile(prefix='resample_',suffix='.tif')
        treesegments = str(treesegments.name)
    
    chm = gdal.Open(chm)
    chmB = chm.GetRasterBand(1)
    chmA = chmB.ReadAsArray()
    transform = chm.GetGeoTransform()

    kernel = gaussian_kernel(size=3,sigma=1)
    chmS = convolve_2d(chmA,kernel=kernel)
    #print (transform.a,transform.c,transform.e,transform.f)
    fmax = calcFocal(chmS,2)
    peakarray = fmax - chmS
    peakarray = np.where((peakarray==0) & (chmA>minheight),chmA,0)
    peakarray = peakarray[0]

   
    # Find labels
    indices = np.argwhere(peakarray>0)
    #print (indices)
    
    # For centering point to middle of pixel
    cpx = transform[1] / 2
    cpy = transform[5] / 2
   
    # Replace each peak with a unique values
    peaks_xy = []
    for i, (row,col) in enumerate(indices):
        peak = (i+1,peakarray[row][col],transform[0] +transform[1] * col+cpx,transform[3] +transform[5] * row+cpy)
        peaks_xy.append(peak)
        peakarray[row,col] = i+1

    #convert peakarray to vector points
    peaks_gdf = gpd.GeoDataFrame({'zone': [i for i,v,x,y in peaks_xy],
                                'height': [round(v,2) for i,v,x,y in peaks_xy],
                                'geometry': [Point(x,y) for i,v,x,y in peaks_xy]})
    
    # Set the CRS to match the raster's CRS
    crs = chm.GetProjection()

    peaks_gdf.set_crs(crs, inplace=True)
    # Directions for rolling (up, down, left, right)
    shifts = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    image = -chmS
    #heights  = image
    labels = peakarray
    #directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    #labels = np.argwhere(labels>0)
    #front = [(i,j) for i,j in front]
    changed = True
    while changed:
        changed = False
        for shift in shifts:
            rolled_labels = np.roll(labels, shift, axis=(0, 1))
            rolled_image = np.roll(image,shift,axis=(0,1))
            mask = (rolled_labels > 0) & (labels == 0) & (abs(image)>minheight)
            labels[mask] = rolled_labels[mask]
            changed |= mask.any()  # Check if any changes were made

    gdal_array.SaveArray(labels.astype("int32"),treesegments,"GTiff",chm)

    chm = None
    chmB = None
    chmA = None

    return treesegments,peaks_gdf



# %%
testCHM = "testfiles/testCHM.tif"
testRe = "testfiles/testResamp.tif"
testResample = "testfiles/watershed_resample.tif"
treeSegments = "testfiles/watershed_result5.tif"
testOrto = "testfiles/testOrto.tif"
testNDVI = "testfiles/testNDVI23.tif"
testpack = "testfiles/testPack.gpkg"

#test = resample_raster(testCHM,None,2,'nearest')
#test = gdal.Open(test)
#print (test)
#test = simple_treesegmentation2(testCHM,2,treeSegments)

#test = zonal_stastics2(testCHM,treeSegments,1)
#print (test)
trees = singleTreeMapping(testCHM,testOrto,2,None)
#trees.to_file(testpack,layer='trese24',driver="GPKG")
print (trees)
#
# test = calculateNDVI2(testOrto,testNDVI)

#simple_treesegmentation(testCHM,2,treeSegments)


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
test = "C:\Users\mjkesala\AppData\Local\Temp\resample_hwocq8oc.tif"

# %%
temp = tempfile.NamedTemporaryFile(prefix='pre_', suffix='.tif')

gdal.Open(str(temp.name))
#print(temp.name)
# %%
