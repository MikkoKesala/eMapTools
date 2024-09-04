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
from geotools2 import gpd2qgis,valuesByDistance,singleTreeMapping,resample_raster,optimized_treesegmentation
from geostats import getLstats
import matplotlib.pyplot as plt
import rasterio
import tempfile
from decimal import Decimal
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
def statisticInValue(zone,zone_values):
    mean = np.mean(zone_values)
    std = np.std(zone_values)
    count_values = len(zone_values)
    min_value = np.min(zone_values)
    max_value = np.max(zone_values) 
    
    return mean,std,count_values,min_value,max_value
    
# %%

def zonal_stastics2(refraster,zoneraster,refband=1):
    """
    This calculates zonal statistic based on raster zones

    Inputs:
    refraster: reference raster which statistic are calculating within zones
    zoneraster: zone areas as raster format. If the cellsize is not same with reference layer, this performs rastersampling
    refband: band number of the refenrece layer. First band is number 1.

    Outputs:
    zonalstastics: zonalstatisc as pandas dataframe format with zoneindex from the zonelayer.
    """
    refra = gdal.Open(refraster)
    reflayer = refra.GetRasterBand(refband).ReadAsArray()
    
    zonera = gdal.Open(zoneraster)
    zonelayer = zonera.GetRasterBand(1).ReadAsArray()
    nodata_value = zonera.GetRasterBand(1).GetNoDataValue()

    
    reftransform = refra.GetGeoTransform()
    zonetransform = zonera.GetGeoTransform()
    
    if reftransform[1]!=zonetransform[1]:
        if reftransform[1]<zonetransform[1]:
            temp = resample_raster(refraster,None,zonetransform[1],'average')
            refra = gdal.Open(temp)
            reflayer = refra.GetRasterBand(refband).ReadAsArray()

        else:
            temp = resample_raster(zoneraster,None,reftransform[1],'nearest')
            zonera = gdal.Open(temp)
            zonelayer = zonera.GetRasterBand(1).ReadAsArray()
            


    # List of unique zones
    zones = np.unique(zonelayer[(zonelayer != nodata_value) & (zonelayer>0)])
    zonera = None
    refra = None
    
    decimal_places = np.max([abs(Decimal(str(reflayer[i][i])).as_tuple().exponent) for i in range(100)])
    if decimal_places >3:
        decimal_places==3
    #cf = np.power(10,decimal_places)
    #reflayer = np.int32(reflayer*cf)
    # Calculate statistics for each zone
    zonal_statistics = {}
    for zone in zones:
        #zone_values = np.extract(zonelayer==zone,reflayer)
        #mask = zone_a == zone  # Create a mask for the current zone
        zone_values = reflayer[zonelayer==zone]  # Extract values from the reference raster for this zone
        #reflayer = np.where(zonelayer!=zone,reflayer,np.NaN)
        #statvalues = np.fromfunction(statisticInValue,zone_values)
        # Calculate statistics
        mean = np.mean(zone_values)
        std = np.std(zone_values)
        count_values = len(zone_values)
        min_value = np.min(zone_values)
        max_value = np.max(zone_values) 
        
        # Store results in a dictionary
        zonal_statistics[zone] = {
            "zoneindex":zone,
            "mean":mean,
            "std":std,
            "count": count_values,
            "min": min_value,
            "max": max_value
        }
    zonal_statistics = pd.DataFrame(zonal_statistics).T
  
    return zonal_statistics



# %%
testCHM = "testfiles/testCHM.tif"
testRe = "testfiles/testResamp.tif"
testResample = "testfiles/watershed_resample.tif"
treeSegments = "testfiles/treesegments_big.tif"
testOrto = "testfiles/testOrto.tif"
testNDVI = "testfiles/testNDVI23.tif"
testpack = "testfiles/testPack.gpkg"

orto = 'Z:/aineistoja/TreeMapping/L3334HCIR_1m.tif'
chm = 'Z:/aineistoja/TreeMapping/L3334HCHM.tif'

#test = resample_raster(testCHM,None,2,'nearest')
#test = gdal.Open(test)
#print (test)
test = optimized_treesegmentation(chm,2,treeSegments)
trees = test[1]
#test = zonal_stastics2(testCHM,treeSegments,1)
#print (test)
#trees = singleTreeMapping(testCHM,testOrto,2,None)
#trees['height'] = np.int16(trees['height']*10) / 10
trees.to_file(testpack,layer='trees_big2',driver="GPKG")
#print (trees)
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
data_folder = r'C:/Users/mjkesala/OneDrive - University of Helsinki/QGIS/ReTreeT/Aineistot/Results2.gpkg'
scenarios = ['waterprotection','biodiversity','climate','equal','dwp','ds','dtw','pret','realized']

#points
all = gpd.read_file(data_folder,layer = 'all2', fid_as_index=True)


# %%
allNb = valuesByDistance(all,"CHM",10)
allNb
# %%
print (np.average([len(i) for i in allNb['list_CHM']]))
# %%
