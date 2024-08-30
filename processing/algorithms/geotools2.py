from osgeo import gdal,gdal_array
import os,tempfile
from math import sqrt
import pandas as pd
import numpy as np
from shapely.geometry import Point
import geopandas as gpd
from qgis.core import QgsVectorLayer,QgsField,QgsFields,QgsFeature,QgsGeometry,QgsProcessingFeedback,QgsMessageLog
from qgis.core import Qgis
from PyQt5.QtCore import QVariant
from datetime import datetime

def calcFocal(in_array,distance):
    """
    This loops raster array and looks all cell values of each cell which are within distance values from cell
    """

    dat =pd.DataFrame(in_array)
    vert = dat
    ijlist = []
    for i in range(0-distance,distance):
        for j in range(0-distance,distance):
            e = sqrt(pow(i,2)+pow(j,2))
            if e <=distance:
                ijlist.append((i,j))
    
    #print (ijlist)

    for i in ijlist:
        df = dat.shift(i[0],axis=0)
        df = df.shift(i[1],axis=1)
        #vert = pd.concat([vert,df]).max(level=0)
        vert = np.maximum(df,vert)
    t = []
    t.append(vert)
    t = np.array(t)
    del dat,vert,df

    return t

def convolve_2d(image, kernel):
    """
    Perform a 2D convolution.
    
    Parameters:
    - image: 2D numpy array, the input image (raster data)
    - kernel: 2D numpy array, the convolution kernel
    
    Returns:
    - output: 2D numpy array, the convolved image
    """
    kernel_height, kernel_width = kernel.shape
    pad_height = kernel_height // 2
    pad_width = kernel_width // 2
    
    # Pad the image with zeros on all sides
    padded_image = np.pad(image, ((pad_height, pad_height), (pad_width, pad_width)), mode='constant')
    
    # Initialize the output array
    output = np.zeros_like(image)
    
    # Perform the convolution
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            # Extract the region of interest
            region = padded_image[i:i + kernel_height, j:j + kernel_width]
            # Perform element-wise multiplication and sum the result
            output[i, j] = np.sum(region * kernel)
    del image,padded_image

    return output

def gaussian_kernel(size, sigma):
    """
    Create a 2D Gaussian kernel.
    
    Parameters:
    - size: int, the size of the kernel (size x size)
    - sigma: float, the standard deviation of the Gaussian distribution
    
    Returns:
    - kernel: 2D numpy array, the Gaussian kernel
    """
    feedback = QgsProcessingFeedback()
    feedback.pushConsoleInfo("Implementing gaussian kernel filtering for CHM")

    kernel = np.fromfunction(
        lambda x, y: (1 / (2 * np.pi * sigma**2)) * np.exp(
            -((x - (size - 1) / 2)**2 + (y - (size - 1) / 2)**2) / (2 * sigma**2)
        ), 
        (size, size)
    )
    return kernel / np.sum(kernel)

def treeTops(chm,minheight=2):
    """
    This detect individual tree tops from CHM

    Parameters:
    - chm: canopy heigh model
    - minheight: minimum height of tree to be considering

    Returns:
    - peaks: peaks as vector points of geopandas dataframe 
    """

    chm = gdal.Open(chm)
    
    chmB = chm.GetRasterBand(1)
    chmA = chmB.ReadAsArray()
    
    transform = chm.GetGeoTransform()
    
    # smoothing chm
    kernel = gaussian_kernel(size=3,sigma=1)
    chmS = convolve_2d(chmA,kernel=kernel)
    
    #get peaks by focal analysis
    fmax = calcFocal(chmS,2)
    peakarray = fmax - chmS
    peakarray = np.where((peakarray==0) & (chmA>minheight),chmA,np.NaN)

    peaks = np.argwhere(peakarray>0)
    peakarray = peakarray[0]
    #print (peaks)
    
    #convert peakarray to vector points
    peak_xy = np.array([(peakarray[y][x],transform[0] +transform[1] * x+0.5,transform[3] +transform[5] * y-0.5) for v,y,x in peaks])
    peaks_gdf = gpd.GeoDataFrame({'geometry': [Point(x,y) for v,x,y in peak_xy], 'height': [round(v,2) for v,x,y in peak_xy]})
    
    # Set the CRS to match the raster's CRS
    crs = chm.GetProjection()

    peaks_gdf.set_crs(crs, inplace=True)
    chm = None
    chmB = None
    chmA = None

    return peaks_gdf

def optimized_treesegmentation(chm,minheight,treesegments):
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


def simple_treesegmentation(chm,minheight,treesegments):
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

    image = -chmS
    labels = peakarray
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    front = np.argwhere(labels>0)
    front = [(i,j) for i,j in front]
    while front:
        new_front = []
        for (i, j) in front:
            treeheight = abs(image[i,j])
            for direction in directions:
                ni, nj = i + direction[0], j + direction[1]
                if ni >= 0 and ni < image.shape[0] and nj >= 0 and nj < image.shape[1]:
                    if labels[ni, nj] == 0 and image[ni, nj] != 0 and abs(image[ni,nj]) > 1 and abs(image[ni,nj])>treeheight /2:
                        labels[ni, nj] = labels[i, j]
                        new_front.append((ni, nj))
        front = new_front
    
    gdal_array.SaveArray(labels.astype("int32"),treesegments,"GTiff",chm)

    chm = None
    chmB = None
    chmA = None

    return treesegments,peaks_gdf


def resample_raster(input,output,pixel_size:int,method):
    """
    This performs raster sampling by using gdal library

    Inputs:
    input_raster: input_raster as raster format
    output_path: output filename of the resampled raster
    new_cell_size: cell size of resampled raster as float or int format
    method: remapling method. Neareast is default. You can change it to the other rasterio sampling method.

    """

    # Define output
    if output is None:
        output = tempfile.TemporaryFile(prefix='resample_',suffix='.tif')
        output = output.name
        #output = os.path()
    # Open the input raster
    dataset = gdal.Open(input)
    
    # Define the output resolution
    x_res = pixel_size  # Pixel size in x direction (e.g., 10 meters)
    y_res = pixel_size  # Pixel size in y direction (e.g., 10 meters)

    # Resample the raster
    gdal.Warp(
        output,
        dataset,
        xRes=x_res,
        yRes=y_res,
        resampleAlg=method  # Use 'nearest', 'bilinear', 'cubic', etc.
        )
    #Close the dataset
    dataset=None

    return output

def zonal_stastics(refraster,zoneraster,refband=1):
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
    zones = np.unique(zonelayer[zonelayer != nodata_value])
    zonera = None
    refra = None
    
    # Calculate statistics for each zone
    zonal_statistics = {}
    for zone in zones:
        zone_values = np.extract(zonelayer==zone,reflayer)
        #mask = zone_a == zone  # Create a mask for the current zone
        #zone_values = reflayer[mask]  # Extract values from the reference raster for this zone
        
        # Calculate statistics
        mean = np.mean(zone_values)
        std = np.std(zone_values)
        count_values = len(zone_values)
        min_value = np.min(zone_values)
        max_value = np.max(zone_values)
        
        # Store results in a dictionary
        zonal_statistics[zone] = {
            "zoneindex":zone,
            "mean": mean,
            "std":std,
            "count": count_values,
            "min": min_value,
            "max": max_value
        }
    zonal_statistics = pd.DataFrame(zonal_statistics).T
  
    return zonal_statistics


def calculateNDVI(falsecolor_orto,output_ndvi):
    """
    Calculating NDVI from falsecolor orthophoto
    Paremeters
    - falsecolor_orto: Colorinfrared (CIR) aerial orthophoto

    Outputs
    - NDVI: normalized  difference vegetation index (NDVI)
    """
    # Define output_ndvi    
    if output_ndvi is None:
        output_ndvi = tempfile.NamedTemporaryFile(prefix='ndvi_',suffix='.tif')
        output_ndvi = str(output_ndvi.name)
    
    # Step 1: Open the dataset
    dataset = gdal.Open(falsecolor_orto)

    # Step 2: Read the NIR and Red bands
    # Assuming that NIR is in band 4 and Red is in band 3 (this may vary)
    nir_band = dataset.GetRasterBand(1)
    red_band = dataset.GetRasterBand(2)

    # Step 3: Convert the bands to numpy arrays
    nir = nir_band.ReadAsArray().astype(float)
    red = red_band.ReadAsArray().astype(float)

    # Step 4: Calculate NDVI
    ndvi = (nir - red) / (nir + red)

    # Step 5: Set up the output image
    driver = gdal.GetDriverByName('GTiff')
    out_raster = driver.Create(output_ndvi, dataset.RasterXSize, dataset.RasterYSize, 1, gdal.GDT_Float32)

    # Step 6: Set the geo-transform and projection
    out_raster.SetGeoTransform(dataset.GetGeoTransform())
    out_raster.SetProjection(dataset.GetProjection())

    # Step 7: Write the NDVI band
    out_band = out_raster.GetRasterBand(1)
    out_band.WriteArray(ndvi)

    # Step 8: Set NDVI band to have NoData value
    out_band.SetNoDataValue(-9999)

    # Step 9: Clean up and close datasets
    out_band.FlushCache()
    out_raster = None
    dataset = None

    return output_ndvi

def singleTreeMapping(chm,orto,minheight,feedback):
    """
    This detect individual trees by minimum curvature method
    
    step 1: calculate tree tops
    step 2: calculate treesegments
    step 3: calculate NDVI
    step 4: calculate zonal statistic of NDVI and CHM within treesegment zones
    step 5: join treesegment id to tree tops by rastersampling
    step 6: join zonal stastic to point

    Parameteres:
    chm: canopy heigh model as raster format
    orto: false-color otrhophoto (band1 = infrared band2 = red band3 = green)
    minheight: minimum height of tree (m) that will convert to single-tree point

    outputs:
    treepoints: single-tree points with statistic of NDVI and CHM
    """
    if feedback is None:
        feedback = QgsProcessingFeedback()

    now = datetime.now().strftime("%H:%M:%S")
    feedback.pushInfo(now+"\tDetecting individual trees from CHM and segmenting tree canopies")
    
    ts = optimized_treesegmentation(chm,minheight,None)
    
    treepoints = ts[1]
    ts = ts[0]

    now = datetime.now().strftime("%H:%M:%S")
    feedback.pushInfo(now+"\tCalculating NDVI")
    ndvi = calculateNDVI(orto,None)

    now = datetime.now().strftime("%H:%M:%S")
    feedback.pushInfo(now+"\tCalculating zonal statistics of CHM and NDVI")
    ndviz = zonal_stastics(ndvi,ts,1)
    ndviz = ndviz.rename(columns={"mean": "ndvi_mean", "std": "ndvi_std","count":"crownsize","max":"ndvi_max","min":"ndvi_min"})
    
    chmz = zonal_stastics(chm,ts,1)
    chmz = chmz.rename(columns={"mean": "chm_mean", "std": "chm_std","max":"chm_max","min":"chm_min"})
    chmz = chmz.drop(columns=["count"])

    now = datetime.now().strftime("%H:%M:%S")
    feedback.pushInfo(now+"\tCleaning results and writing data to output")

    treepoints = treepoints.join(ndviz.set_index('zoneindex'), on='zone')
    treepoints = treepoints.join(chmz.set_index('zoneindex'), on='zone')
    
    # Calculate crownsize to m2 unit
    src = gdal.Open(orto)
    src2 = gdal.Open(chm)
    pixel_size=np.max([src.GetGeoTransform()[1],src2.GetGeoTransform()[1]])

    treepoints["crownsize"] = pixel_size**2 * treepoints["crownsize"]
    
    src = None
    src2 = None

    return treepoints

def gpd2qgis(gdf):
    """
    This convert geopandas dataframe to QgsVectorLayer
    Parameters
    gdf: geopandas dataframe with geometry field

    returns memory layer in QgsVectorLayer format
    """

    # Define mapping from GeoPandas dtypes to PyQt QVariant types
    dtype_to_qvariant = {
        'int64': QVariant.Int,
        'float64': QVariant.Double,
        'object': QVariant.String,  # This usually corresponds to strings
        'bool': QVariant.Bool,
        # Add other mappings as needed
    }
    
    # Determine the coordinate system
    if gdf.crs is None:
        gdf.set_crs('EPSG:4326', inplace=True)
    crs = str(gdf.crs)
    
    # Determine the geometry type
    geometry_type = gdf.geom_type.iloc[0]
    if geometry_type in ('Point','LineString','Polygon'):
        qgis_type = geometry_type+"?crs="+crs
    else:
        raise ValueError(f"Unsupported geometry type: {geometry_type}")
    
    # Create a memory layer
    layer = QgsVectorLayer(qgis_type, "temp_layer", "memory")
    provider = layer.dataProvider()
    
    # Loop through GeoPandas DataFrame columns
    fields = QgsFields()
    for col_name in gdf.columns:
        if col_name == 'geometry':
            continue  # Skip the geometry column
        
        # Determine the type of the field
        col_type = gdf[col_name].dtype.name
        
        # Get the corresponding QVariant type
        qvariant_type = dtype_to_qvariant.get(col_type, QVariant.String)
        
        # Create and add the QgsField
        field = QgsField(col_name, qvariant_type)
        fields.append(field)
    
    provider.addAttributes(fields)
    layer.updateFields()
    
    for index, row in gdf.iterrows():
        #   Create a new QgsFeature
        feature = QgsFeature()
        feature.setGeometry(QgsGeometry.fromWkt(row['geometry'].wkt))
        feature.setAttributes([row[col] for col in gdf.columns if col != 'geometry'])
        provider.addFeature(feature)

    layer.updateExtents()
    
    return layer

def valuesByDistance(gdf,fieldname:str,radius:int):
    """
    This lists values and distances of nearby features that fall in given radius from the feature

    Paramerets:
        gdf : geopandas dataframe
        fieldname: the field name of geodataframe whose values ​​to list
        radius : search radius of features
    Outputs:
        gdf: geodataframe with new fields (list_<fieldname> and distances)
    """
    # Buffering and spatial join
    gdf['buffer'] = gdf.geometry.buffer(radius)
    gdf['geom2'] = gdf['geometry']
    joined_gdf = gpd.sjoin(gdf, gpd.GeoDataFrame(gdf[['buffer',fieldname,'geom2']]).set_geometry('buffer'), how='left', op='intersects')
    
    # clean up columns
    joined_gdf = joined_gdf.rename(columns={fieldname+'_left':fieldname,fieldname+'_right':'nearvalue','geom2_right':'geom2'})

    # Compute distances
    joined_gdf['distance'] = joined_gdf.geometry.distance(gpd.GeoDataFrame(joined_gdf[['geom2']]).set_geometry('geom2'))
    joined_gdf['distance'] = joined_gdf['distance'].round(2)
    
    # Aggregating the attribute values
    aggregated = joined_gdf.groupby(joined_gdf.index).agg({
        'nearvalue': lambda x: list(x),
        'distance': lambda x: list(x)})
    
    # Passing fieldvalues to original and cleaning fields
    gdf['list_'+fieldname] = gdf.index.map(aggregated['nearvalue'])
    gdf['distances'] = gdf.index.map(aggregated['distance'])
    gdf = gdf.drop(['buffer','geom2'],axis=1)
    
    return gdf

def layer2gpd(layer):
    """
    This convert QGIS layer to geopandas dataframe
    Parameters
        layer: qgis vector layer (QgsVectorLayer)

    Outputs
        gdf: geopandas dataframe
    """

    features = layer.getFeatures()
    
    # Extract geometries and attributes
    crs = layer.crs().authid()
    data = []
    for feature in features:
        geom = feature.geometry().asWkt()
        attrs = feature.attributes()
        data.append({'geometry': geom, **dict(zip(layer.fields().names(), attrs))})
    
    # Convert to GeoPandas GeoDataFrame
    gdf = gpd.GeoDataFrame(data, geometry=gpd.GeoSeries.from_wkt([f['geometry'] for f in data]))
    gdf = gdf.set_crs(crs)
    
    return gdf
