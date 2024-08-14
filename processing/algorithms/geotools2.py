from osgeo import gdal,gdal_array
import os,tempfile
from math import sqrt
import pandas as pd
import numpy as np
import rasterio
from rasterio.enums import Resampling
from shapely.geometry import Point
import geopandas as gpd
from qgis.core import QgsVectorLayer,QgsField,QgsFields,QgsFeature,QgsGeometry
from PyQt5.QtCore import QVariant

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
    
    print ("detecting individual trees from CHM and false-color ortophoto")

    
    with rasterio.open(chm) as src:
        chmA = src.read(1)
        #cs = src.cellsize
        transform = src.transform
    
    # smoothing chm
    kernel = gaussian_kernel(size=3,sigma=1)
    chmS = convolve_2d(chmA,kernel=kernel)
    
    #get peaks by focal analysis
    fmax = calcFocal(chmS,2)
    peakarray = fmax - chmS
    peakarray = np.where((peakarray==0) & (chmA>minheight),chmA,np.NaN)

    peaks = np.argwhere(peakarray>0)
    peakarray = peakarray[0]
    
    #convert peakarray to vector points
    peak_xy = np.array([(peakarray[y][x],transform.c +transform.a * x+0.5,transform.f +transform.e * y-0.5) for v,y,x in peaks])
    peaks_gdf = gpd.GeoDataFrame({'geometry': [Point(x,y) for v,x,y in peak_xy], 'height': [round(v,2) for v,x,y in peak_xy]})
    
    # Set the CRS to match the raster's CRS
    peaks_gdf.set_crs(src.crs.to_string(), inplace=True)
    
    return peaks_gdf

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
    if treesegments is None:
        treesegments = tempfile.NamedTemporaryFile(prefix='resample_',suffix='.tif')
        treesegments = str(treesegments.name)
    
    with rasterio.open(chm) as src:
        chmA = src.read(1)
        #cs = src.cellsize
        transform = src.transform
        crs = src.crs
    
    #chmS = chmA
    kernel = gaussian_kernel(size=3,sigma=1)
    chmS = convolve_2d(chmA,kernel=kernel)
    #print (transform.a,transform.c,transform.e,transform.f)
    fmax = calcFocal(chmS,2)
    peakarray = fmax - chmS
    peakarray = np.where((peakarray==0) & (chmA>minheight),1,0)
    peakarray = peakarray[0]
    # Find labels
    indices = np.argwhere(peakarray==1)
    #print (indices)
   
    # Replace each peak with a unique values
    for i, (row,col) in enumerate(indices):
        peakarray[row,col] = i+1

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
    with rasterio.open(
        treesegments,
        'w',
        driver='GTiff',
        height=labels.shape[0],
        width=labels.shape[1],
        count=1,
        dtype=labels.dtype,
        crs=crs,
        transform=transform,
        ) as dst:
            dst.write(labels, 1)

    return treesegments

def resample_raster(input_raster, output_path, new_cell_size,method=Resampling.nearest):
    """
    This performs raster sampling by using rasterio module

    Inputs:
    input_raster: input_raster as raster format
    output_path: output filename of the resampled raster
    new_cell_size: cell size of resampled raster as float or int format
    method: remapling method. Neareast is default. You can change it to the other rasterio sampling method.

    """
    with rasterio.open(input_raster) as src:
        # Calculate the new transform and dimensions
        transform = src.transform
        
        new_transform = transform * transform.scale(
            new_cell_size,
            new_cell_size
        )

        new_width = int((src.width * transform[0]) / new_cell_size)
        new_height = int((src.height * abs(transform[4])) / new_cell_size)

        # Resample the raster data
        data = src.read(
            out_shape=(src.count, new_height, new_width),
            resampling=method
        )

        # Update metadata
        profile = src.profile
        profile.update({
            'transform':new_transform,
            'width': new_width,
            'height': new_height
        })

        # Write the resampled raster to a new file
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(data)

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
    
    with rasterio.open(refraster) as ref:
        ref_a = ref.read(refband)
        
    with rasterio.open(zoneraster) as zone:
        zone_a = zone.read(1)
        
    if ref.transform[0]!=zone.transform[0]:
        temp = tempfile.NamedTemporaryFile(prefix='resample_',suffix='.tif')
        temp = str(temp.name)
        zoneresample = resample_raster(zoneraster,temp,ref.transform[0])
        #temp.close()
        with rasterio.open(temp) as zoner:
            zoner_a = zoner.read(1)

        zonelayer = zoner_a
    else:
        zonelayer = zone_a
    
    reflayer = ref_a
    
    indices = []
    meanl=[]
    stdl=[]
    countl=[]
    maxlist = []
    minlist = []

    for i in range(1,np.max(zonelayer)+1):
        slicevalues = np.extract(zonelayer==i,reflayer)
        mean = np.mean(slicevalues)
        std = np.std(slicevalues)
        count = len(slicevalues)
        maxvalue = np.max(slicevalues)
        minvalue = np.min(slicevalues)
        
        indices.append(i)
        meanl.append(mean)
        stdl.append(std)
        countl.append(count)
        maxlist.append(maxvalue)
        minlist.append(minvalue)
    
    zonalstastics = pd.DataFrame({"zoneindex":indices,
                                "mean":meanl,
                                "std":stdl,
                                "count":countl,
                                "max":maxlist,
                                "min":minlist})
    
    return zonalstastics

def calculateNDVI(falsecolor_orto,output_ndvi):
    """
    Calculating NDVI from falsecolor orthophoto
    Paremeters
    - falsecolor_orto: Colorinfrared (CIR) aerial orthophoto

    Outputs
    - NDVI: normalized  difference vegetation index (NDVI)
    """
    if output_ndvi is None:
        output_ndvi = tempfile.NamedTemporaryFile(prefix='ndvi_',suffix='.tif')
        output_ndvi = str(output_ndvi.name)
        
    ndvif = lambda c,r:  (c - r) / (c + r)

    with rasterio.open(falsecolor_orto) as src:
        nir = src.read(1)
        red = src.read(2)
        crs = src.crs
        transform = src.transform
    
    #if nir - red 
    nir = nir.astype('f')
    red = red.astype('f')
    ndvia = ndvif(nir,red) 

    with rasterio.open(
        output_ndvi,
        'w',
        driver='GTiff',
        height=ndvia.shape[0],
        width=ndvia.shape[1],
        count=1,
        dtype=ndvia.dtype,
        crs=crs,
        transform=transform,
        ) as dst:
            dst.write(ndvia, 1)
    
    return output_ndvi

def singleTreeMapping(chm,orto,minheight=2):
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
    
    print ("detecting individual trees from CHM and false-color ortophoto")

    treepoints = treeTops(chm,minheight=minheight)
    ts = simple_treesegmentation(chm,minheight,None)
    ndvi = calculateNDVI(orto,None)

    ndviz = zonal_stastics(ndvi,ts,1)
    ndviz = ndviz.rename(columns={"mean": "ndvi_mean", "std": "ndvi_std","count":"crownsize","max":"ndvi_max","min":"ndvi_min"})
    
    chm = zonal_stastics(chm,ts,1)
    chm = chm.rename(columns={"mean": "chm_mean", "std": "chm_std","max":"chm_max","min":"chm_min"})
    chm = chm.drop(columns=["count"])
    #blue = zonal_stastics(orto,ts,3)

    # Get the coordinates of the points
    point_coords = [(point.x, point.y) for point in treepoints.geometry]
    sa = rasterio.open(ts)

    # Sample the raster at these coordinates
    raster_values = list(sa.sample(point_coords))

    # Flatten the raster values if it's a single-band raster
    raster_values = [value[0] for value in raster_values]
    treepoints['zone'] = raster_values

    treepoints = treepoints.join(ndviz.set_index('zoneindex'), on='zone')
    treepoints = treepoints.join(chm.set_index('zoneindex'), on='zone')
    
    # Calculate crownsize to m2 unit
    with rasterio.open(orto) as src:
        transform = src.transform
        #arr = src.read(1)
    pixel_size = transform[0]
    treepoints["crownsize"] = pixel_size**2 * treepoints["crownsize"]
    

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