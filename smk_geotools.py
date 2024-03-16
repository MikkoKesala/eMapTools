from osgeo import gdal,gdal_array
import os,sys,tempfile
from math import sqrt
import pandas as pd
import numpy as np
from qgis.PyQt.QtCore import QVariant
from qgis import processing
#import processing
from qgis.core import QgsVectorLayer,QgsField,QgsFeature,edit,QgsRectangle,QgsRasterLayer,QgsRasterPipe,QgsRasterProjector,QgsRasterFileWriter
from qgis.analysis import QgsInterpolator,QgsIDWInterpolator,QgsGridFileWriter

"""
#for developing
QgsApplication.setPrefixPath(QgsApplication.prefixPath(), True)
qgs = QgsApplication([], False)
qgs.initQgis()
sys.path.append(os.path.join(QgsApplication.prefixPath(),"python\plugins"))
#import processing
from processing.core.Processing import Processing
Processing.initialize()
#import processing
#from processing.core.Processing import Processing
#Processing.initialize()
#import processing
from qgis.analysis import QgsNativeAlgorithms"""


def feature2Layer(feat,buffer):
    """
    This creates vectorlayer from feature. You can buffer feature at specific distance
    """
    vl = QgsVectorLayer("Polygon", "temporary_points", "memory")
    pr = vl.dataProvider()

    # add fields
    vl.startEditing()
    pr.addAttributes(feat.fields())
    vl.updateFields() # tell the vector layer to fetch changes from the provider

    fet = QgsFeature()
    fet.setGeometry(feat.geometry().buffer(buffer,5))
    #fet.setGeometry(feat.geometry())
    fet.setAttributes(feat.attributes())
    pr.addFeatures([fet])

    vl.updateExtents()
    vl.commitChanges()

    return vl

def copyRaster2(inp,outp):
    """
    This is simple raster copy
    """
    os.popen('copy '+inp+' '+outp)


def focalMaximaCHM(input_raster,distance,convert):
    """
    This calculates focal maximum value by specific search distance
    """
    
    rastOut = input_raster[0:-4]+"hh.tif"
    chm = gdal.Open(input_raster)
    
    chmB = chm.GetRasterBand(1)
    chmA = chmB.ReadAsArray()
    if convert == True:
        chmA = -0.118*chmA+30.1567 #vaihe 1

    focal = calcFocal(chmA,distance) #vaihe 2
    huip = focal - chmA
    huip = np.where(focal-chmA==0,chmA,np.NaN)
    huip = np.where(huip>=5,huip*10,np.NaN) #vaihe3
    gdal_array.SaveArray(huip.astype("float32"),rastOut,"GTiff",chm)

    return rastOut

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
        vert = pd.concat([vert,df]).max(level=0)
    t = []
    t.append(vert)
    t = np.array(t)

    return t

def delNulls(input_vector):
    """
    This deletes null-values of treemap  
    """
    input_vector = QgsVectorLayer(input_vector, "puukartta", "ogr")
    

    feats = input_vector.getFeatures()
    dfeat=[]
    
    for feat in feats:
        if feat['CHM'] < 0:
            dfeat.append(feat.id())

    input_vector.dataProvider().deleteFeatures(dfeat)


def createTreeMap(input_chm,distance,convert):
    """
    This creates tree map as point layer from chm raster. Algorithm is based on local maxima at specific search distance
    """
    #tempd = tempfile.TemporaryFile()
    #tempd = tempd.name+'.shp'
    focalMax = focalMaximaCHM(input_chm,distance,convert)
    tempd = processing.run("gdal:polygonize", {'INPUT':focalMax,'BAND':1,'FIELD':'CHM','EIGHT_CONNECTEDNESS':False,'EXTRA':'','OUTPUT':'TEMPORARY_OUTPUT'})
    
    delNulls(tempd['OUTPUT'])
    tempd = processing.run("native:centroids", {'INPUT':tempd['OUTPUT'],'ALL_PARTS':False,'OUTPUT':'TEMPORARY_OUTPUT'})
    processing.run("native:createspatialindex", {'INPUT':tempd['OUTPUT']})

    return tempd['OUTPUT']

def addFieldValue(in_feat:QgsVectorLayer,fieldname:str,fieldvalue:float):
    """
    This adds field with given value to the vector layer
    """
    fix = processing.run("native:fixgeometries", {'INPUT':in_feat,'OUTPUT':'TEMPORARY_OUTPUT'})
    fix['OUTPUT'].dataProvider().addAttributes([QgsField(fieldname,QVariant.Double)])
    fix['OUTPUT'].updateFields()
    with edit(fix['OUTPUT']):
        for feat in fix['OUTPUT'].getFeatures():
            feat['leimikko']=fieldvalue

            fix['OUTPUT'].updateFeature(feat)
    
    return fix['OUTPUT']

def joinIntersection(inlayer,joinlayer,joinfields,drop):
    """
    This join by spatial intersection two layers
    """
    joined = processing.run("native:joinattributesbylocation", {'INPUT':inlayer,'JOIN':joinlayer,'PREDICATE':[0],'JOIN_FIELDS':joinfields,'METHOD':0,'DISCARD_NONMATCHING':drop,'PREFIX':'','OUTPUT':'TEMPORARY_OUTPUT'})

    return joined['OUTPUT']

def hsAnalysis(in_feat,fieldname):
    """
    This interpolate field value of point layer. Input layer need to be QgsVectorLayer format
    and fieldname string format. The analysis save interpolated values to new field of point layer with prefix 'HS_'
    """
    #in_feat = QgsVectorLayer(in_feat,"tt","ogr")

    tempd = tempfile.TemporaryFile()
    tempd = tempd.name+'.tif'
    in_feat.updateExtents()
    ext = in_feat.extent()
    idx = in_feat.dataProvider().fieldNameIndex(fieldname)

    layer_data = QgsInterpolator.LayerData()
    layer_data.source = in_feat 
    layer_data.zCoordInterpolation = False
    layer_data.interpolationAttribute = idx
    layer_data.mInputType = 1


    idw_interpolator = QgsIDWInterpolator([layer_data])
    res = 1.0
    
    ncols = int((ext.xMaximum() -ext.xMinimum()) / res )
    nrows = int((ext.yMaximum() - ext.yMinimum()) / res)

    out  = QgsGridFileWriter(idw_interpolator,tempd,ext,ncols,nrows)
    out.writeFile()

    out = processing.run("native:rastersampling", {'INPUT':in_feat,'RASTERCOPY':tempd,'COLUMN_PREFIX':'HS_','OUTPUT':'TEMPORARY_OUTPUT'})
 
    #hs = processing.run("qgis:idwinterpolation",{'INTERPOLATION_DATA':in_name+"::~::0::~::"+str(idx)+"::~::0",'DISTANCE_COEFFICIENT':2,'EXTENT':ext,'PIXEL_SIZE':1,'OUTPUT':'TEMPORARY_OUTPUT'})
    
    return out['OUTPUT']


def copyVector(layer):
    
    feats = [feat for feat in layer.getFeatures()]
    crs = str(layer.crs().authid())
    mem_layer = QgsVectorLayer("Point?crs="+crs, "copy", "memory")

    mem_layer_data = mem_layer.dataProvider()
    attr = layer.dataProvider().fields().toList()
    mem_layer_data.addAttributes(attr)
    mem_layer.updateFields()
    mem_layer_data.addFeatures(feats)

    return mem_layer

    
def point2area(input_points,fname,value):
    """
    This converts points to area. You can filter points specific field value.
    """
    copylayer = copyVector(input_points)
    
    if len(fname) > 0 and value is not None:
        NoLeim = [feat.id() for feat in copylayer.getFeatures() if feat[fname]!=value]
        
    
        copylayer.dataProvider().deleteFeatures(NoLeim)
        copylayer.updateFields()


    buf = processing.run("native:buffer", {'INPUT':copylayer,
                                    'DISTANCE':5,
                                    'SEGMENTS':5,
                                    'END_CAP_STYLE':0,
                                    'JOIN_STYLE':0,
                                    'MITER_LIMIT':2,
                                    'DISSOLVE':True,
                                    'OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']

    
    part = processing.run("native:multiparttosingleparts", {'INPUT':buf,'OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']

    buf1 = processing.run("native:buffer", {'INPUT':part,
                                    'DISTANCE':10,
                                    'SEGMENTS':5,
                                    'END_CAP_STYLE':0,
                                    'JOIN_STYLE':0,
                                    'MITER_LIMIT':2,
                                    'DISSOLVE':False,
                                    'OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']

    buf = processing.run("native:buffer", {'INPUT':buf1,
                                    'DISTANCE':-12,
                                    'SEGMENTS':5,
                                    'END_CAP_STYLE':0,
                                    'JOIN_STYLE':0,
                                    'MITER_LIMIT':2,
                                    'DISSOLVE':False,
                                    'OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']

    part = processing.run("native:multiparttosingleparts", {'INPUT':buf,'OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']
    
    return part

def clipRaster2(input_raster,clip_vector):
    xsize =  int(round(input_raster.rasterUnitsPerPixelX(),0))
    input_raster.setExtent(clip_vector.extent())
    #ysize = int(round(input_raster.rasterUnitsPerPixelY(),0))
    #orig_name = input_raster.name()
    #print (chm.name())
    #input_raster.setName("raster")
    alg_params = {
                'BURN': 1,
                'DATA_TYPE': 5,  # Float32
                'EXTENT': None,
                'EXTRA': '',
                'FIELD': '',
                'HEIGHT': 1,
                'INIT': None,
                'INPUT': clip_vector,
                'INVERT': False,
                'NODATA': 0,
                'OPTIONS': '',
                'UNITS': 0,  # Soluina (pikseleinä)
                'USE_Z': False,
                'WIDTH': 1,
                'OUTPUT':'TEMPORARY_OUTPUT'}

    #raster = processing.run('gdal:rasterize', alg_params)
    #raster = QgsRasterLayer(raster['OUTPUT'],"extent","gdal")
    raster_extent = roundExtent(clip_vector,0)

    alg_params = {
                'CELLSIZE': xsize,
                'CRS': None,
                'EXPRESSION':input_raster.name()+"@1",
                'EXTENT': raster_extent,
                'LAYERS': input_raster,
                'OUTPUT': 'TEMPORARY_OUTPUT'}

    raster_clip = processing.run('qgis:rastercalculator', alg_params)
    
    #nput_raster.setName(orig_name)
    
    return raster_clip['OUTPUT']

def roundExtent(layer,decimals:int):


    if layer.isValid():
        # Get the extent of the raster layer
        extent = layer.extent()

        # Set the number of decimal places to round to
        decimal_places = decimals  # Adjust as needed

        # Round the extent coordinates to the desired precision
        rounded_extent = QgsRectangle(
            round(extent.xMinimum(), decimal_places),
            round(extent.yMinimum(), decimal_places),
            round(extent.xMaximum(), decimal_places),
            round(extent.yMaximum(), decimal_places)
        )

        # Print the rounded extent
        print(f"Rounded Extent: {rounded_extent}")
    else:
        print("Invalid raster layer")

    return rounded_extent

def clipRaster3(rlayer,mask_layer):
    renderer = rlayer.renderer()
    provider = rlayer.dataProvider()
    crs = rlayer.crs()

    pipe = QgsRasterPipe()
    projector = QgsRasterProjector()
    projector.setCrs(provider.crs(), provider.crs())

    if not pipe.set(provider.clone()):
        print("Cannot set pipe provider")

    # Commented for extract raw data
    # if not pipe.set(renderer.clone()):
        # print("Cannot set pipe renderer")

    if not pipe.insert(2, projector):
        print("Cannot set pipe projector")

    
    out_file = tempfile.TemporaryFile()
    out_file = out_file.name+'.tif'
    #out_file = 'D:/temp/temporal.tif'
    file_writer = QgsRasterFileWriter(out_file)
    file_writer.Mode(1)

    print ("Saving")

    extent = mask_layer.extent()

    opts = ["COMPRESS=LZW"]
    file_writer.setCreateOptions(opts)
    error = file_writer.writeRaster(
        pipe,
        extent.width (),
        extent.height(),
        extent,
        crs)

    if error == QgsRasterFileWriter.NoError:
        print ("Raster was saved successfully!")
        #layer = QgsRasterLayer(out_file, "result")
        
    else:
        print ("Raster was not saved!")

    return out_file

def clipRaster4(rlayer,mask_layer,outfolder,rastername):
    renderer = rlayer.renderer()
    provider = rlayer.dataProvider()
    crs = rlayer.crs()

    pipe = QgsRasterPipe()
    projector = QgsRasterProjector()
    projector.setCrs(provider.crs(), provider.crs())

    if not pipe.set(provider.clone()):
        print("Cannot set pipe provider")

    # Commented for extract raw data
    # if not pipe.set(renderer.clone()):
        # print("Cannot set pipe renderer")

    if not pipe.insert(2, projector):
        print("Cannot set pipe projector")

    
    #out_file = tempfile.TemporaryFile()
    #out_file = out_file.name+'.tif'
    #out_file = 'D:/temp/temporal.tif'
    out_file = os.path.join(outfolder,"raster_"+str(rastername))
    file_writer = QgsRasterFileWriter(out_file)
    file_writer.Mode(0)

    print ("Saving")

    extent = mask_layer.extent()

    opts = ["COMPRESS=LZW"]
    file_writer.setCreateOptions(opts)
    error = file_writer.writeRaster(
        pipe,
        extent.width (),
        extent.height(),
        extent,
        crs)

    if error == QgsRasterFileWriter.NoError:
        print ("Raster was saved successfully!")
        #layer = QgsRasterLayer(out_file, "result")
        
    else:
        print ("Raster was not saved!")

    chm = gdal.Open(out_file)
    
    chmB = chm.GetRasterBand(1)
    chmA = chmB.ReadAsArray()
    
    chmA = -0.118*chmA+30.1567 #vaihe 1

    gdal_array.SaveArray(chmA.astype("float32"),out_file,"GTiff",chm)

    return out_file

def go_gaussian(input_points,id_field,weight_field,bandwidth):
    # Distance matrix
    #feedback = QgsProcessingFeedback
    #context = QgsProcessingContext
    #results = {}
    #outputs = {}
    alg_params = {
        'INPUT': input_points,
        'INPUT_FIELD': id_field,
        'MATRIX_TYPE': 0,  # Linear (N*k x 3) distance matrix
        'NEAREST_POINTS': 40,
        'TARGET': input_points,
        'TARGET_FIELD': weight_field,
        'OUTPUT': 'TEMPORARY_OUTPUT'
    }
    dmatrix = processing.run('qgis:distancematrix', alg_params)['OUTPUT']

    # Aggregate
    alg_params = {
        'AGGREGATES': [{'aggregate': 'first_value','delimiter': ',','input': '"InputID"','length': 0,'name': 'InputID','precision': 0,'sub_type': 0,'type': 4,'type_name': 'int8'},
                        {'aggregate': 'array_agg','delimiter': ',','input': '"TargetID"','length': 0,'name': 'TargetID','precision': 0,'sub_type': 10,'type': 11,'type_name': 'stringlist'},
                        {'aggregate': 'array_agg','delimiter': ',','input': '"Distance"','length': 0,'name': 'Distance','precision': 0,'sub_type': 10,'type': 11,'type_name': 'stringlist'}],
        'GROUP_BY': '"InputID"',
        'INPUT': dmatrix,
        'OUTPUT':'TEMPORARY_OUTPUT'
    }
    agg = processing.run('native:aggregate', alg_params)['OUTPUT']

    # Join attributes by field value
    alg_params = {
        'DISCARD_NONMATCHING': False,
        'FIELD': id_field,
        'FIELDS_TO_COPY': [''],
        'FIELD_2': 'InputID',
        'INPUT': input_points,
        'INPUT_2': agg,
        'METHOD': 1,  # Take attributes of the first matching feature only (one-to-one)
        'PREFIX': '',
        'OUTPUT': 'TEMPORARY_OUTPUT'
    }
    epoints = processing.run('native:joinattributestable', alg_params)['OUTPUT']
    #calculating kernel_density with weights
    #points = QgsProcessingUtils.mapLayerFromString(outputs['JoinAttributesByFieldValue']['OUTPUT'],context)
    epoints.dataProvider().addAttributes([QgsField('getisord',QVariant.Double)])
    epoints.updateFields()

    # Calculate global mean and variance
    all_values = [f[weight_field] for f in epoints.getFeatures() if type(f[weight_field]) in (float,int)]
    x_mean = np.mean(all_values)
    n = len(all_values)
    s = np.sqrt((np.sum(a**2 for a in all_values) / n)- (x_mean**2))
    b = bandwidth
    
    with edit(epoints):
        for feat in epoints.getFeatures():
            values = [eval(i) for i in feat['TargetID']]
            values.append(feat[weight_field])
            distances = [eval(i) for i in feat['Distance']]
            distances.append(0.0)
            #(value,distance,gausssian weight) in below
            values_filt = [(values[c],d,1/(np.sqrt(2*np.pi)*b) * np.exp(-0.5*(d/b)**2)) for c,d in enumerate(distances) if d<=b and type(values[c]) in (float,int)]
            
            s_conf = np.sqrt(((n*np.sum([v[1]**2 for v in values_filt]))-np.sum([v[1] for v in values_filt])**2) / (n-1))
            z_score = (np.sum([v[1]*v[0] for v in values_filt])- np.sum([v[1]*x_mean for v in values_filt])) / (s * s_conf)
            z_score = float(np.round(z_score,3))
            
            if type(z_score) in (float,int):
                feat['getisord']=z_score
            else:
                0
            epoints.updateFeature(feat)
    
    del_attributes = ['TargetID','Distance']
    idx = [epoints.fields().indexFromName(i) for i in del_attributes]
    epoints.dataProvider().deleteAttributes(idx)
    epoints.updateFields()

    return epoints

def biodiversity(input_points,id_field,species_field,bandwidth):
    alg_params = {
        'INPUT': input_points,
        'INPUT_FIELD': id_field,
        'MATRIX_TYPE': 0,  # Linear (N*k x 3) distance matrix
        'NEAREST_POINTS': 40,
        'TARGET': input_points,
        'TARGET_FIELD': species_field,
        'OUTPUT': 'TEMPORARY_OUTPUT'
    }
    dmatrix = processing.run('qgis:distancematrix', alg_params)['OUTPUT']

    # Aggregate
    alg_params = {
        'AGGREGATES': [{'aggregate': 'first_value','delimiter': ',','input': '"InputID"','length': 0,'name': 'InputID','precision': 0,'sub_type': 0,'type': 4,'type_name': 'int8'},
                        {'aggregate': 'array_agg','delimiter': ',','input': '"TargetID"','length': 0,'name': 'TargetID','precision': 0,'sub_type': 10,'type': 11,'type_name': 'stringlist'},
                        {'aggregate': 'array_agg','delimiter': ',','input': '"Distance"','length': 0,'name': 'Distance','precision': 0,'sub_type': 10,'type': 11,'type_name': 'stringlist'}],
        'GROUP_BY': '"InputID"',
        'INPUT': dmatrix,
        'OUTPUT':'TEMPORARY_OUTPUT'
    }
    agg = processing.run('native:aggregate', alg_params)['OUTPUT']

    # Join attributes by field value
    alg_params = {
        'DISCARD_NONMATCHING': False,
        'FIELD': id_field,
        'FIELDS_TO_COPY': [''],
        'FIELD_2': 'InputID',
        'INPUT': input_points,
        'INPUT_2': agg,
        'METHOD': 1,  # Take attributes of the first matching feature only (one-to-one)
        'PREFIX': '',
        'OUTPUT': 'TEMPORARY_OUTPUT'
    }
    epoints = processing.run('native:joinattributestable', alg_params)['OUTPUT']
    #calculating kernel_density with weights
    #points = QgsProcessingUtils.mapLayerFromString(outputs['JoinAttributesByFieldValue']['OUTPUT'],context)
    epoints.dataProvider().addAttributes([QgsField('biod',QVariant.Double)])
    epoints.updateFields()
    proportions = lambda s_list : [s/np.sum(s_list) for s in s_list]
    b = bandwidth
    
    
    with edit(epoints):
        for feat in epoints.getFeatures():
            values = [eval(i) for i in feat['TargetID']]
            values.append(feat[species_field])
            distances = [eval(i) for i in feat['Distance']]
            distances.append(0.0)
            #(value,distance,gausssian weight) in below
            #g_score = (num_points - 1) * ((local_means - global_mean) ** 2) / ((num_points * local_variances) * ((num_points - 1) / (num_points - 2)))
            values_filt = [values[c] for c,d in enumerate(distances) if d<=b]
            
            species_amount = pd.Series(values_filt).explode().value_counts().tolist()
            prop = proportions(species_amount)
            simpsons_index = 0
            for p in prop:
                simpsons_index += p ** 2
            
            simpsons = (1 - simpsons_index) / (1-1/len(species_amount)) if simpsons_index > 0 and len(species_amount)>1 else 0
            simpsons = float(np.round(simpsons,3))
            #feedback.pushInfo(str(simpsons))
            if type(simpsons) in (int,float):
                feat['biod']=simpsons
            else:
                0
            epoints.updateFeature(feat)
        
    del_attributes = ['TargetID','Distance']
    idx = [epoints.fields().indexFromName(i) for i in del_attributes]
    epoints.dataProvider().deleteAttributes(idx)
    epoints.updateFields()

    return epoints

def distance2location(input_layer,dest_layer):
    # Shortest line between features
    alg_params = {
        'DESTINATION': dest_layer,
        'DISTANCE': None,
        'METHOD': 0,  # Distance to Nearest Point on feature
        'NEIGHBORS': 1,
        'SOURCE': input_layer,
        'OUTPUT':"TEMPORARY_OUTPUT"
    }
    slines = processing.run('native:shortestline', alg_params)['OUTPUT']


    # Join attributes by field value
    alg_params = {
        'DISCARD_NONMATCHING': False,
        'FIELD': 'fid',
        'FIELDS_TO_COPY': ['distance'],
        'FIELD_2': 'fid',
        'INPUT': input_layer,
        'INPUT_2': slines,
        'METHOD': 1,  # Take attributes of the first matching feature only (one-to-one)
        'PREFIX': '',
        'OUTPUT': "TEMPORARY_OUTPUT"
    }
    
    distance_points = processing.run('native:joinattributestable', alg_params)['OUTPUT']
    
    for field in distance_points.fields():
        if field.name() == 'distance':

            with edit(distance_points):
                idx = distance_points.fields().indexFromName(field.name())
                distance_points.renameAttribute(idx, 'wdistance')
                distance_points.updateFields()

    return distance_points