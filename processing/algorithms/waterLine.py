from osgeo import gdal,gdal_array
import os
import numpy as np
import pandas as pd
from math import sqrt,pow
from qgis.core import *
from qgis.PyQt.QtCore import QVariant
import processing
import tempfile
from .hydrologytools import calculate_massflux_dinf
import requests

def getWater(input_polygon:QgsVectorLayer,taso):
    tempd = tempfile.TemporaryFile()
    tempd = tempd.name+str(taso)+'.tif'
    bbox = getBboxWmsFormat(input_polygon)
    
    ss = bbox[0].split(',')
    
    wmsurl = 'https://aineistot.metsakeskus.fi/metsakeskus/rest/services/Vesiensuojelu/'+taso+'/ImageServer/exportImage?'
    params = {"bbox":str(round(int(ss[0])-100,-1))+","+str(round(int(ss[1])-100,-1))+","+str(round(int(ss[2])+100,-1))+","+str(round(int(ss[3])+100,-1)),
                "bboxSR":3067,
                "size":str((round(int(ss[2])+100,-1)-round(int(ss[0])-100,-1))/2)+","+str((round(int(ss[3])+100,-1) - round(int(ss[1])-100,-1))/2),
                "imageSR":3067,
                "format":'tiff',
                "pixelType":"F32",
                "noData":-9999,
                "noDataInterpretation":"esriNoDataMatchAny",
                "interpolation":"+RSP_BilinearInterpolation",
                "f":"image"}

    try:
        respo= requests.get(wmsurl,params,allow_redirects=True)
        
        if respo.status_code != 200:
           info = "Cannot connect to "+str(taso)+ " data: "+str(wmsurl)
           infolevel = 3
        else:
            open(tempd,'wb').write(respo.content)
    
    except:
        info = "Cannot connect to "+str(taso)+ " data: "+str(wmsurl)
        infolevel = 3
    
        
    try:
        test = gdal.Open(tempd)
        test_b = test.GetRasterBand(1)
        test_a = test_b.ReadAsArray()
        if np.max(test_a) > 1:
            info = str(taso)+" data is ok!"
            infolevel = 1
            del test,test_b,test_a
        else:
           info = "Not able find "+str(taso)+" data from area: "+str(bbox[0])
           infolevel = 3
    except:
        info = "Not able find "+str(taso)+" data from area: "+str(bbox[0])
        infolevel = 3

    return tempd,info,infolevel

def feature2layer(feature):
    vl = QgsVectorLayer("Polygon", "temporary_pol","memory")
    pr = vl.dataProvider()
    # Enter editing mode
    vl.startEditing()
    pr.addAttributes( [QgsField("id",  QVariant.Int)])

    fet = QgsFeature()
    fet.setGeometry( feature.geometry())
    fet.setAttributes(feature.attributes())
    pr.addFeatures( [ fet ] )
    
    vl.commitChanges()
    vl.updateExtents()
    
    return vl

def clipRaster(in_raster,band,clip_raster,band_clip):
    """clip raster by other raster"""
    output = tempfile.TemporaryFile(prefix='riparian_',suffix='.tif')
    output = output.name
    
    
    in_arr = raster2Array(in_raster,band)
    cl_arr = raster2Array(clip_raster,band_clip)
    
    in_arr = np.where((cl_arr>0) & (in_arr>0),in_arr,0)
    
    gdal_array.SaveArray(in_arr.astype("float32"),output,"GTiff",in_raster)
    
    return output

def raster2Array(in_raster,band):
    rast = gdal.Open(in_raster)
    rastB = rast.GetRasterBand(band)
    rastA = rastB.ReadAsArray()

    return rastA

def array2raster(in_array,map_raster):
    output = tempfile.TemporaryFile(prefix='riparian_',suffix='.tif')
    output = output.name
    #tempd = tempfile.gettempdir()
    #tempdname = tempfile.gettempprefix()
    #tempd = os.path.join(tempd,tempdname+'.tif')
    
    gdal_array.SaveArray(in_array.astype("float32"),output,"GTiff",map_raster)

    return output

def cleanGeom(vector):
    """cleans input vector geometry by buffering algorithm"""
    #vector = QgsVectorLayer(vector,"vect","ogr")

    with edit(vector):
        for feat in vector.getFeatures():
            geom = feat.geometry()
            buffer = geom.buffer(10,5)
            buffer = buffer.buffer(-9,5)
            feat.setGeometry(buffer)

            vector.updateFeature(feat)

def calcFocal(in_array,dist):
    """calculates focal statistic retangle maximum from raster array """
    dat =pd.DataFrame(in_array)
    vert = dat
    ijlist = []
    for i in range(0-dist,dist):
        for j in range(0-dist,dist):
            e = sqrt(pow(i,2)+pow(j,2))
            if e <=dist:
                ijlist.append((i,j))
    
    for i in ijlist:
        df = dat.shift(i[0],axis=0)
        df = df.shift(i[1],axis=1)
        vert = np.maximum(df,vert)
    t = []
    t.append(vert)
    t = np.array(t)

    """return focal array of input raster array"""
    return t

def processRaster(input):
    """calculates areas of where max(raster) - raster == 2"""
    output = tempfile.TemporaryFile(prefix='riparian_',suffix='.tif')
    output = output.name
    #output= os.path.join(tempd,tempdname+'.tif')
    
    rast = gdal.Open(input)
    
    rastb = rast.GetRasterBand(1)
    rastA = rastb.ReadAsArray()
    rastA = np.where(rastA>0,rastA,0)
    
    focal = calcFocal(rastA,2)
    huip = np.where(focal-rastA==2,1,0)
    gdal_array.SaveArray(huip.astype("float32"),output,"GTiff",rast)
    
    return output


def snap2water(waterraster,area):
    """snap the input vector to water raster"""

    waterarr = raster2Array(waterraster,1)
    waterarr = np.where(waterarr>1,1,0)

    water = array2raster(waterarr,waterraster)

    data = {'cutarea':[1]}
    dat = pd.DataFrame(data)
    water = raster2vector(water,dat)
    
    snapped = processing.run("native:snapgeometries",
                             {'INPUT':area,
                              'REFERENCE_LAYER':water,'TOLERANCE':10,
                              'BEHAVIOR':1,
                              'OUTPUT':'TEMPORARY_OUTPUT'})
    
    cleanGeom(snapped['OUTPUT'])

    return snapped['OUTPUT']

def rasterizeVector(in_layer,gdal_extent,cells):
    """transform vector to raster"""
    output = tempfile.TemporaryFile(prefix='riparian_',suffix='.tif')
    output = output.name


    processing.run("gdal:rasterize",
                        {'INPUT':in_layer,
                        'FIELD':'',
                        'BURN':1,
                        'USE_Z':False,
                        'UNITS':1,
                        'WIDTH':cells,
                        'HEIGHT':cells,
                        'EXTENT':gdal_extent,
                        'NODATA':0,
                        'OPTIONS':'',
                        'DATA_TYPE':5,
                        'INIT':None,
                        'INVERT':False,
                        'EXTRA':'',
                        'OUTPUT':output})

    return output


def getWaterline(rasters,leimikko):
    """give waterline raster of input area if the area is within 10 meter to waterbody"""

    bbox = getBboxWmsFormat(leimikko)
    ss = bbox[0].split(',')
    extent = str(round(int(ss[0])-100,-1))+","+str(round(int(ss[2])+100,-1))+","+str(round(int(ss[1])-100,-1))+","+str(round(int(ss[3])+100,-1))+" ["+str(bbox[1])+"]"
    
    leimikko = snap2water(rasters[2],leimikko)
    vraster = processRaster(rasters[2]) #vesistörajan määritys

    leimraster = rasterizeVector(leimikko,extent,2)
    vrast_clip = clipRaster(vraster,1,leimraster,1) # rajataan leimikkoon
    #leimraster = clipRaster()
    leimraster = clipRaster(leimraster,1,rasters[0],2)

    return vrast_clip,leimraster

def getBboxWmsFormat(in_feat:QgsVectorLayer):
    desc=in_feat.extent()
    x_min=int(desc.xMinimum())
    y_min=int(desc.yMinimum())
    x_max=int(desc.xMaximum())+1
    y_max=int(desc.yMaximum())+1
    srid=str(in_feat.crs().authid())
    exte = str(x_min)+","+str(y_min)+","+str(x_max)+","+str(y_max)
    witdth = x_max - x_min
    height = y_max - y_min
    
    return exte,srid,witdth,height

def raster2vector(in_rast,data):
    """transfrom input buffer zone raster to vector"""
    vectn = processing.run("gdal:polygonize", 
        {'INPUT':in_rast,
        'BAND':1,
        'FIELD':'DN',
        'EIGHT_CONNECTEDNESS':False,
        'EXTRA':'',
        'OUTPUT':"TEMPORARY_OUTPUT"})
    
    vect = QgsVectorLayer(vectn['OUTPUT'],"vyohyke","ogr")
    arealist = [feat.geometry().area() for feat in vect.getFeatures() if feat['DN']==1]
    namelist = list(data.columns)
    for i in namelist:
        vect.dataProvider().addAttributes([QgsField(i,QVariant.Double)])
        vect.updateFields()
    
    with edit(vect):
        for feat in vect.getFeatures():
            #raster value 0 means out
            if feat['DN'] == 0:
                vect.deleteFeature(feat.id())
            
            #delete small parts
            if max(arealist) - feat.geometry().area() > max(arealist) /1.4:
                vect.deleteFeature(feat.id())

            for i in namelist:
                datac = data[[i]]
                #print (datac.iloc[0,0])
                feat[i] = float(datac.iloc[0,0])
            
            geom = feat.geometry()
            buffer = geom.buffer(10, 5)
            buffer = buffer.buffer(-10,5)
            feat.setGeometry(buffer)

            vect.updateFeature(feat)
    
    return vect


def fillSink(elev):
    """repair input elevation raster that flows goes to lowest areas"""
    tempd = tempfile.TemporaryFile()
    tempd = tempd.name+'.tif'
    
    processing.run("wbt:FillDepressions", 
                   {'dem':elev,
                    'fix_flats':True,
                    'flat_increment':None,
                    'max_depth':None,
                    'output':tempd})

    """
    processing.run("saga:fillsinkswangliu",
                   {'ELEV':elev,
                    'FILLED':tempd,
                    'FDIR':'TEMPORARY_OUTPUT',
                    'WSHED':'TEMPORARY_OUTPUT',
                    'MINSLOPE':0.1})
    """
    return tempd

def calcMassFlux(elev,rusle,ls,water):
    """calculates massflux of input parameters. This will be change
    to saga-gis flowaccumulation massflux algorithm when published"""
    output = tempfile.TemporaryFile()
    output = output.name+'.tif'
    
    mf = calculate_massflux_dinf(elev,rusle,ls,None)

    gdal_array.SaveArray(mf.astype("float32"),output,"GTiff",elev)
    
    """
    processing.run("wbt:DInfMassFlux",
                   {'dem':elev,
                    'loading':rusle,
                    'efficiency':ls,
                    'absorption':water,
                    'output':tempd})

    
    mf = processing.run("saga:flowaccumulationrecursive", 
                   {'ELEVATION':elev,
                    'SINKROUTE':None,
                    'WEIGHTS':ls,
                    'FLOW':'TEMPORARY_OUTPUT',
                    'VAL_INPUT':None,
                    'VAL_MEAN':'TEMPORARY_OUTPUT',
                    'ACCU_MATERIAL':rusle,
                    'ACCU_TARGET':water,
                    'ACCU_TOTAL':tempd,
                    'ACCU_LEFT':'TEMPORARY_OUTPUT',
                    'ACCU_RIGHT':'TEMPORARY_OUTPUT',
                    'FLOW_UNIT':1,
                    'TARGETS':None,
                    'FLOW_LENGTH':'TEMPORARY_OUTPUT',
                    'WEIGHT_LOSS':'TEMPORARY_OUTPUT',
                    'METHOD':2,'CONVERGENCE':1.1,
                    'NO_NEGATIVES':False})"""

    return output
    

def getMassSum(mf,waterborder):
    """calculates massflux sum of the input waterline raster"""
    mf_array = raster2Array(mf,1)
    water_arr = raster2Array(waterborder,1)
    
    mf_array = np.where((water_arr>0) & (mf_array>0),mf_array,0)
    massSum = np.sum(mf_array)
    
    return massSum


def getEffect(mf,mfmin,mfmax):
    """"calculates water protection attributes"""
    ret_max = (mfmax - mfmin) / 1000
    added_material = (mf - mfmin) / 1000
    reserved_material = (mfmax - mf) / 1000
    cost = round((reserved_material / ret_max) * 100,1)
    
    return ret_max,added_material,reserved_material,cost

def getBufferzone(rasters,clipraster,waterborder,dist,target):
    """Increase buffer zone area when distance and target are satisfied"""
    #change the rasters to numpy array
    demfill = rasters[3]

    # demfill = fillSink(rasters[3])
    zraster = raster2Array(rasters[0],1)
    eucarr = raster2Array(rasters[0],2)
    cuttarr = raster2Array(clipraster,1)
    lsarr = raster2Array(rasters[0],3)

    #change necessary value units and filter to clip raster
    zzone = np.where(cuttarr==1,zraster,0)
    z = zzone[zzone>0]
    
    lsarr = np.where(lsarr>0,lsarr / 100.0,0) #ls facto only to cliparea
    ls_max = np.where(cuttarr==1,1,lsarr)
    ls_min = lsarr #min and max scenarios
    
    #print ("ls_max: "+str(np.sum(ls_max)))

    rusarr = raster2Array(rasters[1],1)
    rusarr = np.where(rusarr>0,rusarr/10000*4,0.01)
    rus = array2raster(rusarr,rasters[3])

    ls_max = array2raster(ls_max,rasters[3])
    ls_min = array2raster(ls_min,rasters[3])

    
    
    mfmin = calcMassFlux(demfill,rus,ls_min,rasters[2])
    mfmaxr = calcMassFlux(demfill,rus,ls_max,rasters[2])
    mfmin = round(getMassSum(mfmin,waterborder),2)
    mfmax = round(getMassSum(mfmaxr,waterborder),2)
    #shutil.
    
    
    for i in range(0,100,5):
        #wraster = waterborder
        zp = np.percentile(z,i)
        ls_fact = np.where((cuttarr==1) & (zraster>zp) & (eucarr>=dist[0]),1,lsarr)
        #ls_fact = np.where((cuttarr==1) & (ls_fact<1) & (eucarr>=dist[2]),1,ls_fact)
    
        mdist = np.where((cuttarr==1) & (ls_fact<1),eucarr,0)
        mdist = mdist[mdist>0]
        mdist = round(np.mean(mdist)*2,1)
        #rus = np.where(rusarr>0,rusarr/10000*4*ls_fact,0.01)
        #rus = array2raster(rus,rasters[4])
        ls = array2raster(ls_fact,rasters[3])
        #print (mdist)
        #showRaster(ls)
        
        mfr = calcMassFlux(demfill,rus,ls,rasters[2])
        mf = round(getMassSum(mfr,waterborder),2)
        
        effect = getEffect(mf,mfmin,mfmax)
        if (mdist >= dist[1] and target[0] == False):
            for j in range(i-4,i+1,1):
                zp = np.percentile(z,j)
                ls_fact = np.where((cuttarr==1) & (zraster>zp) & (eucarr>=dist[0]),1,lsarr)
                #ls_fact = np.where((cuttarr==1) & (ls_fact<1) & (eucarr>=dist[2]),1,ls_fact)
                print ("ls: "+str(np.sum(ls_fact)))
                mdist = np.where((cuttarr==1) & (ls_fact<1),eucarr,0)
                mdist = mdist[mdist>0]
                mdist = round(np.mean(mdist)*2,1)
                #rus = np.where(rusarr>0,rusarr/10000*4*ls_fact,0.01)
                #rus = array2raster(rus,rasters[4])
                ls = array2raster(ls_fact,rasters[3])
                mfr = calcMassFlux(demfill,rus,ls,rasters[2])
                mf = round(getMassSum(mfr,waterborder),2)
                effect = getEffect(mf,mfmin,mfmax)
                if (mdist >= dist[1] and target[0] == False):
                    break
            
            break
            
        elif (target[0]==True and mdist>=dist[1] and effect[3]>=target[1]):
            for j in range(i-4,i+1,1):
                zp = np.percentile(z,j)
                ls_fact = np.where((cuttarr==1) & (zraster>zp) & (eucarr>=dist[0]),1,lsarr)
                #ls_fact = np.where((cuttarr==1) & (ls_fact<1) & (eucarr>=dist[2]),1,ls_fact)
    
                mdist = np.where((cuttarr==1) & (ls_fact<1),eucarr,0)
                mdist = mdist[mdist>0]
                mdist = round(np.mean(mdist)*2,1)
                #rus = np.where(rusarr>0,rusarr/10000*4*ls_fact,0.01)
                #rus = array2raster(rus,rasters[4])
                ls = array2raster(ls_fact,rasters[3])
                mfr = calcMassFlux(demfill,rus,ls,rasters[2])
                mf = round(getMassSum(mfr,waterborder),2)
                effect = getEffect(mf,mfmin,mfmax)
                if (mdist >= dist[1] and target[0] == False):
                    break
            break

    bzone = np.where((cuttarr==1) & (ls_fact<1),1,0)
    bzone = array2raster(bzone,rasters[1])

    dataset = {'massatase_max':[mfmax],
                'massatase':[mf],
                'pidatetty':[effect[3]],
                'keskileveys':[mdist]}
    
    df = pd.DataFrame(dataset)
    res = raster2vector(bzone,df)
    #mf = array2raster(mf,rasters[2])
    #mfmax= array2raster(mfmax,rasters[2])
    """return buffer zone as vector layer with dataset attributes"""
    return res