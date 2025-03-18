# -*- coding: utf-8 -*-

"""
This is part of eMapTool plugin. 

This one plans retention tree areas within planning area which is polygon input of the tool.
Retention tree planning is based on ecological values and weighing of the values.

"""

__author__ = 'Mikko Kesälä'
__date__ = '2024-01-01'
__copyright__ = '(C) 2024 by eMap modeling'

# This will get replaced with a git SHA1 when you do a git archive


from stat import S_ISLNK
from qgis import processing
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtCore import QCoreApplication,QVariant
from qgis.core import (QgsProcessing,
                       QgsField,
                       QgsFeatureSink,QgsProcessingParameterField,
                       QgsProcessingParameterFeatureSource,QgsProcessingParameterRasterDestination,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterMapLayer,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterDefinition,
                       QgsProcessingUtils,QgsRasterLayer,QgsVectorLayer)
import os,time,sys,tempfile

#from PIL import Image
from .algorithms.geotools import feature2Layer,createTreeMap,addFieldValue,joinIntersection,point2area,clipRaster3,distance2location
from .algorithms.ecomodels import runEssModel
from .algorithms.geotools2 import clipByFeature,clipBigVRTrasterByFeature,layer2gpd


class planReTreeAreas(QgsProcessingAlgorithm):


    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    OUTPUT = 'OUTPUT'
    AREAS = 'AREAS'
    FOSFORI = 'FOSFORI'
    DTW1 = 'DTW1'
    BIOD = 'BIOD'
    LAHOP = 'LAHOP'
    PUUM = 'PUUM'
    INPUT = 'INPUT'

    
    delfields = ['fid','OBJECTID','layer','path','leimikko','SPECIALFEATURECODE', 'SPECIALFEATUREADDITIONALCODE', 'DEVELOPMENTCLASS', 'STEMCOUNTPINE', 'STEMCOUNTDECIDUOUS', 'STEMCOUNTSPRUCE', 'PaajakoNro', 'Nimi_2','MEANDIAMETERDECIDUOUS', 'MEANDIAMETERPINE', 'MEANDIAMETERSPRUCE', 'MEANHEIGHTDECIDUOUS', 'MEANHEIGHTPINE', 'MEANHEIGHTSPRUCE']
    
    grid_fields = 'GEOMETRY,FERTILITYCLASS,MEANDIAMETERDECIDUOUS,MEANDIAMETERPINE,MEANDIAMETERSPRUCE,MEANHEIGHTDECIDUOUS,MEANHEIGHTPINE,MEANHEIGHTSPRUCE'
    mkasviv_fields = 'PaajakoNro,Nimi'
 
    def initAlgorithm(self, config):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        #inputs
        self.addParameter(QgsProcessingParameterFeatureSource(self.INPUT,self.tr('Planning area'),[QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterField('id_field', 'ID field', type=QgsProcessingParameterField.Any, parentLayerParameterName=self.INPUT, allowMultiple=False))
        self.addParameter(QgsProcessingParameterMapLayer('chm', 'CHM', types=[QgsProcessing.TypeRaster]))
        self.addParameter(QgsProcessingParameterMapLayer('dtw', 'DTW',types=[QgsProcessing.TypeRaster]))
        self.addParameter(QgsProcessingParameterMapLayer('vegetationzone', 'Vegetation zones', types=[QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterMapLayer('waterbody', 'Water body', types=[QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterMapLayer('forestgrid', 'Forest grid',types=[QgsProcessing.TypeVectorPolygon]))
 
        
        #parameters
        params = []
        params.append(QgsProcessingParameterEnum(self.FOSFORI,self.tr('Phosphorus retention (Pret)'),options=['No weighting','Low','Moderate','Significant'],defaultValue=1))
        params.append(QgsProcessingParameterEnum(self.DTW1,self.tr('Depth to water -index (DTW)'),options=['No weighting','Low','Moderate','Significant'],defaultValue=1))
        params.append(QgsProcessingParameterEnum(self.BIOD,self.tr('Tree diversity (Ds)'),options=['No weighting','Low','Moderate','Significant'],defaultValue=1))
        params.append(QgsProcessingParameterEnum(self.LAHOP,self.tr('Dead Wood Potential (DWP)'),options=['No weighting','Low','Moderate','Significant'],defaultValue=1))
        
        for p in params:
            p.setFlags(p.flags() | QgsProcessingParameterDefinition.FlagAdvanced) 
            self.addParameter(p)

        self.addParameter(QgsProcessingParameterNumber(self.PUUM,self.tr('Retention tree count (pcs. / ha)'),type=QgsProcessingParameterNumber.Integer,minValue=5,maxValue=30,defaultValue=10))
        #outputs
        self.addParameter(QgsProcessingParameterFeatureSink(self.OUTPUT,self.tr('Trees')))
        self.addParameter(QgsProcessingParameterFeatureSink(self.AREAS,self.tr('Retention tree groups')))

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        results = {}
        feedb = {1:feedback.setProgressText,
                2:feedback.pushWarning,
                3:feedback.reportError}
        # Retrieve the feature source and sink. The 'dest_id' variable is used
        # to uniquely identify the feature sink, and must be included in the
        # dictionary returned by the processAlgorithm function.
        source = self.parameterAsSource(parameters, self.INPUT, context)
        if source.featureCount() > 50:
            feedback.reportError("Input layer has too many features. 20 features is maximum. Process failed.")
            sys.exit()

        #feedback.pushInfo(str(parameters['INPUT']))
        source = processing.run("native:dissolve", {'INPUT':parameters['INPUT'],'FIELD':[str(parameters['id_field'])],'OUTPUT':'TEMPORARY_OUTPUT'})
        source = processing.run("native:multiparttosingleparts", {'INPUT':source['OUTPUT'],'OUTPUT':'TEMPORARY_OUTPUT'})
        source = source['OUTPUT']
       
        
        # Compute the number of steps to display within the progress bar and
        #slayer = QgsVectorLayer(parameters['INPUT'],"slayer","ogr")
        # get features from source
        
        total = 100.0 / source.featureCount() if source.featureCount() else 0
        
        features = source.getFeatures()
        for current, feature in enumerate(features):
            #feedback.setProgressText("testia ja "+str(out))
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break
            
            #out.select(feature.id())
            #uri = "polygon?crs="+str(source.sourceCrs())+"&field=id:integer&"
            #QgsVectorLayer(uri,"mem","memory")
            #feedback.pushInfo(parameters['dtw'].valueAsPythonString())
            feedback.setProgressText("Site: "+str(feature[parameters['id_field']]))
            feedback.setProgressText("Cutting data by planning area")

            #leimFeat = feature.getFeatures()
            leimArea = [feature.geometry().area()/10000]
            out = feature2Layer(feature,100)
            out1 = feature2Layer(feature,0)
            out.setCrs(source.crs())
            out1.setCrs(source.crs())
            leim = addFieldValue(out1,"leimikko",1)

            #out_gpd = layer2gpd(out)
            
            try:
                chm = QgsProcessingUtils.mapLayerFromString(parameters['chm'],context)
                chm = clipRaster3(chm,out)
                
                dtw = QgsProcessingUtils.mapLayerFromString(parameters['dtw'],context)
                dtw = clipRaster3(dtw,out)
                
                #fgrid = QgsProcessingUtils.mapLayerFromString(parameters['forestgrid'],context)
                #biogeo = QgsProcessingUtils.mapLayerFromString(parameters['vegetaionzone'],context)


                feedback.setProgressText("Creating tree map using file")

                outChm = createTreeMap(chm,3,False)

                feedback.setProgressText("Spatial joins and raster sampling data to trees")
                
                if dtw:
                    outChm = processing.run("native:rastersampling", {'INPUT':outChm,'RASTERCOPY':dtw,'COLUMN_PREFIX':'DTW_','OUTPUT':'TEMPORARY_OUTPUT'})
                    outChm = outChm['OUTPUT']
                else:
                    outChm.dataProvider().addAttributes([QgsField("DTW_1",QVariant.Double)])
                    outChm.updateFields()
                
                
                fgrid = processing.run('native:clip',{'INPUT':parameters['forestgrid'],
                                                    'OVERLAY':out,
                                                    'OUTPUT':'TEMPORARY_OUTPUT'}, context=context, feedback=feedback, is_child_algorithm=False)
                
                fgrid = fgrid['OUTPUT']
                biogeo = processing.run('native:clip',{'INPUT':parameters['vegetationzone'],
                                                    'OVERLAY':out,
                                                    'OUTPUT':'TEMPORARY_OUTPUT'}, context=context, feedback=feedback, is_child_algorithm=False)
                
                biogeo = biogeo['OUTPUT']

                #outChm = joinIntersection(outChm,fgrid,list(self.grid_fields.split(",")),False)
                outChm = processing.run("native:joinattributesbylocation", {'INPUT':outChm,'JOIN':fgrid,'PREDICATE':[0],'JOIN_FIELDS':list(self.grid_fields.split(",")),'METHOD':0,'DISCARD_NONMATCHING':False,'PREFIX':'','OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']
                outChm = processing.run("native:joinattributesbylocation", {'INPUT':outChm,'JOIN':biogeo,'PREDICATE':[0],'JOIN_FIELDS':['paajakonro'],'METHOD':0,'DISCARD_NONMATCHING':False,'PREFIX':'','OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']
                outChm = processing.run("native:joinattributesbylocation", {'INPUT':outChm,'JOIN':leim,'PREDICATE':[0],'JOIN_FIELDS':['leimikko',str(parameters['id_field'])],'METHOD':0,'DISCARD_NONMATCHING':True,'PREFIX':'','OUTPUT':'TEMPORARY_OUTPUT'})['OUTPUT']
                
                #outChm = joinIntersection(outChm,biogeo,['paajakonro'],False)
                #outChm = joinIntersection(outChm,leim,['leimikko'],True)
                
                feedback.setProgressText("Calculating distance to waterbody")
                outChm = distance2location(outChm,parameters['waterbody'])
                feedback.setProgressText("Calculating ecosystem services values")

                fosf = self.parameterAsInt(parameters,self.FOSFORI,context)
                dtw1 = self.parameterAsInt(parameters,self.DTW1,context)
                biod = self.parameterAsInt(parameters,self.BIOD,context)
                lahop = self.parameterAsInt(parameters,self.LAHOP,context)
                weights ={"NP":float(fosf),"BIO":float(biod),"LP":float(lahop),"DTW":float(dtw1)}
                puuMaara = self.parameterAsInt(parameters,self.PUUM,context)

                out = runEssModel(outChm,weights,puuMaara,leimArea[0],"paajakonro")
                #feedback.pushInfo(str(out.featureCount()))
                #out = outChm
            except Exception as e:
                feedback.pushWarning(e)
   

            if current == 0:
                (sink, dest_id) = self.parameterAsSink(parameters, self.OUTPUT,context,
                    out.fields(), out.wkbType(), out.crs())
                
            #feedback.pushInfo(str(out.fields().names()))
            outFeats = out.getFeatures()
            for outFeat in outFeats:
                #feedback.pushInfo(str(outFeat['CHM']))
                sink.addFeature(outFeat, QgsFeatureSink.FastInsert)
        
        style = os.path.join(os.path.dirname(__file__),"styles/reTree4.qml")
        style2 = os.path.join(os.path.dirname(__file__),"styles/retreet_areas.qml")

        layer = QgsProcessingUtils.mapLayerFromString(dest_id, context)
        layer.loadNamedStyle(style)
        
        reareas = point2area(layer,'reTree',1)
        
        (sink, area_id) = self.parameterAsSink(parameters, self.AREAS,context,
                    reareas.fields(), reareas.wkbType(), reareas.crs())
        outFeats = reareas.getFeatures()
        
        for outFeat in outFeats:
                #feedback.pushInfo(str(outFeat['CHM']))
                sink.addFeature(outFeat, QgsFeatureSink.FastInsert)
        
        layer2 = QgsProcessingUtils.mapLayerFromString(area_id, context)
        layer2.loadNamedStyle(style2)

        return {self.OUTPUT: dest_id}


    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'planretreeareas'


    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return 'Plan retention tree groups'

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'ReTreeT planning'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)
    
    def shortHelpString(self):
        helpfile = open(os.path.dirname(__file__) + '/descriptions/planretreetareas.html',encoding="utf-8")
        help = helpfile.read()
        return help


    def createInstance(self):
        return planReTreeAreas()
