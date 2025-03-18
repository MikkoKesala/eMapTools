# -*- coding: utf-8 -*-

"""
This is part of eMapTool plugin. 

This one calculates diversity index from species field within search radius

"""

__author__ = 'Mikko Kesälä'
__date__ = '2024-08-14'
__copyright__ = '(C) 2024 by eMap modeling'

from qgis.core import QgsProcessing
from qgis.core import QgsProcessingAlgorithm
from qgis.core import QgsProcessingMultiStepFeedback
from qgis.core import QgsProcessingParameterFeatureSource
from qgis.core import QgsProcessingParameterFeatureSink
from qgis.core import QgsProcessingParameterField
from qgis.core import QgsProcessingParameterNumber
from qgis.core import QgsProcessingParameterEnum
from qgis.core import QgsProcessingUtils,QgsFeature,QgsGeometry,QgsFeatureSink,QgsField,QgsFields,QgsVectorLayer,QgsWkbTypes
from PyQt5.QtCore import QVariant
from qgis.PyQt.QtCore import QCoreApplication
import geopandas as gpd
import os
from .algorithms.geotools2 import gpd2qgis,valuesByDistance,layer2gpd
from .algorithms.ecoindices import diversityIndices



class BiodiversityIndices(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource('points','Vector points',[QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterField('species', 'Species field', type=QgsProcessingParameterField.Numeric, parentLayerParameterName='points', allowMultiple=False))
        self.addParameter(QgsProcessingParameterNumber('sdistance', 'Search distance (m)', type=QgsProcessingParameterNumber.Integer, minValue=5, maxValue=50, defaultValue=30))
        self.addParameter(QgsProcessingParameterEnum('method', 'method', options=['simpson','shannon','simpson_evenness','pielou'], allowMultiple=True, usesStaticStrings=False, defaultValue=[]))
        self.addParameter(QgsProcessingParameterFeatureSink('output', 'Output', type=QgsProcessing.TypeVectorPoint, createByDefault=True, defaultValue=None))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(2, model_feedback)
        results = {}
        outputs = {}
        
        # detect diversity method
        methods_dict = {0:'simpson',
                1:'shannon',
                2:'simpson_evenness',
                3:'pielou'}
        #get data
        points = QgsProcessingUtils.mapLayerFromString(parameters['points'],context)
        #methods = QgsProcessingUtils.mapLayerFromString(parameters['method'],context)

        methods = [methods_dict[m] for m in parameters['method']]
        points_gdf = layer2gpd(points)
        
        points_gdf = diversityIndices(points_gdf,parameters['species'],parameters['sdistance'],methods)
        points_gdf = points_gdf.drop(["distances","list_"+parameters['species']],axis=1)
        layer = gpd2qgis(points_gdf)

        

        (sink, dest_id) = self.parameterAsSink(parameters, 'output', context,
                            layer.fields(), QgsWkbTypes.Point, layer.crs())
        
        for feature in layer.getFeatures():
            sink.addFeature(feature, QgsFeatureSink.FastInsert)
        
        results['output'] = dest_id
        
        return results

    def name(self):
        return 'Biodiversity'

    def displayName(self):
        return 'Biodiversity within a radius'

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
        return 'Ecological indices'
    
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)
    

    def shortHelpString(self):
        helpfile = open(os.path.dirname(__file__) + '/descriptions/lstatistics.html',encoding="utf-8")
        help = helpfile.read()
        return help


    def createInstance(self):
        return BiodiversityIndices()