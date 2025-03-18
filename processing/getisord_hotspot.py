# -*- coding: utf-8 -*-

"""
This is part of eMapTool plugin. 

We calculate here Getis-Ord z-score value, which show hotspot areas on data

"""

__author__ = 'Mikko Kesälä'
__date__ = '2024-09-11'
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
from .algorithms.geostats import getisord


class GetisOrdHotspot(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource('vectorpoint','Vector point',[QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterField('field', 'Field for statistic', type=QgsProcessingParameterField.Numeric, parentLayerParameterName='vectorpoint', allowMultiple=False))
        self.addParameter(QgsProcessingParameterNumber('sdistance', 'Search distance (m)', type=QgsProcessingParameterNumber.Integer, minValue=5, maxValue=50, defaultValue=30))
        self.addParameter(QgsProcessingParameterEnum('method', 'method', options=['linear','quadric','gaussian','fixed'], allowMultiple=False, usesStaticStrings=False, defaultValue=[]))
        self.addParameter(QgsProcessingParameterFeatureSink('output', 'Output', type=QgsProcessing.TypeVectorPoint, createByDefault=True, defaultValue=None))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(2, model_feedback)
        results = {}
        outputs = {}

        #get data
        vpoints = QgsProcessingUtils.mapLayerFromString(parameters['vectorpoint'],context)
        
        # detect spatial relationship method
        methods = {0:'linear',
                1:'quadric',
                2:'gaussian',
                3:'fixed'}

        method = methods[parameters['method']]

        # converting to geopandas dataframe
        vpoints_gdf = layer2gpd(vpoints)

        # Calculating getisord value
        vpoints_gdf = getisord(vpoints_gdf,parameters['field'],parameters['sdistance'],method)
        
        layer = gpd2qgis(vpoints_gdf)

        (sink, dest_id) = self.parameterAsSink(parameters, 'output', context,
                            layer.fields(), QgsWkbTypes.Point, layer.crs())
        
        for feature in layer.getFeatures():
            sink.addFeature(feature, QgsFeatureSink.FastInsert)
        
        results['output'] = dest_id
        
        return results

    def name(self):
        return 'getisord'

    def displayName(self):
        return 'Getis-Ord hotspot analysis'

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
        return 'Geostatistics'
    
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)
    

    def shortHelpString(self):
        helpfile = open(os.path.dirname(__file__) + '/descriptions/lstatistics.html',encoding="utf-8")
        help = helpfile.read()
        return help


    def createInstance(self):
        return GetisOrdHotspot()