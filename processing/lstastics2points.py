# -*- coding: utf-8 -*-

"""
This is part of eMapTool plugin. 

This one detects singletrees from canopy height model and enriched attribute data from false-color ortophoto

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
from qgis.core import QgsProcessingUtils,QgsFeature,QgsGeometry,QgsFeatureSink,QgsField,QgsFields,QgsVectorLayer,QgsWkbTypes
from PyQt5.QtCore import QVariant
from qgis.PyQt.QtCore import QCoreApplication
import geopandas as gpd
import os
from .algorithms.geotools2 import gpd2qgis,valuesByDistance,layer2gpd
from .algorithms.geostats import getLstats



class Lstatistics2points(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource('trees','Trees',[QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterField('field', 'Field for statistic', type=QgsProcessingParameterField.Numeric, parentLayerParameterName='trees', allowMultiple=False))
        self.addParameter(QgsProcessingParameterNumber('sdistance', 'Search distance (m)', type=QgsProcessingParameterNumber.Integer, minValue=5, maxValue=50, defaultValue=30))
        self.addParameter(QgsProcessingParameterFeatureSink('output_trees', 'Output trees', type=QgsProcessing.TypeVectorPoint, createByDefault=True, defaultValue=None))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(2, model_feedback)
        results = {}
        outputs = {}

        #get data
        trees = QgsProcessingUtils.mapLayerFromString(parameters['trees'],context)
        
        trees_gdf = layer2gpd(trees)

        trees_gdf = valuesByDistance(trees_gdf,parameters['field'],parameters['sdistance'])
        
        trees_gdf['lstat']= trees_gdf['list_'+parameters['field']].apply(lambda x: getLstats(x))
        trees_gdf[parameters['field']+'_lcv'] =trees_gdf['lstat'].apply(lambda x: x[0])
        trees_gdf[parameters['field']+'_lskew'] =trees_gdf['lstat'].apply(lambda x: x[1])
        
        trees_gdf = trees_gdf.drop(['lstat','list_'+parameters['field'],'distances'],axis=1)
        
        layer = gpd2qgis(trees_gdf)

        (sink, dest_id) = self.parameterAsSink(parameters, 'output_trees', context,
                            layer.fields(), QgsWkbTypes.Point, layer.crs())
        
        for feature in layer.getFeatures():
            sink.addFeature(feature, QgsFeatureSink.FastInsert)
        
        results['output_trees'] = dest_id
        
        return results

    def name(self):
        return 'lstatistitics'

    def displayName(self):
        return 'Calculate L-statistics'

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
        return Lstatistics2points()