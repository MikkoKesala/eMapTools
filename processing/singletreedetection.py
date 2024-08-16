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
from qgis.core import QgsProcessingParameterRasterLayer
from qgis.core import QgsProcessingParameterFeatureSink
from qgis.core import QgsProcessingUtils,QgsFeature,QgsGeometry,QgsFeatureSink,QgsField,QgsFields,QgsVectorLayer,QgsWkbTypes
from PyQt5.QtCore import QVariant
from qgis.PyQt.QtCore import QCoreApplication
import geopandas as gpd
import os
from .algorithms.geotools2 import gpd2qgis,singleTreeMapping


class Singletree_base(QgsProcessingAlgorithm):

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer('chm', 'chm', defaultValue=None))
        self.addParameter(QgsProcessingParameterRasterLayer('cir', 'cir', defaultValue=None))
        self.addParameter(QgsProcessingParameterFeatureSink('trees', 'trees', type=QgsProcessing.TypeVectorPoint, createByDefault=True, defaultValue=None))

    def processAlgorithm(self, parameters, context, model_feedback):
        # Use a multi-step feedback, so that individual child algorithm progress reports are adjusted for the
        # overall progress through the model
        feedback = QgsProcessingMultiStepFeedback(2, model_feedback)
        results = {}
        outputs = {}

        #get rasters
        chm = QgsProcessingUtils.mapLayerFromString(parameters['chm'],context)
        chmname = chm.source()

        orto = QgsProcessingUtils.mapLayerFromString(parameters['cir'],context)
        ortoname = orto.source()
        trees = singleTreeMapping(chmname,ortoname,2)

        
        layer = gpd2qgis(trees)

        (sink, dest_id) = self.parameterAsSink(parameters, 'trees', context,
                            layer.fields(), QgsWkbTypes.Point, layer.crs())
        
        for feature in layer.getFeatures():
            sink.addFeature(feature, QgsFeatureSink.FastInsert)
        
        results['trees'] = dest_id
        
        return results

    def name(self):
        return 'singletrees'

    def displayName(self):
        return 'Detect single trees'

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
        return 'Tree mapping'
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)
    

    def shortHelpString(self):
        helpfile = open(os.path.dirname(__file__) + '/descriptions/singletreedetection.html',encoding="utf-8")
        help = helpfile.read()
        return help


    def createInstance(self):
        return Singletree_base()