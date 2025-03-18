# -*- coding: utf-8 -*-

"""
/***************************************************************************
 emaptools
                                 A QGIS plugin
        begin                : 2025-01-24
        copyright            : (C) 2025 by emaptools
        email                : mikko.kesala@helsinki.fi
 ***************************************************************************/

"""

__author__ = 'Mikko Kesälä'
__date__ = '2025-01-24'
__copyright__ = '(C) 2025 by emaptools'

__revision__ = '$Format:%H$'

import os,sys
import numpy as np
from osgeo import gdal
from qgis.PyQt.QtCore import QCoreApplication
from qgis.PyQt.QtGui import QIcon
from qgis.core import (
                       QgsProcessingUtils,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterRasterDestination)
from .algorithms.hydrologytools import d_inf2

class slope(QgsProcessingAlgorithm):
    """
    This is an algorithm that takes a vector layer and
    get background rasters via interface the vector area
    and creates a new buffer zone along the waterbody.
    Input layer should be near waterbody

    """


    def initAlgorithm(self, config):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
        self.addParameter(QgsProcessingParameterRasterLayer('dem', 'DEM', defaultValue=None))
        self.addParameter(QgsProcessingParameterRasterDestination('slope', 'Slope', createByDefault=True, defaultValue=None))

       
    def processAlgorithm(self, parameters, context, feedback):
        
        feedb = {1:feedback.pushInfo,
                2:feedback.pushWarning,
                3:feedback.reportError}


        feedback.pushInfo("Calculating flow accumulation")

        results = {}

        #get rasters
        dem = QgsProcessingUtils.mapLayerFromString(parameters['dem'],context)
        demname = dem.source()

        lsout = self.parameterAsOutputLayer(parameters,"slope",context)
        slope = d_inf2(demname)[1]
        
        ds = gdal.Open(demname)
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
        [rows, cols] = arr.shape
        driver = gdal.GetDriverByName("GTiff")
        outdata = driver.Create(lsout, cols, rows, 1, gdal.GDT_Float32)
        outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
        outdata.SetProjection(ds.GetProjection())##sets same projection as input
        outdata.GetRasterBand(1).WriteArray(slope)
        outdata.GetRasterBand(1).SetNoDataValue(-9999) ##if you want these values transparent
        outdata.FlushCache() ##saves to disk!!
        outdata = None
        ds = None
        band = None
        arr = None

        results['slope'] = lsout

        return results

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Slope (Tarbaton 1997)'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr(self.name())

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
        return 'Hydrology'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return slope()

