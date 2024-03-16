# -*- coding: utf-8 -*-

"""
/***************************************************************************
 eMapTools
                                 A QGIS plugin
 This plugin propose retention trees and riparian buffer zones based on ecological values
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2024-03-16
        copyright            : (C) 2024 by eMap modeling
        email                : kesalamj@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'eMap modeling'
__date__ = '2024-03-16'
__copyright__ = '(C) 2024 by eMap modeling'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'
import os
from qgis.PyQt.QtGui import QIcon
from qgis.core import QgsProcessingProvider
from .processing.planretreetareas import planReTreeAreas
from .processing.points2retreeareas import points2retreeareas


class eMapToolsProvider(QgsProcessingProvider):

    def __init__(self):
        """
        Default constructor.
        """
        QgsProcessingProvider.__init__(self)

    def unload(self):
        """
        Unloads the provider. Any tear-down steps required by the provider
        should be implemented here.
        """
        pass

    def loadAlgorithms(self):
        """
        Loads all algorithms belonging to this provider.
        """
        self.addAlgorithm(planReTreeAreas())
        self.addAlgorithm(points2retreeareas())
        # add additional algorithms here
        # self.addAlgorithm(MyOtherAlgorithm())

    def id(self):
        """
        Returns the unique provider id, used for identifying the provider. This
        string should be a unique, short, character only string, eg "qgis" or
        "gdal". This string should not be localised.
        """
        return 'eMapTools'

    def name(self):
        """
        Returns the provider name, which is used to describe the provider
        within the GUI.

        This string should be short (e.g. "Lastools") and localised.
        """
        return self.tr('eMapTools')

    def icon(self):
        """
        add icon
        """
        pluginPath = os.path.dirname(__file__)
        iconPath = os.path.join(pluginPath, 'emap.png')

        return QIcon(os.path.join(iconPath))

    def longName(self):
        """
        Returns the a longer version of the provider name, which can include
        extra details such as version numbers. E.g. "Lastools LIDAR tools
        (version 2.2.1)". This string should be localised. The default
        implementation returns the same string as name().
        """
        return self.name()
