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
 This script initializes the plugin, making it known to QGIS.
"""

__author__ = 'eMap modeling'
__date__ = '2024-03-16'
__copyright__ = '(C) 2024 by eMap modeling'


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load eMapTools class from file eMapTools.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .ReTreeT import eMapToolsPlugin
    return eMapToolsPlugin()