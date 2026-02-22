# -*- coding: utf-8 -*-
"""
/***************************************************************************
 IndicatrixMapper
                                 A QGIS plugin
The plugin introduces ellipsoidal or spherical caps which can give a
down-to-earth Tissot-indicatrix realization. The Tissot-indicatrix uses constant
radius ellipsoidal caps (calculated by Vincenty's formula) instead of the original
infinitesimal small Tissot circles. These caps are able to transform on the fly
from a reference ellipsoid to a selected project coordinate reference system.
The user can study the distortions of the caps in a blink, such conclusions
can be drawn as the projection is conformal or equal-area.
                             -------------------
        begin                : 2015-03-15
        git sha              : $Format:%H$
        copyright            : (C) 2015 by Ervin Wirth
        email                : wirth.ervin@mailbox.org
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


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load IndicatrixMapper class from file indicatrix_mapper.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .indicatrix_mapper import IndicatrixMapper
    return IndicatrixMapper(iface)
