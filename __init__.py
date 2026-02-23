# -*- coding: utf-8 -*-
"""
/***************************************************************************
 IndicatrixMapper
                                 A QGIS plugin
The plugin generates Quasi-Tissot Indicatrices as dynamic ellipses to visualize
map projection distortions. Using numerical derivatives (Quasi Indicatrix method),
it calculates precise scale factors and rotation angles for any project CRS.
Unlike traditional finite caps, these point-based indicatrices provide an
analytically exact representation of local distortion, updating in real-time
as the project projection changes. Users can instantly analyze conformality,
equivalence, and angular deformation through both visual ellipses and
detailed attribute data.
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
