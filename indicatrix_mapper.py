# -*- coding: utf-8 -*-
"""
/***************************************************************************
 IndicatrixMapperDialog
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
        copyright            : (C) 2015
        email                : 
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
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from qgis.core import *
# Initialize Qt resources from file resources.py
from . import resources_rc
# Import the code for the dialog
from .indicatrix_mapper_dialog import IndicatrixMapperDialog, calculate_indicatrix_params
import os.path


class IndicatrixMapper:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'indicatrix_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = IndicatrixMapperDialog()

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr('&Indicatrix mapper')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar('IndicatrixMapper')
        self.toolbar.setObjectName('IndicatrixMapper')

        # Connect to CRS change signal
        QgsProject.instance().crsChanged.connect(self.refresh_all_indicatrix_layers)

    def refresh_all_indicatrix_layers(self):
        """Refresh distortion attributes for all indicatrix layers in the project."""
        target_crs = QgsProject.instance().crs()
        
        for layer in QgsProject.instance().mapLayers().values():
            if layer.customProperty("indicatrix") == "caps":
                self.refresh_layer_distortion(layer, target_crs)

    def refresh_layer_distortion(self, layer, target_crs):
        """Update h, k, s, omega for a specific indicatrix layer."""
        a = float(layer.customProperty("ellipsoid_a"))
        b = float(layer.customProperty("ellipsoid_b"))
        delta = float(layer.customProperty("numerical_delta", 0.0001))
        
        # Source CRS is always the geographic one the layer was created with
        # Replicate WKT logic to ensure name matches
        if a == b:
            inv_f = 0
        else:
            inv_f = a / (a - b)
            
        wkt = (
            'GEOGCRS["Indicatrix Source Ellipsoid", '
            'DATUM["Custom Datum", '
            'ELLIPSOID["Custom Spheroid", {a}, {inv_f}, LENGTHUNIT["metre", 1]]], '
            'PRIMEM["Greenwich", 0, ANGLEUNIT["degree", 0.0174532925199433]], '
            'CS[ellipsoidal, 2], '
            'AXIS["geodetic latitude (Lat)", north, ORDER[1], ANGLEUNIT["degree", 0.0174532925199433]], '
            'AXIS["geodetic longitude (Lon)", east, ORDER[2], ANGLEUNIT["degree", 0.0174532925199433]]]'
        ).format(a=a, inv_f=inv_f)
        
        source_crs = QgsCoordinateReferenceSystem()
        source_crs.createFromWkt(wkt)
        
        transform = QgsCoordinateTransform(source_crs, target_crs, QgsProject.instance())
        
        # Get field indices
        h_idx = layer.fields().indexOf("h")
        k_idx = layer.fields().indexOf("k")
        s_idx = layer.fields().indexOf("s")
        omega_idx = layer.fields().indexOf("omega")
        a_ind_idx = layer.fields().indexOf("a_ind")
        b_ind_idx = layer.fields().indexOf("b_ind")
        angle_idx = layer.fields().indexOf("angle")
        lon_idx = layer.fields().indexOf("lon")
        lat_idx = layer.fields().indexOf("lat")
        
        if h_idx == -1 or k_idx == -1:
            return

        layer.startEditing()
        for feature in layer.getFeatures():
            lon = feature.attribute(lon_idx)
            lat = feature.attribute(lat_idx)
            
            h, k, s, omega, a_ind, b_ind, angle = calculate_indicatrix_params(lon, lat, transform, a, b, delta)
            
            layer.changeAttributeValue(feature.id(), h_idx, h)
            layer.changeAttributeValue(feature.id(), k_idx, k)
            layer.changeAttributeValue(feature.id(), s_idx, s)
            layer.changeAttributeValue(feature.id(), omega_idx, omega)
            if a_ind_idx != -1:
                layer.changeAttributeValue(feature.id(), a_ind_idx, a_ind)
                layer.changeAttributeValue(feature.id(), b_ind_idx, b_ind)
                layer.changeAttributeValue(feature.id(), angle_idx, angle)
            
        layer.commitChanges()

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('IndicatrixMapper', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToVectorMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/IndicatrixMapper/icon.png'
        self.add_action(
            icon_path,
            text=self.tr('Indicatrix mapper'),
            callback=self.run,
            parent=self.iface.mainWindow())


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        # Disconnect the CRS change signal
        try:
            QgsProject.instance().crsChanged.disconnect(self.refresh_all_indicatrix_layers)
        except:
            pass

        for action in self.actions:
            self.iface.removePluginVectorMenu(
                self.tr('&Indicatrix mapper'),
                action)
            self.iface.removeToolBarIcon(action)


    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass
