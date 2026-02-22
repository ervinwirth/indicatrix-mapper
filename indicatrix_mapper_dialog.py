# -*- coding: utf-8 -*-
"""
/***************************************************************************
 IndicatrixMapperDialog
                                 A QGIS plugin
The plugin introduces ellipsoidal or spherical caps which can give a
down-to-earth Tissot-indicatrix realization. These caps use constant
radius ellipsoidal caps (calculated by Vincenty's formula) instead of the original
infinitesimal small Tissot circles. These caps are able to transform on the fly
from a reference ellipsoid to a selected project coordinate reference system.
The user can study the distortions of the caps in a blink, such conclusions
can be drawn as the projection is conformal or equal-area.
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
from qgis.gui import *
from .indicatrix_mapper_dialog_base import Ui_IndicatrixMapper
import numpy
import qgis.utils

def calculate_indicatrix_params(lon, lat, transform, a, b, delta=0.0001):
    """Calculate distortion parameters using numerical derivatives.
    Based on the Quasi Indicatrix method by Bildirici, Ö. G. (2015).
    Uses Central Difference for O(delta^2) accuracy.
    """
    try:
        # Increment points for Central Difference
        p_lat_plus = QgsPointXY(lon, lat + delta)
        p_lat_minus = QgsPointXY(lon, lat - delta)
        p_lon_plus = QgsPointXY(lon + delta, lat)
        p_lon_minus = QgsPointXY(lon - delta, lat)
        
        # Project points
        pp_lat_p = transform.transform(p_lat_plus)
        pp_lat_m = transform.transform(p_lat_minus)
        pp_lon_p = transform.transform(p_lon_plus)
        pp_lon_m = transform.transform(p_lon_minus)
        
        # Ellipsoid parameters
        e2 = 1 - (b**2 / a**2)
        lat_rad = numpy.radians(lat)
        
        # Radii of curvature
        M = a * (1 - e2) / (1 - e2 * numpy.sin(lat_rad)**2)**1.5
        N = a / numpy.sqrt(1 - e2 * numpy.sin(lat_rad)**2)
        
        # Geodesic distances on ellipsoid for the 2*delta increments
        ds_meridian_ellips = numpy.radians(2 * delta) * M
        ds_parallel_ellips = numpy.radians(2 * delta) * N * numpy.cos(lat_rad)
        
        # Numerical derivatives using Central Difference
        dx_dphi = (pp_lat_p.x() - pp_lat_m.x()) / ds_meridian_ellips
        dy_dphi = (pp_lat_p.y() - pp_lat_m.y()) / ds_meridian_ellips
        dx_dlam = (pp_lon_p.x() - pp_lon_m.x()) / ds_parallel_ellips
        dy_dlam = (pp_lon_p.y() - pp_lon_m.y()) / ds_parallel_ellips
        
        # Scale factors h and k
        h = numpy.sqrt(dx_dphi**2 + dy_dphi**2)
        k = numpy.sqrt(dx_dlam**2 + dy_dlam**2)
        
        # Area distortion s = h * k * sin(theta')
        s = numpy.abs(dx_dphi * dy_dlam - dx_dlam * dy_dphi)
        
        # Angular distortion omega
        a_prime = numpy.sqrt(h**2 + k**2 + 2*s)
        b_prime = numpy.sqrt(max(0, h**2 + k**2 - 2*s))
        a_ind = (a_prime + b_prime) / 2
        b_ind = (a_prime - b_prime) / 2
        
        omega = 2 * numpy.degrees(numpy.arcsin((a_ind - b_ind) / (a_ind + b_ind)))
        
        return float(h), float(k), float(s), float(omega)
    except:
        return 1.0, 1.0, 1.0, 0.0


# FORM_CLASS, _ = uic.loadUiType(os.path.join(
# os.path.dirname(__file__), 'indicatrix_mapper_dialog_base.ui'))


class IndicatrixMapperDialog(QDialog, Ui_IndicatrixMapper):

    def __init__(self):
        """Constructor."""
        QDialog.__init__(self)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)
        self.ui = Ui_IndicatrixMapper()
        
        # Resize dialog to accommodate new fields and avoid overlap
        self.resize(605, 260)
        self.setMinimumSize(QtCore.QSize(600, 260))

        self.all_signals()
        self.min_lat_change()
        self.max_lat_change()
        self.min_lon_change()
        self.max_lon_change()

        # Inject Numerical Delta UI
        self.lblDelta = QLabel(self)
        self.lblDelta.setGeometry(QRect(405, 135, 70, 20))
        self.lblDelta.setText("Delta [deg]")
        self.lblDelta.setFont(self.label_12.font())
        self.lblDelta.setAlignment(Qt.AlignRight | Qt.AlignTrailing | Qt.AlignVCenter)

        self.spbDelta = QDoubleSpinBox(self)
        self.spbDelta.setGeometry(QRect(480, 135, 111, 22))
        self.spbDelta.setDecimals(9)
        self.spbDelta.setMinimum(0.000000001)
        self.spbDelta.setMaximum(1.0)
        self.spbDelta.setSingleStep(0.1)
        self.spbDelta.setProperty("value", 0.0001)
        self.spbDelta.setStepType(QAbstractSpinBox.AdaptiveDecimalStepType)
        self.spbDelta.setObjectName("spbDelta")

        # Move Auxiliary points checkbox further down to avoid overlap
        self.chkBoxAux.setGeometry(QRect(460, 175, 130, 20))
        
        # Move Run button further down
        self.btnRun.setGeometry(QRect(460, 215, 131, 23))

        self.inputs()

    def get_source_crs(self, a, b):
        """Create a geographic CRS with a custom name using WKT."""
        # Calculate inverse flattening
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
        
        crs = QgsCoordinateReferenceSystem()
        crs.createFromWkt(wkt)
        return crs

    def all_signals(self): 
        self.btnRun.clicked.connect(self.fun_execute)
        self.spbRadiusDg.valueChanged.connect(self.deg_to_km)	
        self.spbRadiusKm.valueChanged.connect(self.km_to_deg)
        self.spbLongMin.valueChanged.connect(self.min_lon_change)
        self.spbLongMax.valueChanged.connect(self.max_lon_change)
        self.spbLatMin.valueChanged.connect(self.min_lat_change)
        self.spbLatMax.valueChanged.connect(self.max_lat_change)		
    def inputs(self):
        # circle parameters
        self.res = self.spbCircleSeg.value()  # circle segments no.
        self.resl = self.spbLineSeg.value()  # line segments no.
        self.radius = self.spbRadiusKm.value()  # Tissot cap radius
        self.delta = self.spbDelta.value() # numerical delta
        # latitude resolution
        self.maxlat = self.spbLatMax.value()
        self.minlat = self.spbLatMin.value()

        self.innerplat = self.spbLatRes.value()

        # longitude resolution
        self.minlon = self.spbLongMin.value()
        self.maxlon = self.spbLongMax.value()

        self.innerplon = self.spbLongRes.value()

        # poles checkbox
        self.a = self.spbSemimajor.value()
        self.b = self.spbSemiminor.value()
        self.auxp = self.chkBoxAux.isChecked()

    def km_to_deg(self):
        dg = self.spbRadiusKm.value() / (2 * (self.a / 1000) * numpy.pi) * 360
        self.spbRadiusDg.setValue(dg)

    def deg_to_km(self):
        km = self.spbRadiusDg.value() / 360 * 2 * (self.a / 1000) * numpy.pi
        self.spbRadiusKm.setValue(km)

    def max_lat_change(self):
        ch = int(self.spbLatMax.value() - 1)
        self.spbLatMin.setMaximum(ch)

    def min_lat_change(self):
        ch = int(self.spbLatMin.value() + 1)
        self.spbLatMax.setMinimum(ch)

    def max_lon_change(self):
        ch = self.spbLongMax.value() - 1
        self.spbLongMin.setMaximum(ch)

    def min_lon_change(self):
        ch = self.spbLongMin.value() + 1
        self.spbLongMax.setMinimum(ch)

    def fun_execute(self):
        self.inputs()
        self.cap_layer()
        self.graticule_layer()

    def cap_layer(self):
        source_crs = self.get_source_crs(self.a, self.b)
        self.vl = QgsVectorLayer("Polygon?crs=" + source_crs.toWkt(), "caps", "memory")
        self.pr = self.vl.dataProvider()
        self.vl.startEditing()
        self.pr.addAttributes([
            QgsField("lon", QVariant.Double),
            QgsField("lat", QVariant.Double),
            QgsField("h", QVariant.Double),
            QgsField("k", QVariant.Double),
            QgsField("s", QVariant.Double),
            QgsField("omega", QVariant.Double)
        ])
        # check auxiliary points
        if self.auxp:
            self.vl2 = QgsVectorLayer("MultiPoint?crs=" + crs.toWkt(), "auxpoints", "memory")
            self.pr2 = self.vl2.dataProvider()
            self.vl2.startEditing()
            self.pr2.addAttributes([QgsField("lon", QVariant.Double), QgsField("lat", QVariant.Double)])

        self.generate_caps()
        self.vl.commitChanges()
        self.vl.updateExtents()
        mNoSimplification = QgsVectorSimplifyMethod()
        mNoSimplification.setSimplifyHints(QgsVectorSimplifyMethod.NoSimplification)
        self.vl.setSimplifyMethod(mNoSimplification)
        self.vl.setCustomProperty("indicatrix", "caps")
        self.vl.setCustomProperty("ellipsoid_a", self.a)
        self.vl.setCustomProperty("ellipsoid_b", self.b)
        self.vl.setCustomProperty("numerical_delta", self.delta)
        
        # Setup labeling
        settings = QgsPalLayerSettings()
        settings.fieldName = "format('h: %1\\nk: %2', round(h, 3), round(k, 3))"
        settings.isExpression = True
        settings.placement = Qgis.LabelPlacement.AroundPoint
        
        text_format = QgsTextFormat()
        text_format.setSize(8)
        settings.setFormat(text_format)
        
        layer_settings = QgsVectorLayerSimpleLabeling(settings)
        self.vl.setLabeling(layer_settings)
        self.vl.setLabelsEnabled(True)
        
        QgsProject.instance().addMapLayer(self.vl)

        if self.auxp:
            self.vl2.commitChanges()
            self.vl2.updateExtents()
            mNoSimplification = QgsVectorSimplifyMethod()
            mNoSimplification.setSimplifyHints(QgsVectorSimplifyMethod.NoSimplification)
            self.vl2.setSimplifyMethod(mNoSimplification)
            QgsProject.instance().addMapLayer(self.vl2)

    def graticule_layer(self):
        source_crs = self.get_source_crs(self.a, self.b)
        self.vl = QgsVectorLayer("LineString?crs=" + source_crs.toWkt(), "graticule", "memory")
        self.pr = self.vl.dataProvider()
        # changes are only possible when editing the layer
        self.vl.startEditing()
        self.pr.addAttributes([QgsField("lon", QVariant.Double), QgsField("lat", QVariant.Double)])
        # calculate graticule lines
        self.generate_graticule()
        # commit to stop editing the layer
        self.vl.commitChanges()
        # update layer's extent when new features have been added
        # because change of extent in provider is not propagated to the layer
        self.vl.updateExtents()
        # no simplify
        mNoSimplification = QgsVectorSimplifyMethod()
        mNoSimplification.setSimplifyHints(QgsVectorSimplifyMethod.NoSimplification)
        self.vl.setSimplifyMethod(mNoSimplification)
        # add layer to the legend
        QgsProject.instance().addMapLayer(self.vl)

    def cap_calc(self):
        lon1 = self.lon * .1745329251994329577e-1 # intial longitude in radians
        lat1 = self.lat * .1745329251994329577e-1  #intial latitude in radians
        s = self.radius * 1000
        a = self.a
        b = self.b
        f = (a - b) / a
        # set up inner points
        pointlist = []
        # necessary azimuth domains
        alphalist = list(numpy.linspace( -numpy.pi, numpy.pi, num = self.res, endpoint = False))
        # credit Michael Kleder, 2007
        for alpha1 in alphalist:
            sinAlpha1 = numpy.sin(alpha1)
            cosAlpha1 = numpy.cos(alpha1)
            tanU1 = (1 - f) * numpy.tan(lat1)
            cosU1 = 1 / numpy.sqrt(1 + tanU1 * tanU1)
            sinU1 = tanU1 * cosU1
            sigma1 = numpy.arctan2(tanU1, cosAlpha1)
            sinAlpha = cosU1 * sinAlpha1
            cosSqAlpha = 1 - sinAlpha * sinAlpha
            uSq = cosSqAlpha * (a * a - b * b) / (b * b)
            A = 1 + uSq / 16384 * (4096 + uSq * ( -768 + uSq * (320 - 175 * uSq)))
            B = uSq / 1024 * (256 + uSq * ( -128 + uSq * (74 - 47 * uSq)))
            sigma = s / (b * A)
            sigmaP = 2 * numpy.pi
            while (numpy.abs(sigma - sigmaP) > 1e-12):
                cos2SigmaM = numpy.cos(2 * sigma1 + sigma)
                sinSigma = numpy.sin(sigma)
                cosSigma = numpy.cos(sigma)
                deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * ( - 3 + 4 * sinSigma * sinSigma) * ( - 3 + 4 * cos2SigmaM * cos2SigmaM)))
                sigmaP = sigma
                sigma = s / (b * A) + deltaSigma
            tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1
            lat2 = numpy.arctan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1, (1 - f) * numpy.sqrt(sinAlpha * sinAlpha + tmp * tmp))
            lam = numpy.arctan2(sinSigma * sinAlpha1, cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1)
            C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
            L = lam - (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * ( - 1 + 2 * cos2SigmaM * cos2SigmaM)))
            lon2 = lon1 + L
            # output degrees
            lat2 = lat2 * 57.295779513082322865
            lon2 = lon2 * 57.295779513082322865
            lon2 = numpy.mod(lon2,360) # follow [0,360] convention
            # longitude correction
            if lon2 > 180:
                lon2 = lon2 - 360
            pointlist.append(QgsPointXY(lon2, lat2))
        
        # Calculate distortion parameters at center
        h, k, s, omega = calculate_indicatrix_params(self.lon, self.lat, self.transform, self.a, self.b, self.delta)
        
        # add a feature
        fet = QgsFeature()
        fet.setGeometry(QgsGeometry.fromPolygonXY([pointlist]))
        fet.setAttributes([float(self.lon), float(self.lat), h, k, s, omega])
        self.pr.addFeatures([fet])
        if self.auxp:
            fet2 = QgsFeature()
            fet2.setGeometry(QgsGeometry.fromMultiPointXY(pointlist))
            fet2.setAttributes([float(self.lon), float(self.lat)])
            self.pr2.addFeatures([fet2])

    def generate_caps(self):
        # Setup transformation to current project CRS for distortion calculation
        target_crs = QgsProject.instance().crs()
        source_crs = self.get_source_crs(self.a, self.b)
        self.transform = QgsCoordinateTransform(source_crs, target_crs, QgsProject.instance())

        # set up inner points
        latlist = numpy.append(numpy.arange(self.maxlat,  self.minlat, -self.innerplat), self.minlat)
        lonlist = numpy.append(numpy.arange(self.minlon, self.maxlon, self.innerplon), self.maxlon)

        if (90 in latlist ): # North pole
            rlist = numpy.array([90])
            latlist = numpy.setdiff1d(latlist,rlist)
            self.lat = 90
            self.lon = 0
            self.cap_calc()

        if (-90 in latlist ): # South pole
            rlist = numpy.array([-90])
            latlist = numpy.setdiff1d(latlist,rlist)
            self.lat = -90
            self.lon = 0
            self.cap_calc()

        for lon in lonlist:
            for lat in latlist:
                self.lat = lat
                self.lon = lon
                self.cap_calc()


    def generate_graticule(self):
        # set up inner points
        latlist = numpy.append(numpy.arange(self.maxlat,  self.minlat, -self.innerplat), self.minlat)
        lonlist = numpy.append(numpy.arange(self.minlon, self.maxlon, self.innerplon), self.maxlon)

        lres = self.resl
        for lon in lonlist:
            vlinelist = list(numpy.linspace(self.minlat, self.maxlat, num=lres))
            vlist = []
            for v in vlinelist:
                vlist.append(QgsPoint(lon, v))
            fetv = QgsFeature()
            fetv.setGeometry(QgsGeometry.fromPolyline(vlist))
            fetv.setAttributes([float(lon), float(0)])
            self.pr.addFeatures([fetv])
        for lat in latlist:
            vlinelist = list(numpy.linspace(self.minlon, self.maxlon, num=lres))
            vlist = []
            for v in vlinelist:
                vlist.append(QgsPoint(v, lat))
            fetv = QgsFeature()
            fetv.setGeometry(QgsGeometry.fromPolyline(vlist))
            fetv.setAttributes([float(0), float(lat)])
            self.pr.addFeatures([fetv])
