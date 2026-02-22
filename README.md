# Indicatrix Mapper (QGIS Plugin)

The plugin introduces ellipsoidal or spherical caps which can give a **down-to-earth** Tissot-indicatrix realization.

Unlike traditional infinitesimal Tissot circles, this plugin uses constant-radius ellipsoidal caps calculated using **Vincenty's formula**. These caps are generated as memory layers and can transform on-the-fly from a reference ellipsoid to any selected project coordinate reference system (CRS). This allows users to visually analyze if a projection is conformal, equal-area, or otherwise distorted.

## Features
- **Geodesic Caps:** Accurate ellipsoidal caps generated using Vincenty's formula.
- **On-the-fly Distortion:** Visualize projection warping in real-time using QGIS's projection engine.
- **Analytical Parameters:** Dynamic calculation of $h, k, s, \omega$ distortion values.
- **Customizable Extent:** Define the geographic range and resolution of the indicatrix grid.
- **Graticule Generation:** Create reference lon/lat lines alongside the caps.

## Scientific Background
The analytical distortion parameters ($h, k, s, \omega$) are calculated using numerical differentiation of the forward projection. This approach is known as the **Quasi Indicatrix** method:
- *Bildirici, Ö. G. (2015). Quasi-indicatrix for distortion analysis of map projections. Survey Review.*

## Installation
1. Download or clone this repository.
2. Copy the folder to your QGIS plugin directory:
   - **Windows:** `%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\indicatrix-mapper`
   - **Linux/macOS:** `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/indicatrix-mapper`
3. Restart QGIS.
4. Enable "Indicatrix mapper" in **Plugins -> Manage and Install Plugins**.

## Documentation
- [Development and Build Instructions](DEVELOPMENT.md)
- [Original Idea and Research](IDEA.md)
