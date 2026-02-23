# Indicatrix Mapper (Quasi-Tissot Indicatrices)

A QGIS plugin to visualize and analyze map projection distortions using **Quasi-Tissot Indicatrices**.

![Indicatrix Mapper Icon](icon.png)

## Overview

The **Indicatrix Mapper** generates dynamic ellipses that represent the local distortion of map projections. Using numerical derivatives (the **Quasi Indicatrix** method), this plugin calculates the exact scale factors ($h, k$), area distortion ($s$), angular deformation ($\omega$), and the principal axes ($a, b$) with their specific orientation.

These indicatrices are rendered using QGIS point symbology, allowing them to:
- Update **instantaneously** when the project CRS is changed.
- Maintain **analytical precision** even for complex or custom projections.
- Be styled and analyzed using standard QGIS layer properties.

## Key Features

- **Real-time Distortions:** Scale, area, and angular distortion attributes are recalculated automatically whenever you change the map projection.
- **Dynamic Symbology:** Ellipses are rotated and scaled via data-defined overrides based on the calculated principal axes of the Tissot indicatrix.
- **Customizable Grids:** Define the extent and resolution of the indicatrix grid.
- **Modern API Compatibility:** Fully compatible with QGIS 3.x (using `QMetaType`, `setStrokeColor`, and other modern classes).
- **Comprehensive Data:** The attribute table provides $h, k, s, \omega, a_{ind}, b_{ind}$ and the rotation angle for every point.

## Theoretical Background

The plugin implements the Quasi Indicatrix method based on:
> Bildirici, Ö. G. (2015). *Numerical calculation of Tissot indicatrix parameters for map projections.*

This approach allows the plugin to handle any projection supported by QGIS/PROJ without requiring complex symbolic differentiation.

## Installation

1. Clone or download this repository.
2. Link or copy the folder into your QGIS profile's plugin directory:
   - `AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\IndicatrixMapper`
3. Restart QGIS and enable the plugin from the **Plugin Manager**.

## Usage

1. Open the **Indicatrix Mapper** from the Vector menu or toolbar.
2. Set your desired extent and resolution (latitude/longitude steps).
3. Define the **radius** (in km or degrees) to control the visual size of the ellipses.
4. Click **Run!** to generate the indicatrix and graticule layers.
5. Change the **Project CRS** (bottom-right of the QGIS window) to see the ellipses rotate and scale in real-time.

---
**Author:** Ervin Wirth  
**Contact:** wirth.ervin@mailbox.org
