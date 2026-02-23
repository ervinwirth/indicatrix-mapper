# Development & Build Instructions

This document provides technical details for building, testing, and contributing to the **Indicatrix Mapper** QGIS plugin.

## Prerequisites
- **QGIS 3.x:** The plugin is written in Python 3 and uses the PyQGIS API.
- **Python 3:** Bundled with QGIS.
- **PyQt5:** For the user interface (included with QGIS).
- **NumPy:** For mathematical calculations (included with QGIS).

## Environment Setup
The easiest way to develop for QGIS is to use the **OSGeo4W Shell** (Windows) or the system terminal with QGIS installed (Linux/macOS).

### Windows (OSGeo4W)
1. Open the OSGeo4W Shell.
2. Run `qgis` to start QGIS or `python-qgis` for a QGIS-aware Python environment.

## Build Process
As a Python-based plugin, there is no "compilation" in the traditional sense. However, the UI and resource files must be compiled for use in Python.

### Compiling UI Files
If you modify the `.ui` file (using Qt Designer), recompile it with:
```bash
pyuic5 indicatrix_mapper_dialog_base.ui -o indicatrix_mapper_dialog_base.py
```

### Compiling Resource Files
If you add new icons or assets to `resources.qrc`, recompile it with:
```bash
pyrcc5 resources.qrc -o resources_rc.py
```

## Testing
Currently, the plugin does not have an automated test suite. (Planned: Add unit tests for `cap_calc`).

### Manual Testing
1. Copy the plugin folder to your QGIS profile:
   - Windows: `%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\tiss`
   - Linux: `~/.local/share/QGIS/QGIS3/profiles\default\python\plugins\tiss`
2. In QGIS, use the **Plugin Reloader** plugin to quickly reload your changes without restarting QGIS.
3. Open the "Indicatrix mapper" dialog from the Vector menu.
4. Set an extent and click **Run!**.
5. Change the Project CRS (bottom-right of the QGIS window) and observe the dynamic refresh of the ellipses.

## Packaging for QGIS Repository
The QGIS Plugin Repository is strict about the zip file structure (it must contain a single root folder, e.g., `tiss`). To ensure maximum compatibility and avoid "Cannot find a folder" errors, use the following Python command to create the release package:

```bash
python -c "import zipfile, os; files = ['__init__.py', 'indicatrix_mapper.py', 'indicatrix_mapper_dialog.py', 'indicatrix_mapper_dialog_base.py', 'metadata.txt', 'resources_rc.py', 'icon.png', 'graticule.png', 'README.md', 'LICENSE']; zip_name = 'tiss-v3.1.0.zip'; [zipf.write(f, os.path.join('tiss', f)) for zipf in [zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED)] for f in files if os.path.exists(f)]"
```

## Code Structure
- `indicatrix_mapper.py`: Main entry point, GUI setup, and CRS change listener.
- `indicatrix_mapper_dialog.py`: Core logic for Quasi-Tissot parameters and layer creation.
- `indicatrix_mapper_dialog_base.py`: Auto-generated UI code.
- `resources_rc.py`: Auto-generated resource code.
