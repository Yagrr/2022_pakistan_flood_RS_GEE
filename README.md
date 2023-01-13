# 2022 Pakistan Flood - Remote Sensing project (GEE)
Agricultural flood impact analysis of Sindh province during flooding in 1st-31st August.

Main results are in the 'Pakistan GEE Results 250m Statistics' folder. 

Results:

Sindh_District_results.gpkg
Sindh_Tehsil_results.gpkg 

All data and software used in this project are open-data or publicly available.

Project files are under MIT license and can be used as you please.
By this license, you must credit me as (Laurent Chan) and link back to this repository.

## Data sources
.tif files are exports from Google Earth Engine (GEE) processed using QGIS.

Pakistan boundary files were obtained from:
https://data.humdata.org/dataset/cod-ab-pak

Waterways and Indus river .shp files were obtained from QGIS' QuickOSM plugin:
OpenStreetMap contributors, ‘Planet dump retrieved from https://planet.osm.org’, 2015. https://planet.openstreetmap.org/

## Source code citations
The Refined Lee filter used in the source code was implemented in GEE by Guido Lemoine:
https://code.earthengine.google.com/2ef38463ebaf5ae133a478f173fd0ab5

Edge Otsu algorithm adapted from Markert et al. (2020) and mygeoblog.com:
https://doi.org/10.3390/rs12152469
https://mygeoblog.com/2021/01/25/edge-otsu-for-surface-water-mapping-detection/

The display product section was adapted from UN Spider:
https://code.earthengine.google.com/f5c2f984c053c8ea574bfcd4040d084e
