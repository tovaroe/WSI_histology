# WSI_histology
MATLAB functionality for working with whole slide images of histological slide scans.

## 1. Converting a WSI (ndpi) to multiple smaller image tiles
1. Run scripts/multi_ndpi_to_tiles.m
2. Select all the ndpis you want to save as tiles in the file selection dialog box
3. Select the output folder, where tiles will be saved

The image tiles will be saved into a new folder located in the output folder for each WSI. If you set makeAnnotationCSV to true, an accompanying csv which specifies the relative area of each annotation for each tile will be saved.


## Dependencies
MATLAB bindings for the openslide library by Daniel Forsberg: https://github.com/fordanic/openslide-matlab
