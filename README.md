# WSI_histology
MATLAB functionality for working with whole slide images of histological slide scans.

## 1. Converting a WSI (ndpi) to multiple smaller image tiles
1. Run scripts/multi_ndpi_to_tiles.m
2. Select all the ndpis you want to save as tiles in the file selection dialog box
3. Select the output folder, where tiles will be saved

The image tiles will be saved into a new folder located in the output folder for each WSI. If you set makeAnnotationCSV to true, an accompanying csv which specifies the relative area of each annotation for each tile will be saved.


## 2. Calculating celluarities and tissue areas for multiple WSIs (ndpi)
1. Run scripts/multi_cellularityWSI.m
2. Select all the ndpi for which you want to calculate the cellularity in the file selection dialog box

The cellularities of all WSIs will be saved as "HE_analysis.csv" in the same folder as the WSIs. Additionally, by default, cellularity heatmaps, cell centroid masks and tissue masks will be saved into an 'Analyses' folder.

## Dependencies
MATLAB bindings for the openslide library by Daniel Forsberg: https://github.com/fordanic/openslide-matlab

For working with annotations: xml2struct https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct

## Citing
If you find any of these function useful, please cite: 
Roetzer-Pejrimovsky, T., Moser, AC., Atli, B. et al. The Digital Brain Tumour Atlas, an open histopathology resource. Sci Data 9, 55 (2022). 
https://doi.org/10.1038/s41597-022-01157-0


