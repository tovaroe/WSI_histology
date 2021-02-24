function [] = ndpi_to_tiles(WSI_dir, magnification, tiles_size, output_dir, varargin)
% NDPI_TO_TILES converts ndpi to image tiles
%   
%   ndpi_to_tiles(WSI_DIR, MAGNIFICATION, TILES_SIZE, OUTPUT_DIR, ...)
%   opens the ndpi-file at WSI_DIR and saves it as multiple image tiles at
%   the given MAGNIFICATION and given TILES_SIZE to the output directory
%   OUTPUT_DIR.
%
%OPTIONAL PARAMETERS
%useAnnotation - if you only want to save tiles that correspond to a given
%       freehand annotation, set this to the corresponding integer (see
%       extract_annotation function) (default NaN)
%makeAnnotationCSV - if set to true, saves a csv with the relative area of
%       every annotation for every saves image tile (default false)
%filterFunction - set this to a custom filter function, that takes a double
%       image as argument and returns true (save the image tile) or false
%       (skip the image tile). It is highly recommended to implement this
%       custom filter function, as the default function might not work well
%       for all purposes and magnifications.
%minRelativeTissueArea - the minimum relative tissue area each image tile
%       needs to have in order to be saved (ignores all tiles with less
%       relative tissue area) (default NaN)
%segmentTissue - custom function for determining which pixels correspond to
%       tissue on each image tile. Should take an image and return a
%       boolean matrix. Defining a custom function is highly recommended.
%overlap - overlap of the image tiles in px (default [0 0])
%splitHE - if true, split image into Haematoxylin and Eosin color channels 
%       before saving (default false)
%
% Relies on Daniel Forsberg's openslide MATLAB implementation: 
% https://github.com/openslide/openslide
%
%
% Copyright (c) 2021 Thomas Roetzer, MedUni Vienna
% thomas.roetzer@meduniwien.ac.at
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.    

%%
% See input parser for arguments

%% Parse inputs
defaultUseAnnotation = NaN;
defaultMakeAnnotationCSV = false;
defaultFilterFunction = @(x) true;
defaultMinRelativeTissueArea = NaN;
defaultSegmentTissue = @(x) segment_tissue(x);
defaultOverlap = [0 0];

p = inputParser;

addParameter(p, 'useAnnotation', defaultUseAnnotation);
addParameter(p, 'makeAnnotationCSV', defaultMakeAnnotationCSV);
addParameter(p, 'segmentTissue', defaultSegmentTissue);
addParameter(p, 'filterFunction', defaultFilterFunction);
addParameter(p, 'minRelativeTissueArea', defaultMinRelativeTissueArea);
addParameter(p, 'overlap', defaultOverlap);
addParameter(p, 'splitHE', false);

parse(p, varargin{:});

% Use annotation
useAnnotation = p.Results.useAnnotation;
makeAnnotationCSV = p.Results.makeAnnotationCSV;
segmentTissue = p.Results.segmentTissue;
minRelativeTissueArea = p.Results.minRelativeTissueArea;
filterFunction = p.Results.filterFunction;
overlap = p.Results.overlap;
splitHE = p.Results.splitHE;

%%

% Get the annotation
if (~isnan(useAnnotation) || makeAnnotationCSV)
    annotation_path = [WSI_dir '.ndpa'];
    if exist(annotation_path, 'file')
        annotation = extract_annotation(annotation_path, magnification);
    else
        annotation = NaN;
    end
    annotationTable = table;
end

% Open whole-slide image
slidePtr = openslide_open(WSI_dir);

% Get filename
path_split = split(WSI_dir, filesep);
filename = char(path_split(end));
savename = strrep(filename, '.ndpi', '');

% Get whole-slide image properties
[~, ~, ~, ~, numberOfLevels, ...
    downsampleFactors, objectivePower] = openslide_get_slide_properties(slidePtr);

% Get available magnifications
WSI_magnifications = objectivePower./downsampleFactors;

% Check if given magnification is available in the WSI
if ~sum(WSI_magnifications == magnification)
    error([num2str(magnification) 'x-magnification not available!']);
end

% Get the 'level' of the WSI with the given magnification
for i=1:numberOfLevels
    if WSI_magnifications(i) == magnification
        level = i-1;
    end
end

% Get the width and the height of the 'level'
[width, height] = openslide_get_level_dimensions(slidePtr,level);

% Extract the tiles
tile_height = tiles_size(1);
tile_width = tiles_size(2);
%tiles_y = floor(double(height)/(tile_height-overlap(1)))-1;
tiles_y = 1+floor((double(height)-tile_height-1)/(tile_height-overlap(1)));
%tiles_x = floor(double(width)/(tile_width-overlap(2)))-1;
tiles_x = 1+floor((double(width)-tile_width-1)/(tile_width-overlap(2)));

id = 0;

for y=1:tiles_y
    for x=1:tiles_x
        
        
        y_pos_start = (y-1)*(tile_height-overlap(1))+1;
        x_pos_start = (x-1)*(tile_width-overlap(2))+1;
        
        
        % Use annotation
        if ~isnan(useAnnotation)
            correct_annotation_matrix = annotation(y_pos_start:y_pos_start+tile_height-1, x_pos_start:x_pos_start+tile_width-1) == useAnnotation;
            if mean(correct_annotation_matrix(:))<.9 
                continue
            end
        end
        
        % Read tile from slidescan
        tile = openslide_read_region(slidePtr,x_pos_start,y_pos_start,tile_width,tile_height, 'level', level);
        tile = tile(:,:,2:4);
                
        if ~filterFunction(double(tile))
            continue
        end
        
        if ~isnan(minRelativeTissueArea) && ~makeAnnotationCSV
            tissue = segmentTissue(tile);
            if mean(tissue(:)) < minRelativeTissueArea
                continue
            end
        end
                
        % Make annotation csv
        if makeAnnotationCSV
            if isnan(annotation)
                annotation_matrix = zeros(tile_height, tile_width);
            else
                annotation_matrix = annotation(y_pos_start:y_pos_start+tile_height-1, x_pos_start:x_pos_start+tile_width-1);
            end

            % 1 ... RedMask
            red = mean(mean(annotation_matrix==1));
            % 2 ... YellowMask
            yellow = mean(mean(annotation_matrix==2));
            % 3 ... GreenMask
            green = mean(mean(annotation_matrix==3));
            % 4 ... BlueMask
            blue = mean(mean(annotation_matrix==4));
            % 5 ... BlackMask
            black = mean(mean(annotation_matrix==5));
            % 6 ... MagentaMask
            magenta = mean(mean(annotation_matrix==6));
            % 7 ... TealMask
            teal = mean(mean(annotation_matrix==7));
            % 8 ... WhiteMask
            white = mean(mean(annotation_matrix==8));
            
            % Tumor from mask
            tissue = segmentTissue(tile);
            tumor = mean(mean(tissue & ~(annotation_matrix)));
            
            if ~isnan(minRelativeTissueArea)
                tissue = tissue | annotation_matrix;
                if mean(tissue(:)) < minRelativeTissueArea
                    continue
                end
            end
                
            id = id + 1;
            annotationTable(size(annotationTable,1)+1, :) = {num2str(id) tumor red yellow green blue black magenta teal white};
        end

        if splitHE
            CD = colour_deconvolution(tile, 'HE');
            H = CD(:,:,1);
            E = CD(:,:,2);
            tilenameH = [savename '_'  num2str(id) '_H' '.png'];
            tilenameE = [savename '_'  num2str(id) '_E' '.png'];
            imwrite(H, [output_dir,filesep,tilenameH]);
            imwrite(E, [output_dir,filesep,tilenameE]);
            
        else                       
            tilename = [savename '_'  num2str(id) '.png'];
            imwrite(tile, [output_dir,filesep,tilename]);
        end
    end
end

% save annotation csv
if makeAnnotationCSV
    if size(annotationTable, 1) > 0
        annotationTable.Properties.VariableNames = [{'id'} {'tissue'} {'red'} {'yellow'} {'green'} {'blue'} {'black'} {'magenta'} {'teal'} {'white'}];
        writetable(annotationTable, [output_dir, filesep, savename, 'annotation.csv']);
    end
end
    
% Clear the pointer
openslide_close(slidePtr)
clear slidePtr

end
