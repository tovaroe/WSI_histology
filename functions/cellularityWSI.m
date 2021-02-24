function [] = cellularityWSI(FileName,PathName, varargin)
%CELLULARITYWSI saves cellularities to a csv
%
%   cellularityWSI(FN, PN, ...) reads the ndpi-WSI saved at the path PN and
%   calculates the cellularity for the whole slide and annotated regions,
%   if present. It saves the data as 'HE_Analysis.csv' in the same folder
%   or appends the data to an already existing 'HE_analysis.csv'.
%
%
%OPTIONAL PARAMETERS
%useAnnotation - if true, performs cellularity calculation for all freehand
%    annotations in the corresponding ndpa-file, but uses up more memory.
%    (default true)
%makeHeatmap - if true, saves a cellularity Heatmap into an
%   Analyses-folder. (default true)
%saveMasks - if true, saves centroid mask and tissue mask into an
%   Analyses-folder
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
%% Parse inputs
defaultUseAnnotation = true;
defaultMakeHeatmap = true;
defaultSaveMasks = true;

p = inputParser;

addParameter(p, 'useAnnotation', defaultUseAnnotation);
addParameter(p, 'makeHeatmap', defaultMakeHeatmap);
addParameter(p, 'saveMasks', defaultSaveMasks);

parse(p, varargin{:});

% Use annotation
useAnnotation = p.Results.useAnnotation;
makeHeatmap = p.Results.makeHeatmap;
saveMasks = p.Results.saveMasks;

%%
% Calulates only cellularity for whole slide + heatmap

if ~exist([PathName filesep 'Analyses'], 'dir')
    mkdir([PathName filesep 'Analyses']);
end

%% Load slide scans and annotation data
% Input file names for image and xml file
NDPAFileName = strrep(FileName, '.ndpi', '.ndpi.ndpa');

% Segment cells and tissue
image = extract_slidescan_image([PathName filesep FileName], 10);
disp('loaded')
[CellMask, TissueMask] = segment_cells(image);
TissueMask = bwareaopen(TissueMask,100000); 
CellMask = CellMask{1};

% Allocate annotation masks
if useAnnotation
    if isfile([PathName filesep NDPAFileName])
        annotation = extract_annotation([PathName filesep NDPAFileName], 10);

        GreenMask = annotation == 3; 
        RedMask = annotation == 1; 
        YellowMask = annotation == 2;  
        BlueMask = annotation == 4;  
        BlackMask = annotation == 5;
        MagentaMask = annotation == 6;
        TealMask = annotation == 7;
        WhiteMask = annotation == 8;
    else
        warning(strcat(FileName, ': no annotation available.'))
        useAnnotation = false;
    end
end
    
% Determine tumor tissue
if useAnnotation
    UndefinedMask = (TissueMask & ~(RedMask | YellowMask | GreenMask | BlueMask | BlackMask | MagentaMask | TealMask | WhiteMask));
end

CellMask(~TissueMask) = 0;

% Areas in mm²
PixelSize = ((1.824e-3)/2)^2;  % Area of one pixel in exported 10x image [mm²]
Atissue = sum(TissueMask(:))*PixelSize;

if useAnnotation
    Ayellow = sum(YellowMask(:))*PixelSize;
    Agreen = sum(GreenMask(:))*PixelSize;
    Ablack = sum(BlackMask(:))*PixelSize;
    Ablue = sum(BlueMask(:))*PixelSize;
    Ared = sum(RedMask(:))*PixelSize;
    Amagenta = sum(MagentaMask(:))*PixelSize;
    Awhite = sum(WhiteMask(:))*PixelSize;
    Ateal = sum(TealMask(:))*PixelSize;
    Aundefined = sum(UndefinedMask(:))*PixelSize;
else
    Ayellow = 0;
    Agreen = 0;
    Ablack = 0;
    Ablue = 0;
    Ared = 0;
    Amagenta = 0;
    Awhite = 0;
    Ateal = 0;
    Aundefined = Atissue;
end

%% Cellularity map generation
centroidMask = bwmorph(imfill(CellMask, 'holes'), 'shrink', Inf);
cellularity_tissue = sum(centroidMask(:))/Atissue;

if makeHeatmap
    heatmap = cellularityHeatmap(centroidMask, TissueMask, 10);
    heatmap(~TissueMask) = 0;
end

if useAnnotation
    cellularity_black = sum(sum(centroidMask & BlackMask))/Ablack;
    cellularity_green = sum(sum(centroidMask & GreenMask))/Agreen;
    cellularity_yellow = sum(sum(centroidMask & YellowMask))/Ayellow;
    cellularity_red = sum(sum(centroidMask & RedMask))/Ared;
    cellularity_blue = sum(sum(centroidMask & BlueMask))/Ablue;
    cellularity_magenta = sum(sum(centroidMask & MagentaMask))/Amagenta;
    cellularity_teal = sum(sum(centroidMask & TealMask))/Ateal;
    cellularity_white = sum(sum(centroidMask & WhiteMask))/Awhite;
    cellularity_undefined = sum(sum(centroidMask & UndefinedMask))/Aundefined;
else
    cellularity_black = NaN;
    cellularity_green = NaN;
    cellularity_yellow = NaN;
    cellularity_red = NaN;
    cellularity_blue = NaN;
    cellularity_magenta = NaN;
    cellularity_teal = NaN;
    cellularity_white = NaN;
    cellularity_undefined = cellularity_tissue;
end


%% Store native nuclei centroids & heatmap
StripName = strrep(FileName, '.ndpi','');
if saveMasks 
    save([PathName filesep 'Analyses' filesep StripName 'centroids' '.mat'], 'centroidMask'); 
    save([PathName filesep 'Analyses' filesep StripName 'tissuemask' '.mat'], 'TissueMask'); 
end

if makeHeatmap 
    imagesc(heatmap), colorbar, saveas(gcf, [PathName filesep 'Analyses' filesep StripName 'nuclei_heatmap' '.png']);
    save([PathName filesep 'Analyses' filesep StripName 'nuclei_heatmap.mat'], 'heatmap');
end

%% Save to csv
id = {FileName};
annotationTable = table(id, cellularity_tissue, cellularity_undefined, cellularity_black, cellularity_blue, cellularity_green, cellularity_magenta, cellularity_red, cellularity_teal, cellularity_white, cellularity_yellow, Atissue, Aundefined, Ablack, Ablue, Agreen, Amagenta, Ared, Ateal, Atissue, Awhite, Ayellow);

if exist([PathName filesep 'HE_analysis.csv'], 'file')
    annotationTable = [readtable([PathName filesep 'HE_analysis.csv']); annotationTable];
end

writetable(annotationTable, [PathName, filesep, 'HE_analysis.csv']);

%% End massage
disp([StripName ' evaluated and saved']);