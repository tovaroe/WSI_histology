function [cell_masks,tissue_mask] = segment_cells(image, varargin)
%SEGMENT_CELLS Returns binary cell image(s) & tissue image
%
%   [cell_masks, tissue_mask] = segment_cells(IM, 'PARAM1', 'VALUE1', ...)
%   processes the image IM to obtain cell masks of the given color channels
%   and a tissue mask of the native input image. In default mode, this
%   function assumes an H&E-staining and segments the Haematoxylin-stained
%   nuclei only. Note that cell_masks is returned as a cell of one or
%   multiple matrices.
%
%OPTIONAL PARAMETERS
%stain - the histological stains that should be used for colour 
%       deconvolution used [ {'HE'} | 'H DAB' | 'H DAB2' | 'H DAB Ventana']
%areaChannels - color channels used to obtain the tissue mask (default 1:2)
%areaThreshold - threshold used for obtaining the tissue mask (default 200)
%areaClosureStrel - strel used for binary closure of the tissue mask
%       (default strel('disk', 30)
%tissueSizeMin - minimum pixel size of individual tissue fragments to be 
%       kept in the tissue mask (default 1000)
%showFigures - boolean if figures for each segmentation step should be
%       showed (default false)
%
%channels - cell array with one cell for each channel (or any combination
%       of channels) to be combined when constructing its respective cell
%       mask. I.e., the number of cells will be equal to the number of
%       returned cell_masks. 
%
%EACH OF THE FOLLOWING PARAMETERS MUST BE GIVEN AS A LIST WITH A LENGTH
%EQUAL TO THE LENGHT OF THE CHANNELS-ARRAY
%
%globalThreshold - global theshold for cell mask generation (default 120)
%otsuThresholdFactor - factor to adapt Otsu thresholding (default 1)
%phansalkarRadius - radius used for calcualation of Phansalkar threshold
%       (default 8)
%phansalkar_k - paramater k for calculation of Phansalkar threshold 
%       (default 0.2)
%cellSizeMin - minimum size of cells to be kept in the cell mask (default 
%       4)
%watershedSigma - sigma used for watershedding (default 1)
%survivalNeighbors - value used for filtering out pixels in the binary 
%       cell_mask that have less than survivalNeighbors neighbouring pixels
%       (default 2)
%survivalNeighborhood - strel used for determining which pixel qualify as
%       'neighboring' (default strel('disk',1))
%
% Copyright (c) 2019 Thomas Roetzer, MedUni Vienna
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
defaultStain = 'HE';
expectedStains = {'HE', 'H DAB', 'H DAB 2', 'H DAB Ventana', 'HE Pigment'};
defaultAreaChannels = 1:2;
defaultAreaThreshold = 220;
defaultAreaClosureStrel = strel('disk', 30);
defaultTissueSizeMin = 1000;
defaultGlobalThreshold = [120 120 120];
defaultOtsuThresholdFactor = [1 1 1];
defaultPhansalkarRadius = [8 8 8];
defaultPhansalkar_k = [0.2 0.2 0.2];
defaultCellSizeMin = [4 4 4];
defaultWatershedSigma = [1 1 1];
defaultVesselStainCorrection = [false false false]; %sometimes, vessels get aberrantly stained like CD34. This attempts a correction.
defaultSurvivalNeighbors = [2 2 2];
defaultSurvivalNeighborhood = [strel('disk', 1) strel('disk', 1) strel('disk', 1)];
defaultChannels = {1}; %e.g. {1 1:2}, {1 3}, etc

p = inputParser;

validImage = @(x) isa(x, 'uint8') && size(x,3)==3;
addRequired(p,'image', validImage);
addParameter(p, 'stain', defaultStain, @(x) any(validatestring(x, expectedStains)));
addParameter(p, 'areaChannels', defaultAreaChannels);
addParameter(p, 'areaThreshold', defaultAreaThreshold);
addParameter(p, 'areaClosureStrel', defaultAreaClosureStrel);
addParameter(p, 'tissueSizeMin', defaultTissueSizeMin);
addParameter(p, 'globalThreshold', defaultGlobalThreshold);
addParameter(p, 'otsuThresholdFactor', defaultOtsuThresholdFactor);
addParameter(p, 'doPhansalkar', [true true true]);
addParameter(p, 'phansalkarRadius', defaultPhansalkarRadius);
addParameter(p, 'phansalkar_k', defaultPhansalkar_k);
addParameter(p, 'cellSizeMin', defaultCellSizeMin);
addParameter(p, 'watershedSigma', defaultWatershedSigma);
addParameter(p, 'vesselStainCorrection', defaultVesselStainCorrection);
addParameter(p, 'survivalNeighbors', defaultSurvivalNeighbors);
addParameter(p, 'survivalNeighborhood', defaultSurvivalNeighborhood);
addParameter(p, 'channels', defaultChannels);
addParameter(p, 'showFigures', false);

parse(p, image, varargin{:});

%% Pre
cell_masks = {};
CD =  colour_deconvolution(image, p.Results.stain);

%% Area 
tissue_mask = sum(CD(:,:,p.Results.areaChannels), 3)/length(p.Results.areaChannels) < p.Results.areaThreshold;
tissue_mask = imclose(tissue_mask, p.Results.areaClosureStrel);
tissue_mask = bwareaopen(tissue_mask, p.Results.tissueSizeMin);

if p.Results.showFigures
    figure(1), imshow(image), title("Original image");
    figure(2), imshow(tissue_mask), title("Tissue mask");
end

for channel=1:length(p.Results.channels)
    I = median(CD(:,:,p.Results.channels{channel}), 3);

    %% Cellularity

    % Global thresholding
    global_bw = I < p.Results.globalThreshold(channel);
    bwareaopen(global_bw,4);
    global_bw = imdilate(global_bw, strel('square', 5));

    % Otsu
    counts = imhist(I(tissue_mask), 64);
    T = otsuthresh(counts)*p.Results.otsuThresholdFactor(channel);
    otsu = ~imbinarize(I,T);
    otsu = imfill(otsu, 'holes');

    % Phansalkar
    if p.Results.doPhansalkar(channel)
        phans = ~phansalkar(I,p.Results.phansalkarRadius(channel),p.Results.phansalkar_k(channel));
        if p.Results.vesselStainCorrection(channel)
            cc = bwconncomp(phans);
            stats = regionprops(cc, 'Area', 'Solidity', 'PixelIdxList');
            idx = find([stats.Solidity] > 0.8 | [stats.Area] < 500); % only keeps small or solid objects
            phans = ismember(labelmatrix(cc), idx);
        end
        phans = watershed_segmentation(phans, p.Results.watershedSigma(channel));
    else 
        phans = ones(size(I));
    end

    % Combine
    cell_mask= phans & otsu & global_bw;

    if p.Results.survivalNeighbors
        cell_mask = bwultsurvive(cell_mask, p.Results.survivalNeighborhood(channel), p.Results.survivalNeighbors(channel)); % This is similar to a custom erosion, where only pixels with at least two neighbours survive
    end

    cell_mask = bwareaopen(cell_mask, p.Results.cellSizeMin(channel));
    cell_masks{channel} = cell_mask;

    % Plot boundaries for test purposes
    if p.Results.showFigures 
        boundaries = bwboundaries(cell_mask);
        image_control = image;
        for k=1:length(boundaries)
            boundary = boundaries{k};
            for j=1:length(boundary)
                px = boundary(j,:);
                image_control(px(1), px(2), :)=255;
            end
        end
        figure((channel-1)*4+3), imshow(image_control), title("Segmentation channel " + num2str(channel));
        figure((channel-1)*4+4), imshow(global_bw), title("Global threshold channel " + num2str(channel));
        figure((channel-1)*4+5), imshow(otsu), title("Otsu channel " + num2str(channel));
        figure((channel-1)*4+6), imshow(phans), title("Phansalkar channel " + num2str(channel));            
    end      
end

