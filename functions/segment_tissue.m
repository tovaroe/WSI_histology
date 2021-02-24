function tissue_mask = segment_tissue(image, varargin)
%SEGMENT_TISSUE returns a binary tissue mask
%
%   tissue_mask = segment_tissue(IM, VARARGIN) performs a quick
%   segmentation of tissue in the input image IM and returns a binary mask.
%   
%OPTIONAL PARAMETERS
%stain - the histological stains that should be used for colour 
%       deconvolution used [ {'HE'} | 'H DAB' | 'H DAB 2' | 'H DAB 
%       Ventana']
%areaChannels - color channels used to obtain the tissue mask (default 1:2)
%areaThreshold - threshold used for obtaining the tissue mask (default 200)
%areaClosureStrel - strel used for binary closure of the tissue mask
%       (default strel('disk', 30)
%tissueSizeMin - minimum pixel size of individual tissue fragments to be 
%       kept in the tissue mask (default 100)
%showFigures - show original image and segmentation side-by-side (default
%       false)

%Copyright (c) 2021 Thomas Roetzer, MedUni Vienna
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
defaultStain = 'HE';
expectedStains = {'HE', 'H DAB', 'H DAB 2', 'H DAB Ventana'};
defaultAreaChannels = 1:2;
defaultAreaThreshold = 200;
defaultAreaClosureStrel = strel('disk', 30);
defaultTissueSizeMin = 100;

p = inputParser;

validImage = @(x) isa(x, 'uint8') && size(x,3)==3;
addRequired(p,'image', validImage);
addParameter(p, 'stain', defaultStain, @(x) any(validatestring(x, expectedStains)));
addParameter(p, 'areaChannels', defaultAreaChannels);
addParameter(p, 'areaThreshold', defaultAreaThreshold);
addParameter(p, 'areaClosureStrel', defaultAreaClosureStrel);
addParameter(p, 'tissueSizeMin', defaultTissueSizeMin);
addParameter(p, 'showFigures', false);

parse(p, image, varargin{:});

%% Pre
CD =  colour_deconvolution(image, p.Results.stain);

%% Area 
tissue_mask = sum(CD(:,:,p.Results.areaChannels), 3)/length(p.Results.areaChannels) < p.Results.areaThreshold;
tissue_mask = imclose(tissue_mask, p.Results.areaClosureStrel);
tissue_mask = bwareaopen(tissue_mask, p.Results.tissueSizeMin);

if p.Results.showFigures
    figure(1), imshow(image), title("Original image");
    figure(2), imshow(tissue_mask), title("Tissue mask");
end

end

