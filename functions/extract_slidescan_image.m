function image = extract_slidescan_image(WSI, magnification)
%EXTRACT_SLIDESCAN_IMAGE returns a uint8 image matrix
%   
%   extract_slidescan_image(WSI, MAG) reads the ndpi at the file path WSI
%   and returns a uint8 image matrix in magnification MAG.
%
% Relies on Daniel Forsberg's openslide MATLAB implementation: 
% https://github.com/openslide/openslide
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
% Open whole-slide image
slidePtr = openslide_open(WSI);

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

% Read the whole image
ARGB = openslide_read_region(slidePtr,0,0,width,height, 'level', level);

% Output only the RGB-part of the image
image = ARGB(:,:,2:4);

% Clear the pointer
clear slidePtr

end

