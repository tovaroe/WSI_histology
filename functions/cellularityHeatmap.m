function heatmap = cellularityHeatmap(centroidMask,tissueMask, magnification)
%CELLULARITYHEATMAP returns a double image
%
%   heatmap = cellularityHeatmap(CM, TM, MAG)
%   estimates a pixel-wise cellularity from gaussian distributions around
%   each cell centroid in the centroid mask CM, with boundaries defined by 
%   the tissue mask TM, and for a magnification (as defined by ndpi scans)
%   of MAG. The gaussian filtering is performed such that 99% of the total 
%   sum lies in a circle with radius 1mm² around each centroid.
%   The returned heatmap consists of the pixel-wise cellularites per mm².

% Copyright (C) 2021 Thomas Roetzer, MedUni Vienna
% thomas.roetzer@meduniwien.ac.at

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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% Only pixels with minNeighbors in the given strel st survive
%%
nm_per_px = 228*40/magnification;

% Testing gaussian approach
radius = sqrt(1e12/pi)/nm_per_px; % radius of circle with 1mm² in px
gauss_std = radius/sqrt(9.21); % (R: pchisq(9.21, df=1:10))

heatmap = imgaussfilt(im2double(centroidMask), gauss_std) / (nm_per_px/1e6)^2;

% Correcting the boarder by calculating applicable area for density
% estimation
tissuearea = imgaussfilt(im2double(tissueMask), gauss_std);
tissuearea(~tissueMask) = 1;

heatmap = heatmap./tissuearea;
end

