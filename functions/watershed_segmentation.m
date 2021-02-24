function I = watershed_segmentation(I, varargin)
%WATERSHED_SEGMENTATION return a binary image with split cells
%
%   I = watershed_segmentation(I) used the watershed algorithm [1] to
%   split closely adjacent cells (blobs) in the input image I.
%
%   I = watershed_segmentation(I, sigma) uses the given sigma for the
%   gauss-filtering (i.e., higher values for sigma perform less aggressive
%   splits).
%
% [1] Meyer, Fernand, "Topographic distance and watershed lines,” Signal 
% Processing , Vol. 38, July 1994, pp. 113-125.
%
% Copyright (C) 2018 Thomas Roetzer, MedUni Vienna
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%
    if (nargin == 1)
        gauss_sigma = 2;
    end
    if (nargin == 2)
        gauss_sigma = varargin{1};
    end
    
    EDM = -single(bwdist(~I));
    EDM = imgaussfilt(EDM, gauss_sigma);
    W = watershed(EDM);
    I(~W) = 0;
end

