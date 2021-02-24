function I_t = phansalkar(I, varargin)
%PHANSALKAR returns a thresholded (binary) image
%
%    I_t = phansalkar(IM, VARARGIN) performs phansalkar thresholding [1] of 
%    the input image IM and returns the thresholded (binary) image.
%
%OPTIONAL PARAMETERS
% are the parameters radius, k, r, p, q (given in this exact order as input 
% arguments) as described in the ImageJ/Fiji implementation (see 
% https://imagej.net/Auto_Local_Threshold).
%
% [1] Based on Phansalkar, Neerad, et al. "Adaptive local thresholding for 
% detection of nuclei in diversity stained cytology images." Communications
% and Signal Processing (ICCSP), 2011 International Conference on. IEEE, 
% 2011.
% Strongly inspired by the ImageJ implementation.

% Copyright (C) 2018 Thomas Roetzer, MedUni Vienna
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
%%
I = im2double(I);
switch nargin
    case 1
        radius = 15;
        k = 0.25;
        r = 0.5;
        p = 2;
        q = 10;
    case 2
        radius = varargin{1};
        k = 0.25;
        r = 0.5;
        p = 2;
        q = 10;
    case 3
        radius = varargin{1};
        k = varargin{2};
        r = 0.5;
        p = 2;
        q = 10;
    case 4
        radius = varargin{1};
        k = varargin{2};
        r = varargin{3};
        p = 2;
        q = 10;
    case 5
        radius = varargin{1};
        k = varargin{2};
        r = varargin{3};
        p = varargin{4};
        q = 10;
    case 6
        radius = varargin{1};
        k = varargin{2};
        r = varargin{3};
        p = varargin{4};
        q = varargin{5};
    otherwise
        error('Invalid number of input arguments')
end

% get a disk-neighbourhood
SE = strel('disk', radius, 0);
disk = SE.Neighborhood;
disk_sum = sum(disk(:));
disk_norm = disk/disk_sum;

I_pad = padarray(I, [radius radius], 'symmetric');
means = conv2(I_pad, disk_norm, 'valid');
stds  = stdfilt(I, disk);

thresholds = means .* (1 + p * exp(-q * means) + k * ((stds / r) - 1));

I_t = (I >= thresholds);

end

