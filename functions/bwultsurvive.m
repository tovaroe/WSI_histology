function output = bwultsurvive(bw, st, minNeighbors)
%BWULTSURVIVE returns a binary image
%
%    output = bwultsurvive(bw, st, minNeighbors) uses the bwsurvive
%    function to filter out positive pixels in bw that do not have at least
%    minNeighbors in the neighborhood given by the strel st repeatedly 
%    until arriving at a stable final binary image.
%
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
% Only pixels with minNeighbors in the given strel st survive
    overlap = false;
    while ~all(overlap(:))
        output = bwsurvive(bw, st, minNeighbors);
        overlap = output == bw; 
        bw = output;
    end
end