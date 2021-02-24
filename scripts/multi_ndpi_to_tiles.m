%% Script: select multiple ndpi files for conversion into image-patches

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

%% Set parameters
magnification = 20;
tiles_size = [1024 1024];
overlap = [64 64];
splitHE = false;
makeAnnotationCSV = true;
minRelativeTissueArea = 0.9;

%% Convert ndpis to tiles
[FileNames,PathName,FilterIndex] = uigetfile('*.ndpi', 'Select HE WSIs!','MultiSelect', 'on');
output_path = uigetdir(pwd, 'Select output folder');
total = length(FileNames);

for k = 1:total
    
    FileName = (FileNames(k));
    FileName = char(FileName);
    
    foldername = strrep(strrep(strrep(FileName, '.ndpi', ''), ' ', '_'), '.', '_');
    outputdir = [output_path filesep foldername '_x' num2str(magnification)];
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    end
    
    ndpi_to_tiles([PathName filesep FileName], magnification, tiles_size, outputdir, ...
        'overlap', overlap, ...
        'splitHE', splitHE, ...
        'makeAnnotationCSV', makeAnnotationCSV, ...
        'minRelativeTissueArea', minRelativeTissueArea);
    
    disp(FileName);
end