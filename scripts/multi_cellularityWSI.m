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

%%
% select multiple images for segmentation
[FileNames,PathName,FilterIndex] = uigetfile('*.ndpi', 'Select HE WSIs!','MultiSelect', 'on');

if ~exist([PathName filesep 'Analyses'], 'dir')
    mkdir([PathName filesep 'Analyses']);
end

% Sort out already analyzed files
if exist([PathName filesep 'HE_analysis.csv'], 'file')
    analyses = readtable([PathName filesep 'HE_analysis.csv']);
    FilesAnalyzed = analyses.id;
    FileNames = setdiff(FileNames, intersect(FileNames, FilesAnalyzed));
end

for k = 1:length(FileNames)
    disp(['Roughly ' num2str(100*k/length(FileNames)) '%'])
    FileName = (FileNames(k));
    FileName = char(FileName);
    cellularityWSI(FileName,PathName)
end
