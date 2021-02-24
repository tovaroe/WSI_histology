function [ annotation_image ] = extract_annotation( annotation_path, magnification )
%EXTRACT_ANNOTATION returns an annotation overlay image
%
%       annotation_image = extract_annotation(ANNOT_PATH, 
%       MAGNIFICATION) extracts FREEHAND annotations from an ndpa-file at 
%       ANNOT_PATH and retruns an annotation overlay image for the 
%       given MAGNIFICATION. The associated ndpi-WSI has to be present in 
%       the same folder and with the same name as the ndpa-annotation.
%
% Output:
% 1 ... RedMask
% 2 ... YellowMask
% 3 ... GreenMask
% 4 ... BlueMask
% 5 ... BlackMask
% 6 ... MagentaMask
% 7 ... TealMask
% 8 ... WhiteMask

% Copyright (C) 2021 Bernhard Baumann & Thomas Roetzer, MedUni Vienna
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

%% Main function

% Read associated WSI to get image size and pixel scaling
WSI_path = strrep(annotation_path, '.ndpa', '');
slidePtr = openslide_open(WSI_path);
[mppX, mppY, width, height, ~, ~, objectivePower] = openslide_get_slide_properties(slidePtr);
width = double(width*magnification/objectivePower);
height = double(height*magnification/objectivePower);

% Define functions for converting the ndpa-nanometer values to pixel
% values
XOffsetFromSlideCentre = openslide_get_property_value(slidePtr, 'hamamatsu.XOffsetFromSlideCentre');
YOffsetFromSlideCentre = openslide_get_property_value(slidePtr, 'hamamatsu.YOffsetFromSlideCentre');
nm2pixX = @(x) nm2pix(x, XOffsetFromSlideCentre, mppX, objectivePower, magnification, width);
nm2pixY = @(x) nm2pix(x, YOffsetFromSlideCentre, mppY, objectivePower, magnification, height);

openslide_close(slidePtr)
clear slidePtr

% Make annotation (output) image
annotation_image = zeros([height width],'uint8');

% Read number of annotations
NDPA = xml2struct(annotation_path);
NumAnnot = length(NDPA.annotations.ndpviewstate); % number of annotations

% Loop through annotations 
if NumAnnot
    for k = 1:NumAnnot

        % Only continue if annotation is freehand
        if NumAnnot==1
            ndpviewstate = NDPA.annotations.ndpviewstate;
        else
            ndpviewstate = NDPA.annotations.ndpviewstate{1,k};
        end

        if strcmp(ndpviewstate.annotation.Attributes.displayname, 'AnnotateFreehand')
            % Determine number of points for k-th annotation
            NumPoints = length(ndpviewstate.annotation.pointlist.point);

            if NumPoints > 1  % catch misannotations (single dots)

                % Allocate X and Y coordinates
                X = zeros(1,NumPoints);
                Y = X;

                % Extract X and Y coordinates from xml struct file
                for m = 1:NumPoints
                    X(m) = nm2pixX(str2double(ndpviewstate.annotation.pointlist.point{1,m}.x.Text));
                    Y(m) = nm2pixY(str2double(ndpviewstate.annotation.pointlist.point{1,m}.y.Text));
                end

                % Extract color of k-th annotation from xml struct file
                AnnotColor = ndpviewstate.annotation.Attributes.color;
                AnnotRGB = [hex2dec(AnnotColor(2:3)) hex2dec(AnnotColor(4:5)) hex2dec(AnnotColor(6:7))]./255;
                AnnotRGBstr = [ '[' num2str(AnnotRGB) ']' ];

                % Add annotation polygon to respective mask
                switch AnnotRGBstr
                    case '[1  0  0]' % red
                        annotation_image(poly2mask(X, Y, height,width)) = 1;                         
                    case '[1  1  0]' % yellow
                        annotation_image(poly2mask(X, Y, height,width)) = 2; 
                    case '[0  1  0]' % green
                        annotation_image(poly2mask(X, Y, height,width)) = 3; 
                    case '[0  0  1]' % blue
                        annotation_image(poly2mask(X, Y, height,width)) = 4; 
                    case '[0  0  0]' % black
                        annotation_image(poly2mask(X, Y, height,width)) = 5; 
                    case '[1  0  1]' % magenta
                        annotation_image(poly2mask(X, Y, height,width)) = 6; 
                    case '[0  1  1]' % teal
                        annotation_image(poly2mask(X, Y, height,width)) = 7;
                    case '[1  1  1]' % white
                        annotation_image(poly2mask(X, Y, height,width)) = 8;
                    otherwise
                        disp([AnnotRGBstr ' is an unknown annotation color'])
                end
            end
        end
    end
end

function pixel_coord = nm2pix(nm, OffsetFromSlideCenter, mpp, objectivePower, magnification, image_size)
    nm_coord = (nm - str2double(OffsetFromSlideCenter)); % Calculate coords in nm relative to image center
    pixel_coord = nm_coord/(mpp*1000); % Convert to coords in pixels
    pixel_coord = pixel_coord/(objectivePower/magnification); % Convert to current magnification
    pixel_coord = pixel_coord + round(image_size/2); % convert to absolute pixel_coords
    pixel_coord = round(pixel_coord);
end
end