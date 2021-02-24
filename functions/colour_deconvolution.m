function image = colour_deconvolution(image, stains)
%COLOUR_DECONVOLUTION returns a 3-channel image
%
%   image = segment_cells(IM, stain)
%   performs colour deconvolution of the input image IM using the given
%   stain. Valid values for stain are 'HE', 'KB', 'LFB NFR', 'H DAB', 'H
%   DAB 2', and 'H DAB Ventana'.
%
% This method is based on Ruifrok's and Johnston's method (Ruifrok AC, 
% Johnston DA. Quantification of histochemical staining by color 
% deconvolution. Anal Quant Cytol Histol 23: 291-299, 2001.) and strongly
% inspired by the ImageJ implementation of Gabriel Landini
% (https://blog.bham.ac.uk/intellimic/g-landini-software/colour-deconvolution/).
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

if strcmp(stains, 'HE')
    % Color deconvolution of H&E stains
    % RGB OD matrix for hematoxylin, eosin and background
    v1 = [0.644; 0.72; 0.27];
    v2 = [0.07; 0.99; 0.11]; 
    v3 = [0.636; 0.001; 0.772];

elseif strcmp(stains, 'KB')
    v1 = [0.953; 0.303; -0.003];
    v2 = [0.004; 0.999; -0.006];
    v3 = [0.290; 0.001; 0.957];

elseif strcmp(stains, 'LFB NFR')
    v1 = [0.74;0.606;0.267];
    v2 = [0.214;0.851;0.478];
    v3 = [0.26800000;0.57000000;0.7760000];

elseif strcmp(stains, 'H DAB')
    v1 = [0.6500286; 0.704031; 0.2860126];
    v2 = [0.26814753; 0.57031375; 0.77642715];
    v3 = [0.7110272; 0.42318153; 0.5615672];

elseif strcmp(stains, 'H DAB 2')
    v1 = [0.7219949; 0.61885273; 0.30942637];
    v2 = [0.38015208; 0.58023214; 0.72028816];
    v3 = [0.57810706; 0.5294827; 0.62083834];

elseif strcmp(stains, 'H DAB Ventana')
    v1 = [0.7219949; 0.61885273; 0.30942637];
    v2 = [0.48616475; 0.6279628; 0.60770595];
    v3 = [0.49230802; 0.47189406; 0.7314019];
    
elseif strcmp(stains, 'HE Pigment')
    v1 = [0.56907135; 0.7671424; 0.29605788];
    v2 = [0.29937613; 0.8949957; 0.33069116];
    v3 = [0.41451415; 0.63914895; 0.6478168];    
end

% Calculate ODmatrix
ODmatrix = [v1/norm(v1) v2/norm(v2) v3/norm(v3)]';    

% Deconvolute the images
image = -log(im2single(image));
image = reshape(reshape(image, [], 3) / ODmatrix, size(image));
image = exp(-image);
image = im2uint8(image);

end

