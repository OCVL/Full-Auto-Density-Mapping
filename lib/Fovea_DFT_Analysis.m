function [] = Fovea_DFT_Analysis(fovea_path, fovea_file, scaling, unit, lut, smallestscale)
% function [] = Fovea_DFT_Analysis(fovea_path, fovea_file, scaling, unit, lut, smallestscale)
%
%
% Created 2/20/24
%
% Copyright (C) 2024 Robert F Cooper
% 
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

    if ~exist('smallestscale','var')
       smallestscale=false; 
    end
    
    
    if ~isempty(lut) % If we didn't directly input a scale,
        % ask if we want to scale the montage to a common (smallest) scale from given LUT.
        pressedbutton=[];
        if ~exist('smallestscale','var')
            pressedbutton = questdlg('Scale montage to the common (smallest) scale from the selected LUT?',...
                                     'Scale to common size?', 'No');
        end
        
        if strcmp(pressedbutton,'Yes') || smallestscale
    
            all_scales = nan(length(lut{1}),1);
    
            for i=1:length(lut{1})
                % Calculate the scale for each identifier.                                
    
                axiallength = lut{2}(i);
                pixelsperdegree = lut{3}(i);
    
                micronsperdegree = (291*axiallength)/24;
    
                switch unit
                    case 'microns (mm density)'
                        all_scales(i) = 1 / (pixelsperdegree / micronsperdegree);
                    case 'degrees'
                        all_scales(i) = 1/pixelsperdegree;
                    case 'arcmin'
                        all_scales(i) = 60/pixelsperdegree;
                end
            end
    
            global_scale = min(all_scales); % Everything will be scaled to this.
    
            rel_scale = scaling./global_scale; % We will need to scale our images by this.
            scaling = global_scale; % The new scale becomes this.
    
        elseif strcmp(pressedbutton,'Cancel')
            return
        end
    end
    
    foveaim = imread(fullfile(fovea_path, fovea_file));
    imsize = size(foveaim);
    imsize = round(imsize(1:2).*rel_scale);
    
    if contains(fovea_file,'11102_OD')  % This factor was added to account for trial lenses.
        imsize = round( imsize(1:2) * 1.13 );
    elseif contains(fovea_file,'11051_OD')  % This factor was added to account for trial lenses.
        imsize = round( imsize(1:2) * 1.1369 );
    end
    
    foveaim = double(imresize(foveaim, imsize(1:2)));
            
    [~, spacmap, confmap, summap, imbox] = fit_fourier_spacing(foveaim, 96);
            
    imsize = size(foveaim);
    densim = zeros(imsize(1:2));
    
    weightedmap = spacmap./confmap;
    
    density_map = sqrt(3)./ (2*(weightedmap*scaling).^2);
    density_map(isinf(density_map))=0;
    
    if strcmp(unit,'microns (mm density)')
        density_map = (1000^2).*density_map;
    end
    
    densim(imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3)) = density_map;
    
    result_path = fullfile(fovea_path,'Results_fovea_map');
    mkdir(result_path)
           
    lower01 = quantile(density_map(~isnan(density_map)),0.01);
    
    scaled_density = density_map-lower01;
    upper99 = quantile(scaled_density(~isnan(scaled_density)),0.99);
    
    scaled_density = 255.*(scaled_density./ upper99 );
    scaled_density(scaled_density>255) = 255;
        
    fnamesplits = strsplit(fovea_file,'_');
    prefix=[];
    if length(fnamesplits)<4
        numprefixes = length(fnamesplits);
    else
        numprefixes = 4;
    end
    
    for f=1:numprefixes % build our prefix, using the first 4 things separated by underscores, or the full filename, whichever comes first.
        prefix = [prefix fnamesplits{f} '_'];
    end

    imwrite(uint8(scaled_density),parula(256), fullfile(result_path, [prefix 'Fovea.tif']));
    safesave(fullfile(result_path, [prefix 'Fovea.mat']),...
             weightedmap,density_map,densim, spacmap, confmap, summap, imbox);
            
    
end

function []= safesave(newsave, x1, x2, x3, x4, x5, x6, x7)
    % disp(['Saving to:' newsave]);
    weightedmap = x1;
    density_map = x2;
    densim = x3;
    spacmap = x4;
    confmap = x5;
    summap = x6;
    imbox = x7;

    save(newsave, 'weightedmap','density_map','densim', 'spacmap', 'confmap', 'summap', 'imbox', '-v7.3')
end