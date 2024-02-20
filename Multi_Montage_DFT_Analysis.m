% Created 1/13/21
%
% This script's sole job is to run multiple montages in series; to save time, 
% instead of operating on individual files though, it operates on a series
% of folders.
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

clear;
close all;

% Add our support libraries to the path.
basePath = which('Multi_Montage_DFT_Analysis.m');
[basePath ] = fileparts(basePath);
path(path, fullfile(basePath,'lib')); 


liststr = {'microns (mm density)','degrees'};
[selectedunit, oked] = listdlg('PromptString','Select output units:',...
                                      'SelectionMode','single',...
                                      'ListString',liststr);
if oked == 0
    error('Cancelled by user.');
end

unit = liststr{selectedunit};

thisfolder = pwd;

thisfolder = uigetdir(thisfolder, 'Select the folder containing the montage folders you wish to assess.');
    
[folderList, isdir]=read_folder_contents(thisfolder);

[lutfname, lutfolder] = uigetfile(fullfile(pwd,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

[fovfname, fovfolder] = uigetfile(fullfile(pwd,'*.csv'),'Select fovea image list, OR cancel if you do not want to analyze foveas separately.');
if fovfname ~=0    
    fovea_filenames = readcell(fullfile(fovfolder, fovfname), 'Delimiter',',');
end

%%
restartf = 1;
endf = 1; %length(folderList);
%%
for f=restartf:endf
    restartf=f;
    confocal_coords=[];
    something=false;
    if isdir{f}
        confocaldir=fullfile(thisfolder, folderList{f},'confocal');
        splitdir=fullfile(thisfolder, folderList{f},'split detection');
        
        if exist(confocaldir, 'dir')
            disp(['***** Running analysis on confocal data: ' confocaldir ' *****']);
            fNames = read_folder_contents(confocaldir,'tif');

            [ scalinginfo, unit, lut ]=determine_scaling(confocaldir, fNames, fullfile(lutfolder,lutfname) ,unit);


            % Run an analysis on the selected foveal image, if the file
            % exists.
            % splitlist=split(filelist{f,1},'_');
            % filedir=fullfile(lutfolder, [splitlist{2} '_' splitlist{4}],'confocal');

            
            if ~isempty(fovea_filenames)                
                fov_ind = find(cellfun(@any, cellfun(@(s) strcmp(s, fovea_filenames), fNames,'UniformOutput', false)), 1, 'first');
                if isempty(fov_ind)
                    error('Unable to match foveal image to current dataset. Ensure EACH dataset has a designated foveal image or restart without selecting a fovea image file.')
                end
                disp(['***** Running analysis on foveal image: ' fNames{fov_ind} ' *****']);

                Fovea_DFT_Analysis(confocaldir, fNames{fov_ind}, scalinginfo, unit, lut, true);
                
            end

            disp(['***** Running analysis on confocal montage: ' confocaldir ' *****']);
            [confocal_coords, mask]=Foveated_Montage_DFT_Analysis(confocaldir, fNames, scalinginfo, unit, lut, true);
            something=true;
        end
        close all;
        
        if exist(splitdir, 'dir')
            disp(['***** Running analysis on split-detector montage: ' splitdir ' *****']);
            fNames = read_folder_contents(splitdir,'tif');

            [ scalinginfo, ~, lut ]=determine_scaling(splitdir, fNames, fullfile(lutfolder,lutfname) ,unit);

            Foveated_Montage_DFT_Analysis(splitdir, fNames, scalinginfo, unit, lut, true, confocal_coords, mask);
            something=true;
        end
                
        
        if ~something
            disp(['***** ' folderList{f} ' doesn''t contain any analyzable materials. *****']);
        end
    end
end

%% Process each of the above, merging their pieces together
restartf=1;
endf=length(folderList);
%%
for f=restartf:endf
    restartf=f;
    confocal_coords=[];
    something=false;
    if isdir{f}
        confocaldir=fullfile(thisfolder, folderList{f},'confocal','Results_foveated');
        foveadir=fullfile(thisfolder, folderList{f},'confocal','Results_fovea_map');
        splitdir=fullfile(thisfolder, folderList{f},'split detection','Results_foveated');
        
        if exist(confocaldir, 'dir') && exist(splitdir, 'dir')
            fNameC = read_folder_contents(confocaldir,'mat');
            fNameCF = read_folder_contents(foveadir,'mat');
            fNameS = read_folder_contents(splitdir,'mat');
            
            if ~isempty(fNameC) && ~isempty(fNameS)
                something=true;
                
                
                if fovfname ~=0   
                    disp(['***** Merging data from: ' fNameCF{1} ', ' confocaldir ' and\n ' splitdir ' *****']);
                    load(fullfile(foveadir, fNameCF{1}), 'densim', 'imbox', 'confmap', 'summap');
                    
                    densim = densim(:,:,1);
                    densim(isnan(densim)) = 0;
                    blendederrim_fov = zeros(size(densim));
                    blendederrim_fov(imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3)) = confmap./summap;
                    blendederrim_fov(isnan(blendederrim_fov)) = 0;
    
                    foveal_data_mask = blendederrim_fov ~= 0;
                    fovea_blend_range = 64;
                    fovea_merge_loc = floor(min(size(confmap))/2); % Set the max merge loc based on the short edge of the fovea image.
                    density_map_fov = densim;
                    clear densim imbox confmap summap
                else
                    disp(['***** Merging data from: ' confocaldir ' and\n ' splitdir ' *****']);
                end
                
                load(fullfile(confocaldir, fNameC{1}), 'density_map', 'blendederrim','fovea_coords', 'unit', 'scaling');
                blendederrim_conf = blendederrim;
                density_map_conf = density_map;
                confocalmask = ~isnan(density_map_conf);
                
                load(fullfile(splitdir, fNameS{1}), 'density_map', 'blendederrim');
                blendederrim_split = blendederrim.*confocalmask;
                density_map_split = density_map.*confocalmask;
                
                disp('Data loaded, merging...');
                clear density_map blendederrim
            
                errfovpolar = imcart2pseudopolar(blendederrim_fov, 1, .5, fovea_coords,'makima' , 0);
                errconfpolar = imcart2pseudopolar(blendederrim_conf, 1, .5, fovea_coords,'makima' , 0);
                errsplitpolar = imcart2pseudopolar(blendederrim_split, 1, .5, fovea_coords,'makima' , 0);
                errfovpolar(errfovpolar<=0) = NaN;
                errconfpolar(errconfpolar<=0) = NaN;
                errsplitpolar(errsplitpolar<=0) = NaN;
                
                
                fovannuli = zeros(size(blendederrim_fov));
                confannuli = zeros(size(blendederrim_conf));
                splitannuli = zeros(size(blendederrim_conf));

                avgdifferr = (mean(errconfpolar,'omitnan')-mean(errsplitpolar,'omitnan'))./ ...
                             ( (mean(errsplitpolar,'omitnan') + mean(errconfpolar,'omitnan'))/2 );

                % figure(5); plot(mean(errconfpolar,'omitnan')); hold on; plot(mean(errsplitpolar,'omitnan')); plot(mean(errfovpolar,'omitnan')); plot(avgdifferr);  hold off; 
                % axis([0 length(errconfpolar) -2 2]); legend('Confocal Confidence', 'Split Confidence', 'foveal', '% diff');
                % drawnow;
                % saveas(gcf, fullfile(thisfolder, folderList{f},'conf_split_polarerror.png') );
                % 
                
                offset = find((avgdifferr>=0) == 1, 1, 'first');    
                
                MAXOFFSET = 2000;
                
                % We want to make sure that we found a valid blending spot
                % for the confocal to split.
                if ~isempty(offset) && offset < MAXOFFSET
                    
                    confhigh  = find((avgdifferr(:, offset:end)>=-0.2) == 0, 1, 'first') + offset;
                    splithigh = find( (avgdifferr(:, offset:end)<=0.2) == 1, 1, 'first') + offset;

                    if isempty(confhigh)
                        confhigh = find((avgdifferr(:, offset:end)>=0) == 0, 1, 'first') + offset;
                    end
                    
                    if isempty(splithigh)
                        splithigh = find((avgdifferr(:, offset:end)<=0) == 0, 1, 'first') + offset;
                    end
                    
                    blendrange = round((confhigh-splithigh)/2);
                                                    
                    mergeloc = confhigh-blendrange;

                else
                    disp('Warning: unable to find ideal merging location. Guessing from closest values...');
                    [val, offset] = max(avgdifferr(1:MAXOFFSET));
                    confhigh  = offset; %this needs work.
                    mergeloc = confhigh;
                    blendrange = 512;
                end

                % First we merge in our fovea to the confocal map, if it is
                % available.
                if fovfname ~=0 
                    foveadisk = strel('disk',fovea_merge_loc,0);
                    foveadisk = foveadisk.Neighborhood;
                    
                    diskshiftx = floor(fovea_coords(1)-(size(foveadisk,2)/2));
                    disksizex = diskshiftx+size(foveadisk,2)-1;
                    diskshifty = floor(fovea_coords(2)-(size(foveadisk,1)/2));
                    disksizey = diskshifty+size(foveadisk,1)-1;
                    
                    diskstartx = 1;
                    diskstarty = 1;
                    diffx=0;    
                    diffy=0;
                    
                    if diskshiftx<1
                        diskstartx = (-diskshiftx)+1;
                        diskshiftx=1;                    
                    end
                    if diskshifty<1
                        diskstarty = (-diskshifty)+1;
                        diskshifty=1;                    
                    end
                    
                    if disksizex > size(fovannuli,2)
                        diffx = size(fovannuli,2)-disksizex;
                        disksizex = disksizex+diffx;
                    end
                    if disksizey > size(confannuli,1)
                        diffy = size(fovannuli,1)-disksizey;
                        disksizey = disksizey+diffy;
                    end
                    
                    fovannuli(diskshifty:disksizey,...
                               diskshiftx:disksizex) = foveadisk(diskstarty:end+diffy,...
                                                                diskstartx:end+diffx);
                    %Trim the annulus based on the actual presence of foveal
                    %data, so we don't get weird behavior. Erode it by the
                    %blend range so that we have data to blend with (no hard
                    %edges)
                    fovannuli = imerode(foveal_data_mask, ones(2*fovea_blend_range)).*fovannuli;
                    
                    conf_to_fovannuli = abs(1-fovannuli);
                    
    %                 figure(10); subplot(1,2,1); imagesc(fovannuli); axis image;
                    
                    fovannuli = imgaussfilt(fovannuli, fovea_blend_range/2);
                    conf_to_fovannuli = imgaussfilt(conf_to_fovannuli, fovea_blend_range/2);
                    
    %                 figure(10); subplot(1,2,2); imagesc(fovannuli); axis image;
                    
                    density_map_fov = fovannuli.*density_map_fov;
                    fovdensity_map_conf = conf_to_fovannuli.*density_map_conf;
                    
                    density_map_conf = (density_map_fov+fovdensity_map_conf)./(fovannuli+conf_to_fovannuli);
                end

                
                % THEN we merge the confocal to the split.
                confdisk = strel('disk',mergeloc,0);
                confdisk = confdisk.Neighborhood;

                diskshiftx = floor(fovea_coords(1)-(size(confdisk,2)/2));
                disksizex = diskshiftx+size(confdisk,2)-1;
                diskshifty = floor(fovea_coords(2)-(size(confdisk,1)/2));
                disksizey = diskshifty+size(confdisk,1)-1;
                
                diskstartx = 1;
                diskstarty = 1;
                diffx=0;    
                diffy=0;
                
                if diskshiftx<1
                    diskstartx = (-diskshiftx)+1;
                    diskshiftx=1;                    
                end
                if diskshifty<1
                    diskstarty = (-diskshifty)+1;
                    diskshifty=1;                    
                end
                
                if disksizex > size(confannuli,2)
                    diffx = size(confannuli,2)-disksizex;
                    disksizex = disksizex+diffx;
                end
                if disksizey > size(confannuli,1)
                    diffy = size(confannuli,1)-disksizey;
                    disksizey = disksizey+diffy;
                end

                confannuli(diskshifty:disksizey,...
                           diskshiftx:disksizex) = confdisk(diskstarty:end+diffy,...
                                                            diskstartx:end+diffx);
                       
                splitannuli = abs(1-confannuli);
                
                % Blur over the blend range
                confannuli = imgaussfilt(confannuli, blendrange/2);
                splitannuli = imgaussfilt(splitannuli, blendrange/2);
                
                % Weight each map against its filtered annuli
                density_map_conf = confannuli.*density_map_conf;
                density_map_split = splitannuli.*density_map_split;

                
                % figure(2); subplot(1,2,1); imagesc(density_map_split); axis image;
                % subplot(1,2,2); imagesc(density_map_conf);axis image;
                % drawnow;
                % saveas(gcf, fullfile(thisfolder, folderList{f},'conf_split_density_areas.png') );
                % close;
                
                density_map_comb = (density_map_conf+density_map_split)./(confannuli+splitannuli);
                
                clear density_map_conf density_map_split


                blendederrim_conf = confannuli.*blendederrim_conf;
                blendederrim_split = splitannuli.*blendederrim_split;
                
                blendederr_comb = (blendederrim_conf+blendederrim_split)./(confannuli+splitannuli);
                
                clear blendederrim_split blendederrim_conf 
                
                
                % figure(3); imagesc(density_map_comb); axis image;
                % drawnow;
                % saveas(gcf, fullfile(thisfolder, folderList{f},'merged_density.png') );
                % figure(4);
                % imagesc(blendederr_comb); axis image;
                % drawnow;
                % saveas(gcf, fullfile(thisfolder, folderList{f},'merged_error.png') );
                
                
                safesave_sm(fullfile(confocaldir, fNameC{1}), fullfile(thisfolder, 'Aggregation_Analysis', strrep(fNameC{1}, 'confocal', 'merged')), density_map_comb, blendederr_comb);
                clear density_map_comb confannuli splitannuli
            end
        end
        

        if ~something
            disp(['***** ' folderList{f} ' doesn''t contain any analyzable materials. *****']);
        end
    end
end

function []= safesave_sm(baseload, newsave, denscomb, errcomb)
    disp(['Using path: ' baseload ' as base for ' newsave]);
    load(baseload);
    density_map_comb = denscomb;
    blendederr_comb = errcomb;
    save(newsave, 'density_map_comb', 'blendederr_comb', 'fovea_coords', 'imsize', 'scaling',  '-v7.3')
end

function []= safesave(baseload, newsave, denscomb, errcomb, confannuli, splitannuli)
    disp(['Using path: ' baseload ' as base for ' newsave]);
    load(baseload);
    density_map_comb = denscomb;
    blendederr_comb = errcomb;
    confocal_annulus = confannuli;
    split_annulus = splitannuli;
    save(newsave, '-v7.3')
end

