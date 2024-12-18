% Copyright (C) 2019 Robert F Cooper, created 2018-11-07
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
clear;
close all force;

% Add our support libraries to the path.
basePath = which('DFT_Aggregator.m');
[basePath ] = fileparts(basePath);
path(path, fullfile(basePath,'lib'));

[fNames,thispath ]=uigetfile(fullfile(pwd,'*.mat'),'Select all montages you wish to combine.', 'MultiSelect', 'on');

if ~all(iscell(fNames))
    fNames = {fNames};
end

for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'scaling');
    scales(f) = scaling;   
end
scaling = unique(scales);

if length(scaling)>1
    error('The montage scales must be the same!');
end

global_eye = 'OD'; % Everything will be flipped to this.

foveal_coords = nan(length(fNames), 2);
montage_size = nan(length(fNames), 2);
montage_rect = cell(length(fNames), 1);

result_path = fullfile(thispath,'Results');
mkdir(result_path)


for f=1:length(fNames)
    disp(['Calculating Foveal Location: ' num2str(f) ' of ' num2str(length(fNames))])
    load( fullfile(thispath, fNames{f}), 'fovea_coords', 'imsize', 'unit');
    
    montage_size(f,:) = imsize(1:2);
    foveal_coords(f,:) = fovea_coords;
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the fovea coordinates in the x direction.
        foveal_coords(f,1) = montage_size(f,2)-foveal_coords(f,1);
    end
    
    montage_rect{f} = [               1                  1; % TLC
                       montage_size(f,2)                 1; % TRC
                       montage_size(f,2) montage_size(f,1); % BRC
                                      1  montage_size(f,1); % BLC
                                      1                  1];% TLC 
    
                                     
    % Shift all montage rectangles to a null point
    montage_rect{f} = montage_rect{f}-foveal_coords(f,:);    
end
clear fovea_coords
% Find the bounding rectangle
global_fovea_coords = min(cell2mat(montage_rect));
maxglobalbounds = max(cell2mat(montage_rect));

global_dimension = fliplr(ceil(maxglobalbounds-global_fovea_coords))+1; % Flipped so height/width (row/col), not width/height (x/y)

global_fovea_coords = round(-global_fovea_coords);

if ~exist('unit','var')
    liststr = {'microns (mm density)','degrees'};
    [unitind, oked] = listdlg('PromptString','Select output units:',...
                                          'SelectionMode','single',...
                                          'ListString',liststr);

    unit = liststr{unitind};
    if oked == 0
        error('Cancelled by user.');
    end
end

%% Add all of the dft information to the montage, and combine.

avg_density = zeros(global_dimension);
% weighted_avg_spacing = zeros(global_dimension);
avg_confidence = zeros(global_dimension);
combined_sum_map = zeros(global_dimension,'double');

% Find the weighted averages of all of the subjects.
for f=1:length(fNames)
    disp(['Calculating Averages: ' num2str(f) ' of ' num2str(length(fNames))])
    load( fullfile(thispath, fNames{f}), 'density_map_comb', 'blendederr_comb');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))+global_fovea_coords(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))+global_fovea_coords(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        density_map_comb = fliplr(density_map_comb);
        blendederr_comb = fliplr(blendederr_comb);
    end    

    % Lightly filter our data to remove any outliers in confidence or
    % density.
    [X, Y] = meshgrid(1:size(density_map_comb,2), 1:size(density_map_comb,1));

    if strcmp(unit, 'microns (mm density)')
        farrangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 800/scaling; 
    elseif strcmp(unit, 'degrees')
        farrangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 2.5/scaling;
    end
    
    confcombfar = blendederr_comb.*~farrangemask;
    confcombfar(confcombfar==0 | confcombfar>0.9)=NaN;
    lowconffar = quantile(confcombfar(~isnan(confcombfar) ), [0.1]);

    confcombclose = blendederr_comb.*farrangemask;
    confcombclose(confcombclose==0 | confcombclose>0.9)=NaN;
    lowconfclose = quantile(confcombclose(~isnan(confcombclose) ), [0.01]);

    if strcmp(unit, 'microns (mm density)')
        closerangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 400/scaling;
    elseif strcmp(unit, 'degrees')
        closerangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 1.25/scaling;
    end

    denscombclose = density_map_comb.*closerangemask;
    denscombclose(denscombclose==0)=NaN;
    densclose = quantile(denscombclose(~isnan(denscombclose)), [0.5]);

    % If we're looking at data with these filenames, drop their long
    % distance data as its quality is questionable.
    if contains(fNames{f},'11101') || contains(fNames{f},'11092')
        if strcmp(unit, 'microns (mm density)')
            good_data_mask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 2400/scaling; 
        elseif strcmp(unit, 'degrees')
            good_data_mask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 8/scaling;
        end

        density_map_comb = density_map_comb .*good_data_mask;
        blendederr_comb = blendederr_comb .*good_data_mask;
    end

    % Remove the offending values.
    density_map_comb(( (density_map_comb > densclose) & ~closerangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconffar) & ~farrangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconfclose) & farrangemask )) = NaN;
    %Also remove them from our confidence.
    blendederr_comb(isnan(density_map_comb)) = NaN;

     
    avg_density( rowrange, colrange) = sum(cat(3, avg_density( rowrange, colrange), density_map_comb),3,'omitnan');
    avg_confidence( rowrange, colrange) = sum(cat(3, avg_confidence( rowrange, colrange), blendederr_comb),3,'omitnan');
    combined_sum_map( rowrange, colrange) = sum(cat(3, combined_sum_map( rowrange, colrange), ~isnan(density_map_comb)),3,'omitnan');
    

    % figure(2); clf; hold on;
    % subplot(1,2,1); imagesc(density_map_comb); axis image; 
    % subplot(1,2,2); imagesc(combined_sum_map); axis image;
    % drawnow;

    clear density_map_comb blendederr_comb rowrange colrange
end

avg_density = avg_density./combined_sum_map;

avg_confidence = avg_confidence./combined_sum_map;

%% To find global foveal mask
foveal_density_map=avg_density(global_fovea_coords(2)-768:global_fovea_coords(2)+768, global_fovea_coords(1)-768:global_fovea_coords(1)+768);

foveal_density_map = imgaussfilt(foveal_density_map, 8);

polar_foveal_density_map = imcart2pseudopolar(foveal_density_map, 1, 1,[769, 769],'cubic');
polar_foveal_density_map(polar_foveal_density_map==0) = NaN;

ridgeline = nan(size(polar_foveal_density_map,1), 1);
ridgeloc = nan(size(polar_foveal_density_map,1), 1);

% figure; imagesc(polar_foveal_density_map); axis image;
% hold on; 
for p=1:size(polar_foveal_density_map,1)
    if strcmp(unit, 'microns (mm density)')
        [line, loc]=findpeaks(polar_foveal_density_map(p,:),'MinPeakProminence',4000,'SortStr', 'descend','NPeaks', 1);
    elseif strcmp(unit, 'degrees')
        [line, loc]=findpeaks(polar_foveal_density_map(p,:),'MinPeakProminence',750,'SortStr', 'descend','NPeaks', 1);
    end

    if ~isempty(line)
        ridgeline(p)=line;
        ridgeloc(p)=loc;
    end
    % plot(polar_foveal_density_map(p,:)); hold on;
    % plot(ridgeloc(p), p,'r*');
end
hold off;
min_ridge = round(min(ridgeline));


[x,y]=pol2cart((0:2*pi/360:2*pi-(2*pi/360))',ridgeloc);
x=x+769;
y=y+769;

y(isnan(x)) =[];
x(isnan(x)) =[];
% imagesc(foveal_density_map); hold on; plot(x,y)


smfoveamask = poly2mask(x,y,size(foveal_density_map,1),size(foveal_density_map,2));

% Set our tolerance just below our minimum ridge size, vs the foveal min
% centralvals = polar_foveal_density_map(:, floor(min(ridgeloc)/2));
% centralvals = mean(centralvals(:),'omitnan');
% foveal_density_map(isnan(foveal_density_map)) = 0;
% 
% tol = min_ridge-centralvals;
% smfoveamask = ~grayconnected(foveal_density_map, 769, 769, tol);
% 
% smfoveamask=imerode(smfoveamask,strel('disk',50));
% figure; imagesc(smfoveamask)
% 
threshold_mask = true(size(avg_density));

threshold_mask(global_fovea_coords(2)-768:global_fovea_coords(2)+768, global_fovea_coords(1)-768:global_fovea_coords(1)+768) = ~smfoveamask;

% Ensure that we only see the averages of data with more than 20% of our points going into it.
threshold_mask = logical(threshold_mask.*(combined_sum_map>=0.2*max(combined_sum_map(:)) )); 


%% Plot our masked data.

spacingcaxis = quantile(avg_density(threshold_mask(:)),[0.01 0.99]);
figure(1); imagesc(avg_density.*threshold_mask); title('Combined Spacing'); axis image;
clim(spacingcaxis); colorbar;

%% Output confidence image

avg_confidence(avg_confidence==0) = NaN;
lowconffar = quantile(avg_confidence(threshold_mask(:)),0.001);
rescaled_avg_conf = (avg_confidence.*threshold_mask)-lowconffar;
rescaled_avg_conf(rescaled_avg_conf<0) = 0;

highconf = quantile(rescaled_avg_conf(rescaled_avg_conf~=0),0.995); 
rescaled_avg_conf = 255*(rescaled_avg_conf./highconf); 
rescaled_avg_conf(rescaled_avg_conf>255) = 255;

turbfull=flipud(turbo(512));

imwrite(rescaled_avg_conf, turbfull(19:274, :), fullfile(result_path,[num2str(length(fNames)) '_subjects_combined_conf.tif']))

figure(2); clf; imagesc((avg_confidence.*threshold_mask)); axis image; colormap(turbfull(19:274, :)); colorbar; title('Combined Confidence');
clim([lowconffar lowconffar+highconf]); colorbar;
saveas(gcf,fullfile(result_path,[num2str(length(fNames)) '_subjects_combined_conf.svg']));
clear rescaled_avg_conf;
%% Output density image

avg_density(avg_density==0) = NaN;
lowdens = quantile(avg_density(threshold_mask(:)),0.01)
rescaled_avg_density = (avg_density.*threshold_mask)-lowdens;
rescaled_avg_density(rescaled_avg_density<0) = 0;
highdens = quantile(rescaled_avg_density(rescaled_avg_density~=0),0.995)
rescaled_avg_density = 255*(rescaled_avg_density./highdens); 
rescaled_avg_density(rescaled_avg_density>255) = 255;

figure(3); imagesc(avg_density.*threshold_mask); title('Combined Density'); axis image;
clim([lowdens lowdens+highdens]); colorbar;
imwrite(rescaled_avg_density, parula(256), fullfile(result_path,[num2str(length(fNames)) '_subjects_combined_density.tif']))
saveas(gcf,fullfile(result_path,[num2str(length(fNames)) '_subjects_combined_density.svg']));
clear rescaled_avg_density;

%% Output sum map

rescaled_sum_map = uint8( (255*(combined_sum_map./max(combined_sum_map(:)))).*threshold_mask );

figure(4); imagesc(combined_sum_map); title('Sum Map'); axis image; colormap winter; colorbar;
clim([0 max(combined_sum_map(:))]);

saveas(gcf,fullfile(result_path,[num2str(length(fNames)) '_subjects_sum_map.svg']));
imwrite(rescaled_sum_map, winter(256), fullfile(result_path,[num2str(length(fNames)) '_subjects_sum_map.tif']));

%% Determine average/stddev of all data.
value_map_variance = nan(global_dimension);

% For finding outliers:
% figure(1); clf; hold on;
    

% Find the std deviations of all of the subjects.
for f=1:length(fNames)
    disp(['Calculating Standard Deviation: ' num2str(f) ' of ' num2str(length(fNames))])
    load( fullfile(thispath, fNames{f}), 'density_map_comb', 'blendederr_comb');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))+global_fovea_coords(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))+global_fovea_coords(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        density_map_comb = fliplr(density_map_comb);
        blendederr_comb = fliplr(blendederr_comb);
    end
    
    % Lightly filter our data to remove any outliers in confidence or
    % density.
    [X, Y] = meshgrid(1:size(density_map_comb,2), 1:size(density_map_comb,1));

    if strcmp(unit, 'microns (mm density)')
        farrangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 800/scaling; 
    elseif strcmp(unit, 'degrees')
        farrangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 2.5/scaling;
    end
    
    confcombfar = blendederr_comb.*~farrangemask;
    confcombfar(confcombfar==0)=NaN;
    lowconffar = quantile(confcombfar(~isnan(confcombfar)), [0.10])

    confcombclose = blendederr_comb.*farrangemask;
    confcombclose(confcombclose==0)=NaN;
    lowconfclose = quantile(confcombclose(~isnan(confcombclose)), [0.01])

    if strcmp(unit, 'microns (mm density)')
        closerangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 400/scaling;
    elseif strcmp(unit, 'degrees')
        closerangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 1.25/scaling;
    end

    denscombclose = density_map_comb.*closerangemask;
    denscombclose(denscombclose==0)=NaN;
    densclose = quantile(denscombclose(~isnan(denscombclose)), [0.5])

    % If we're looking at data with these filenames, drop their long
    % distance data as its quality is questionable.
    if contains(fNames{f},'11101') || contains(fNames{f},'11092')
        if strcmp(unit, 'microns (mm density)')
            good_data_mask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 2400/scaling; 
        elseif strcmp(unit, 'degrees')
            good_data_mask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 8/scaling;
        end

        density_map_comb = density_map_comb .*good_data_mask;
        blendederr_comb = blendederr_comb .*good_data_mask;
    end

    % Remove the offending values.
    density_map_comb(( (density_map_comb > densclose) & ~closerangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconffar) & ~farrangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconfclose) & farrangemask )) = NaN;
    %Also remove them from our confidence.
    blendederr_comb(isnan(density_map_comb)) = NaN;
    
    density_map_comb = density_map_comb.*threshold_mask( rowrange, colrange);
    density_map_comb(density_map_comb==0) = NaN;

    thisdiff = sum( cat(3, density_map_comb, -avg_density( rowrange, colrange)) ,3).^2; % This line is good. Nans need to be carried through.
    

    value_map_variance( rowrange, colrange) = sum( cat(3,value_map_variance( rowrange, colrange), thisdiff), 3,'omitnan'); 
     clear thisdiff;


    clear blendedim rowrange colrange
end
save( fullfile(result_path,[num2str(length(fNames)) '_aggregate_data.mat']), ...
                            'avg_density', 'avg_confidence','global_fovea_coords',...
                            'global_eye', 'threshold_mask', 'value_map_variance',...
                            'unit', 'scaling', 'montage_rect','combined_sum_map','-v7.3');

%%
value_std_dev = value_map_variance;
value_std_dev = value_std_dev./(combined_sum_map-1);
value_std_dev(isinf(value_std_dev)) =0;
value_std_dev = real(sqrt(value_std_dev));

lowdens = quantile(value_std_dev(threshold_mask(:)),0.01)
rescaled = (value_std_dev.*threshold_mask)-lowdens;
rescaled(rescaled<0) = 0;
highdens = quantile(rescaled(rescaled~=0),0.995)
rescaled = 255*(rescaled./highdens); 
rescaled(rescaled>255) = 255;
imwrite(rescaled, parula(256), fullfile(result_path,[num2str(length(fNames)) '_subjects_combined_stddev.tif']))
 % maskedspac = value_std_dev.*threshold_mask;

figure(3); imagesc(value_std_dev.*threshold_mask); title('Combined Spacing Std dev');
saveas(gcf,[num2str(length(fNames)) 'subjects_combined_stddev.png']);
saveas(gcf,[num2str(length(fNames)) 'subjects_combined_stddev.svg']);
