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
    load( fullfile(thispath, fNames{f}), 'fovea_coords', 'imsize');
    
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

    farrangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 2.5/scaling;   
    
    confcombfar = blendederr_comb.*~farrangemask;
    confcombfar(confcombfar==0)=NaN;
    lowconffar = quantile(confcombfar(~isnan(confcombfar)), [0.10])

    confcombclose = blendederr_comb.*farrangemask;
    confcombclose(confcombclose==0)=NaN;
    lowconfclose = quantile(confcombclose(~isnan(confcombclose)), [0.01])

    closerangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 1.25/scaling;
    denscombclose = density_map_comb.*closerangemask;
    denscombclose(denscombclose==0)=NaN;
    densclose = quantile(denscombclose(~isnan(denscombclose)), [0.5])

    % If we're looking at data with these filenames, drop their long
    % distance data as its quality is questionable.
    if contains(fNames{f},'11101') || contains(fNames{f},'11092')
        good_data_mask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 8/scaling;
        density_map_comb = density_map_comb .*good_data_mask;
        blendederr_comb = blendederr_comb .*good_data_mask;
    end

    % Remove the offending values.
    density_map_comb(( (density_map_comb > densclose) & ~closerangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconffar) & ~farrangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconfclose) & farrangemask )) = NaN;
    %Also remove them from our confidence.
    blendederr_comb(isnan(density_map_comb)) = NaN;

    figure(1); clf; hold on;
    subplot(1,2,1); imagesc(density_map_comb); axis image; 
    subplot(1,2,2); imagesc(blendederr_comb); axis image; 
    title( strrep(fNames{f},'_','\_') ); drawnow; hold off;

    

    avg_density( rowrange, colrange) = sum(cat(3, avg_density( rowrange, colrange), density_map_comb),3,'omitnan');
    avg_confidence( rowrange, colrange) = sum(cat(3, avg_confidence( rowrange, colrange), blendederr_comb),3,'omitnan');
    combined_sum_map( rowrange, colrange) = sum(cat(3, combined_sum_map( rowrange, colrange), ~isnan(density_map_comb)),3,'omitnan');
    clear density_map_comb blendederr_comb rowrange colrange
end

avg_density = avg_density./combined_sum_map;
% avg_spacing= avg_spacing.*2/sqrt(3); % To ICD <-- only needed if row, not cell, spacing.
avg_confidence = avg_confidence./combined_sum_map;

%% To find foveal mask
% avg_density = avg_spacing; %.*sqrt(3)/2; %<-- only needed if row, not cell, spacing.
% avg_density = (1000^2).*sqrt(3)./ (2*(avg_density).^2); %<-- for cells/mm
% avg_density = sqrt(3)./ (2*(avg_spacing).^2); %<-- for cells/deg

foveal_density_map=avg_density(global_fovea_coords(2)-768:global_fovea_coords(2)+768, global_fovea_coords(1)-768:global_fovea_coords(1)+768);

polar_foveal_density_map = imcart2pseudopolar(foveal_density_map, 1, 1,[769, 769],'makima');

ridgeline = nan(size(polar_foveal_density_map,1), 1);

figure; imagesc(polar_foveal_density_map); axis image; hold on; 
for p=1:size(polar_foveal_density_map,1)
    [ridgeline(p), ridgeloc] = findpeaks(polar_foveal_density_map(p,:),'MinPeakProminence',750,'NPeaks', 1);
%     plot(polar_foveal_density_map(p,:)); hold on;
    plot(ridgeloc, p,'r*');
end
hold off;
min_ridge = round(min(ridgeline));

% Set our tolerance just below our minimum ridge size
tol = min_ridge-foveal_density_map(769,769)-100;
smfoveamask = ~grayconnected(foveal_density_map, 769, 769, tol);

smfoveamask=imerode(smfoveamask,strel('disk',50));
figure; imagesc(smfoveamask)


threshold_mask = true(size(avg_density));

threshold_mask(global_fovea_coords(2)-768:global_fovea_coords(2)+768, global_fovea_coords(1)-768:global_fovea_coords(1)+768) = smfoveamask;


% Ensure that we only see the averages of data with more than X points going into it.
threshold_mask = logical(threshold_mask.*(combined_sum_map>=10)); 


%% Plot our masked data.

spacingcaxis = quantile(avg_density(threshold_mask(:)),[0.01 0.99]);
figure(1); imagesc(avg_density.*threshold_mask); title('Combined Spacing'); axis image;
caxis(spacingcaxis); colorbar;

%% Output spacing image
rescaled_avg_spacing = (avg_density.*threshold_mask)-spacingcaxis(1);
rescaled_avg_spacing(rescaled_avg_spacing<0) = 0;
rescaled_avg_spacing = 255*(rescaled_avg_spacing./quantile(rescaled_avg_spacing(rescaled_avg_spacing~=0),0.99)); 
rescaled_avg_spacing(rescaled_avg_spacing>255) = 255;

imwrite(rescaled_avg_spacing, parula(256), fullfile(result_path, [num2str(length(fNames)) 'subjects_combined_spacing.tif']))
clear rescaled_avg_spacing;

%% Output error image
% errquartiles = quantile(avg_error(avg_error~=0 & ~isnan(avg_error) & ~isinf(avg_error)),[0.001 0.01 0.05]);
% [firemap, amap] = firecmap(errquartiles(1), errquartiles(3), errquartiles(2),...
%                            min(avg_error(avg_error~=0 & ~isnan(avg_error)& ~isinf(avg_error))), ...
%                            max(avg_error(avg_error~=0 & ~isnan(avg_error)& ~isinf(avg_error))), 256);
% 
% rescaled_avg_err = avg_error.*threshold_mask;
% rescaled_avg_err = 255*rescaled_avg_err; 
% rescaled_avg_err(rescaled_avg_err>255) = 255;
% imwrite(rescaled_avg_err, firemap, [num2str(length(fNames)) 'subjects_combined_error.png'], 'Transparency',amap)
% 
% figure(2); imagesc(rescaled_avg_err); colormap(firemap); axis image; colorbar; title('Average Error');
% clear rescaled_avg_err;
%% Output confidence image
avg_confidence(avg_confidence==0) = NaN;
lowconffar = quantile(avg_confidence(threshold_mask(:)),0.001);
rescaled_avg_conf = (avg_confidence.*threshold_mask)-lowconffar;
rescaled_avg_conf(rescaled_avg_conf<0) = 0;
highconf = quantile(rescaled_avg_conf(rescaled_avg_conf~=0),0.99); 
rescaled_avg_conf = 255*(rescaled_avg_conf./highconf); 
rescaled_avg_conf(rescaled_avg_conf>255) = 255;

turbfull=flipud(turbo(512));

imwrite(rescaled_avg_conf, turbfull(19:274, :), fullfile(result_path,[num2str(length(fNames)) 'subjects_combined_conf.tif']))

figure(2); clf; imagesc((avg_confidence.*threshold_mask)); axis image; colormap(turbfull(19:274, :)); colorbar; title('Combined Confidence');
caxis([lowconffar lowconffar+highconf]); colorbar;
saveas(gcf,fullfile(result_path,[num2str(length(fNames)) 'subjects_combined_conf.svg']));
clear rescaled_avg_conf;
%% Output density image
% avg_density = avg_density; %.*sqrt(3)/2; %<-- only needed if row, not cell, spacing.
% avg_density = (1000^2).*sqrt(3)./ (2*(avg_density).^2); %<-- for cells/mm
% avg_density = sqrt(3)./ (2*(avg_density).^2); %<-- for cells/deg
avg_density(avg_density==0) = NaN;
lowdens = quantile(avg_density(threshold_mask(:)),0.01);
rescaled_avg_density = (avg_density.*threshold_mask)-lowdens;
rescaled_avg_density(rescaled_avg_density<0) = 0;
highdens = quantile(rescaled_avg_density(rescaled_avg_density~=0),0.99);
rescaled_avg_density = 255*(rescaled_avg_density./highdens); 
rescaled_avg_density(rescaled_avg_density>255) = 255;

figure(3); imagesc(avg_density.*threshold_mask); title('Combined Density'); axis image;
caxis([lowdens lowdens+highdens]); colorbar;
imwrite(rescaled_avg_density, parula(256), fullfile(result_path,[num2str(length(fNames)) 'subjects_combined_density.tif']))
saveas(gcf,fullfile(result_path,[num2str(length(fNames)) 'subjects_combined_density.svg']));
clear rescaled_avg_density;

%% Output sum map

% imwrite(uint8(255*(combined_sum_map./max(combined_sum_map(:)))), parula(256), fullfile(result_path,[num2str(length(fNames)) 'subjects_sum_map.tif']));
rescaled_sum_map = uint8( (255*(combined_sum_map./max(combined_sum_map(:)))).*threshold_mask );

figure(4); imagesc(combined_sum_map); title('Sum Map'); axis image; colormap winter; colorbar;
caxis([0 50]);
saveas(gcf,fullfile(result_path,[num2str(length(fNames)) 'subjects_sum_map.svg']));

imwrite(rescaled_sum_map, winter(256), fullfile(result_path,[num2str(length(fNames)) 'subjects_sum_map.tif']));

%% Determine average/stddev of all data.
value_map_variance = nan(global_dimension);


% For finding outliers:
figure(1); clf; hold on;
    

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

    farrangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 2.5/scaling;   
    
    confcombfar = blendederr_comb.*~farrangemask;
    confcombfar(confcombfar==0)=NaN;
    lowconffar = quantile(confcombfar(~isnan(confcombfar)), [0.10])

    confcombclose = blendederr_comb.*farrangemask;
    confcombclose(confcombclose==0)=NaN;
    lowconfclose = quantile(confcombclose(~isnan(confcombclose)), [0.01])

    closerangemask = sqrt((X-foveal_coords(f,1)).^2 + (Y-foveal_coords(f,2)).^2) <= 1.25/scaling;
    denscombclose = density_map_comb.*closerangemask;
    denscombclose(denscombclose==0)=NaN;
    densclose = quantile(denscombclose(~isnan(denscombclose)), [0.5])

    % Remove the offending values.
    density_map_comb(( (density_map_comb > densclose) & ~closerangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconffar) & ~farrangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconfclose) & farrangemask )) = NaN;
    %Also remove them from our confidence.
    blendederr_comb(isnan(density_map_comb)) = NaN;
    
    density_map_comb = density_map_comb.*threshold_mask( rowrange, colrange);
    density_map_comb(density_map_comb==0) = NaN;

    thisdiff = sum( cat(3, density_map_comb, -avg_density( rowrange, colrange)) ,3).^2; % This line is good. Nans need to be carried through.
    
%     thisdiff =  sum( cat(3, (scaling.*blendedim).*(2/sqrt(3)), -avg_spacing( rowrange, colrange)) ,3).^2;

    value_map_variance( rowrange, colrange) = sum( cat(3,value_map_variance( rowrange, colrange), thisdiff), 3,'omitnan'); 
     clear thisdiff;


    clear blendedim rowrange colrange
end
save( fullfile(result_path,[num2str(length(fNames)) '_aggregate_data.mat']), ...
                            'avg_density', 'avg_confidence','global_fovea_coords',...
                            'global_eye', 'threshold_mask', 'value_map_variance',...
                            'scaling', 'montage_rect','combined_sum_map','-v7.3');

return;
% plot(micron_position,temp_strip,'k', 'LineWidth',3);
%%
 value_std_dev = value_map_variance;
 value_std_dev = value_std_dev./(combined_sum_map-1);
 value_std_dev(isinf(value_std_dev)) =0;
 value_std_dev = sqrt(value_std_dev);
 maskedspac = value_std_dev.*threshold_mask;

figure(3); imagesc(value_std_dev.*threshold_mask); title('Combined Spacing Std dev');
saveas(gcf,[num2str(length(fNames)) 'subjects_combined_stddev.png']);
saveas(gcf,[num2str(length(fNames)) 'subjects_combined_stddev.svg']);

% maxstriplen = min([size(-global_fovea_coords global_fovea_coords]); %Find out the largest strip we can make in pixels.

deg_position = (0:maxstriplen)*scaling;

figure(20); clf; hold on;
% Superior
plot(deg_position,mean(flipud(maskedspac(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-180:global_fovea_coords(1)+180)),2,  'omitnan')); hold on;
% Inferior
plot(deg_position,mean(maskedspac(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-180:global_fovea_coords(1)+180),2, 'omitnan'));
% Nasal
nasal_strip = mean(maskedspac(global_fovea_coords(2)-180:global_fovea_coords(2)+180,global_fovea_coords(1):global_fovea_coords(1)+strip_length), 1,  'omitnan');
plot(deg_position, nasal_strip);
% Temporal
temp_strip = mean(fliplr(maskedspac(global_fovea_coords(2)-180:global_fovea_coords(2)+180,global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1,  'omitnan');
plot(deg_position,temp_strip);
legend('Superior','Inferior','Nasal','Temporal')
xlabel('Radial distance (degrees)')
ylabel('Std Dev of Density (cells/degrees^2)')

saveas(gcf, fullfile(result_path,[num2str(length(fNames)) 'subjects_directionalstddev.png']));
saveas(gcf, fullfile(result_path,[num2str(length(fNames)) 'subjects_directionalstddev.svg']));
clear maskedspac temp_strip nasal_strip

save( fullfile(result_path,[num2str(length(fNames)) '_aggregate_data.mat']), ...
                            'avg_density', 'avg_confidence','global_fovea_coords',...
                            'global_eye', 'threshold_mask', 'value_std_dev','value_map_variance',...
                            'scaling', 'montage_rect','combined_sum_map','-v7.3');

return;


%% Determine weighted average/stddev of all data.
weighted_avg_spacing = weighted_avg_spacing./weighted_avg_error;
weighted_avg_spacing_std_dev = zeros(global_dimension);
weighted_avg_error_std_dev = zeros(global_dimension);


% Find the weighted std deviations of all of the subjects.
for f=1:length(fNames)
    load( fullfile(thispath, fNames{f}), 'blendedim', 'blendederrim');
    
    rowrange = round((montage_rect{f}(1,2):montage_rect{f}(3,2))-global_fovea_coords(2)+1);
    colrange = round((montage_rect{f}(1,1):montage_rect{f}(3,1))-global_fovea_coords(1)+1);
    
    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        blendedim = fliplr(blendedim);
        blendederrim = fliplr(blendederrim);
    end
    
    % Weighted std deviation equation found here
    % https://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
    weighted_avg_spacing_std_dev( rowrange, colrange) = blendederrim.*(sum(cat(3, blendedim, -weighted_avg_spacing( rowrange, colrange)),3,'omitnan').^2);

end

weighted_std_dev_err_sum = ((combined_sum_map-1)./combined_sum_map).*weighted_avg_error;
weighted_avg_spacing_std_dev  = sqrt(weighted_avg_spacing_std_dev./weighted_std_dev_err_sum);
weighted_avg_spacing_std_dev(isinf(weighted_avg_spacing_std_dev)) =0 ;


lower01 = quantile(weighted_avg_spacing(~isnan(weighted_avg_spacing)),0.01);
upper99 = quantile(weighted_avg_spacing(~isnan(weighted_avg_spacing)),0.99);

figure(10); imagesc(weighted_avg_spacing.*threshold_mask); title('Combined Weighted Spacing');
caxis([lower01 upper99]);
figure(11); imagesc(weighted_avg_spacing_std_dev.*threshold_mask); title('Combined Weighted Spacing Weighted Std dev');


