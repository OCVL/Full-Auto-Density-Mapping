clear
close all;

% Find the global dataset
[agg_fName, globalpath ]=uigetfile(fullfile(pwd,'*.mat'),'Select the aggregate file.');

% Individuals will be above this.
individual_path = getparent(globalpath,2,'full');
[fNames] = read_folder_contents(individual_path, 'mat');

% Load the global dataset.
load(fullfile(globalpath, agg_fName));

% Find the micron/pixel equivalent to deg/pix
[~,~,lut] = determine_scaling(individual_path,fNames,[],'degrees');
persub_micronsperdegree = (291*lut{2})/24;



%% Downsample everything first- we have a lot of redundant data.
downsample_factor = 0.25;

original_size = size(avg_density);
downsampled_size = round(original_size*downsample_factor)+1;
downsamp_scaling = scaling/downsample_factor;
global_fovea_coords = round(global_fovea_coords*downsample_factor);
montage_rect = cellfun(@(c) round(c*downsample_factor), montage_rect, 'UniformOutput', false);


avg_density = imresize(avg_density, downsampled_size, 'lanczos3');
avg_confidence = imresize(avg_confidence, downsampled_size, 'lanczos3');
threshold_mask = imresize(threshold_mask, downsampled_size, 'nearest');

value_map_variance = value_map_variance./(combined_sum_map-1);
value_map_variance(isinf(value_map_variance)) =0;

value_map_variance(value_map_variance==0)=NaN;
value_map_stddev = imresize(sqrt(value_map_variance), downsampled_size, 'lanczos3');
clear value_map_variance

combined_sum_map = imresize(combined_sum_map, downsampled_size, 'nearest');

value_map_stddev = value_map_stddev.*threshold_mask;
maskeddensity = avg_density.*threshold_mask;
maskedconf = avg_confidence.*threshold_mask;

value_map_stddev(value_map_stddev==0) = NaN;
maskedconf(maskedconf==0) = NaN;
maskeddensity(maskeddensity==0) = NaN;


%% Obtain global strip averages of density.
maxstriplen = min([downsampled_size-global_fovea_coords global_fovea_coords]); %Find out the largest strip we can make in pixels.

deg_position = (0:maxstriplen)*downsamp_scaling;
strip_radius = round(1/(2*downsamp_scaling)); % Corresponds to 1 degree.
strip_length = length(deg_position)-1;

figure(10); clf; 
% Superior
sup_strip = flipud(maskeddensity(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius));
avg_sup_strip = mean(sup_strip,2, 'omitnan');
plot(deg_position,avg_sup_strip); hold on;
% Inferior
inf_strip = maskeddensity(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius);
avg_inf_strip = mean(inf_strip,2, 'omitnan');
plot(deg_position,avg_inf_strip);
% Nasal
nasal_strip = maskeddensity(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,global_fovea_coords(1):global_fovea_coords(1)+strip_length);
avg_nasal_strip = mean(nasal_strip, 1, 'omitnan');
plot(deg_position, avg_nasal_strip);
% Temporal
temporal_strip = fliplr(maskeddensity(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius, global_fovea_coords(1)-strip_length:global_fovea_coords(1)));
avg_temp_strip = mean(temporal_strip, 1, 'omitnan');
plot(deg_position,avg_temp_strip);
legend('Superior','Inferior','Nasal','Temporal')
xlabel('Radial distance (degrees)')
ylabel('Density (cells/degrees^2)')
hold off;

saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgdens.svg']) );
saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgdens.png']) );


% avg_tempnasal_strip = mean([nasal_strip;temporal_strip], 1, 'omitnan');
% plot(deg_position,avg_tempnasal_strip);
% avg_supinf_strip = mean([sup_strip inf_strip], 2, 'omitnan');
% plot(deg_position,avg_supinf_strip);
% legend('Temporal/Nasal','Superior/Inferior')
% xlabel('Radial distance (degrees)')
% ylabel('Density (cells/degrees^2)')


%% STDDEV
figure(100); clf; hold on;

sup_stddev = mean(flipud(value_map_stddev(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius)), 2, 'omitnan');
plot(deg_position, sup_stddev);
inf_stddev = mean(value_map_stddev(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius), 2, 'omitnan');
plot(deg_position, inf_stddev); 
nasal_stddev = mean(value_map_stddev(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,global_fovea_coords(1):global_fovea_coords(1)+strip_length), 1, 'omitnan');
plot(deg_position,nasal_stddev); 
temp_stddev = mean(fliplr(value_map_stddev(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius, global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1, 'omitnan');
plot(deg_position,temp_stddev); 

legend('Superior','Inferior','Nasal','Temporal')
xlabel('Radial distance (degrees)')
ylabel('Density (cells/degrees^2)')

% upperstd_temp_strip = avg_temp_strip+temp_nasal_stddev;
% plot(deg_position,upperstd_temp_strip);
% lowerstd_temp_strip = avg_temp_strip-temp_nasal_stddev;
% plot(deg_position,lowerstd_temp_strip);
hold off;


%% Obtain global strip averages of confidence.

figure(11); clf; 
% Superior
avg_sup_strip_conf = mean(flipud(maskedconf(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius)),2, 'omitnan');
plot(deg_position,avg_sup_strip_conf); hold on;
% Inferior
avg_inf_strip_conf = mean(maskedconf(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius),2, 'omitnan');
plot(deg_position,avg_inf_strip_conf);
% Nasal
avg_nasal_strip_conf = mean(maskedconf(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,global_fovea_coords(1):global_fovea_coords(1)+strip_length), 1, 'omitnan');
plot(deg_position, avg_nasal_strip_conf);
% Temporal
avg_temp_strip_conf = mean(fliplr(maskedconf(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius, global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1, 'omitnan');
plot(deg_position,avg_temp_strip_conf);
legend('Superior','Inferior','Nasal','Temporal')
xlabel('Radial distance (degrees)')
ylabel('Average confidence (AU)')
axis([0 10 0 1])
hold off;

saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgconf.svg']) );
saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgconf.png']) );

%% Perform a pairwise analysis of the difference of each person from the average strips
figure(12); clf; hold on;
xlabel('Radial distance (degrees)')
ylabel('Density (cells/degrees^2)')
% figure(13); clf; hold on;
% xlabel('Radial distance (degrees)')
% ylabel('Confidence (AU)')

fList = 1:length(fNames);
% fList = [36,40,45,49];
numcones_micrad = cell(size(fList));
numcones_degrad = cell(size(fList));
num_cones_in_annulus = cell(size(fList));


for f=fList
    disp(['Loading, downsampling and prepping: ' fNames{f}])
    load( fullfile(individual_path, fNames{f}), 'density_map_comb', 'blendederr_comb','fovea_coords', 'imsize');
    
    indiv_size = [montage_rect{f}(3,2)-montage_rect{f}(1,2)+1 montage_rect{f}(3,1)-montage_rect{f}(1,1)+1];
    indiv_shifted_density = nan(downsampled_size);
    indiv_shifted_confidence = nan(downsampled_size);

    mm_position = (0:maxstriplen)*(downsamp_scaling*persub_micronsperdegree(f));
    
    rowrange = ceil((montage_rect{f}(1,2):montage_rect{f}(3,2))+global_fovea_coords(2));
    colrange = ceil((montage_rect{f}(1,1):montage_rect{f}(3,1))+global_fovea_coords(1));

    if isempty(strfind(fNames{f}, global_eye)) % If it doesn't match our eye, flip the montage data.
        density_map_comb = fliplr(density_map_comb);
        blendederr_comb = fliplr(blendederr_comb);
        fovea_coords(1) = imsize(2)-fovea_coords(1);
    end    

    % Lightly filter our data to remove any outliers in confidence or
    % density.
    [X, Y] = meshgrid(1:size(density_map_comb,2), 1:size(density_map_comb,1));

    farrangemask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 2.5/scaling;   
    
    confcombfar = blendederr_comb.*~farrangemask;
    confcombfar(confcombfar==0)=NaN;
    lowconffar = quantile(confcombfar(~isnan(confcombfar)), [0.10]);

    confcombclose = blendederr_comb.*farrangemask;
    confcombclose(confcombclose==0)=NaN;
    lowconfclose = quantile(confcombclose(~isnan(confcombclose)), [0.01]);

    closerangemask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 1.25/scaling;
    denscombclose = density_map_comb.*closerangemask;
    denscombclose(denscombclose==0)=NaN;
    densclose = quantile(denscombclose(~isnan(denscombclose)), [0.5]);

    % If we're looking at data with these filenames, drop their long
    % distance data as its quality is questionable.
    if contains(fNames{f},'11101') || contains(fNames{f},'11092')
        good_data_mask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2)  <= 8/scaling;
        density_map_comb = density_map_comb .*good_data_mask;
        blendederr_comb = blendederr_comb .*good_data_mask;
    end
    clear X Y
    
    % Remove the offending values.
    density_map_comb(( (density_map_comb > densclose) & ~closerangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconffar) & ~farrangemask )) = NaN;
    density_map_comb(( (blendederr_comb < lowconfclose) & farrangemask )) = NaN;
    %Also remove them from our confidence.
    blendederr_comb(isnan(density_map_comb)) = NaN;


    if min(rowrange) == 0
        rowrange=rowrange+1;
    end
    if min(colrange) == 0
        colrange=colrange+1;
    end
    if max(rowrange) > downsampled_size(1)
        rowrange=downsampled_size(1);
    end
    if max(colrange) > downsampled_size(2)
        colrange=downsampled_size(2);
    end

    % Downsample the dataset at the same rate the average was, and mask it.
    density_map_comb = imresize(density_map_comb, indiv_size, 'lanczos3');
    blendederr_comb = imresize(blendederr_comb, indiv_size, 'lanczos3');
    figure(121); imagesc(density_map_comb); axis image; hold on;
    plot(fovea_coords(1)*.25, fovea_coords(2)*.25, 'r*')
    hold off;
    
    % Put the dataset it the position we had it before.
    indiv_shifted_density( rowrange, colrange) =  density_map_comb;
    indiv_shifted_confidence( rowrange, colrange) = blendederr_comb;

%     indiv_shifted_density = indiv_shifted_density.*threshold_mask;
%     indiv_shifted_confidence = indiv_shifted_confidence.*threshold_mask;
    
    indiv_shifted_density(indiv_shifted_density==0) = NaN;
    indiv_shifted_confidence(indiv_shifted_confidence==0) = NaN;

    
    
    % Create and extract total cone annuli.    
    [colmesh, rowmesh] = meshgrid(1:size(indiv_shifted_density,2), 1:size(indiv_shifted_density,1));
    
    distmesh = sqrt((colmesh-global_fovea_coords(1)).^2 + ...
                    (rowmesh-global_fovea_coords(2)).^2);

    maxdist = max(distmesh(:));
    radius = 5; %In pixels
    
    r=1;    
    radii = 0:radius:maxdist;
    
    micrad = nan(size(radii));
    degrad = nan(size(radii));
    num_cia = nan(size(radii));
    
    parfor r=1:length(radii)-1
    
%     while inner_rad + radius < maxdist
    
        aoi = double(distmesh > radii(r) & distmesh <= radii(r+1));
        aoi(aoi == 0) = NaN;
        
        
    
        annular_area = sum(aoi(:).*(downsamp_scaling*downsamp_scaling), 'omitnan');
        
        aoi_density = indiv_shifted_density.*aoi;
%         imagesc(aoi_density)
        
        degrad(r) = downsamp_scaling*(radii(r+1)+radii(r))/2;
        micrad(r) = persub_micronsperdegree(f)*degrad(r);
        num_cia(r) = mean(aoi_density(:), 'omitnan').*annular_area;
        
%         r=r+1;
%         inner_rad = inner_rad + radius;
    end
    
    numcones_degrad{f} = degrad;
    numcones_micrad{f} = micrad;
    num_cones_in_annulus{f} = num_cia;
    
    tot_cones = cumsum(num_cones_in_annulus{f},'omitnan');
    tot_cones(isnan(num_cia)) = NaN;
    
    num_cones_in_annulus{f} = tot_cones;
    
    figure(12); 
    subplot(2,1,1);hold on;
    plot(numcones_degrad{f}, tot_cones);
    subplot(2,1,2);hold on;
    plot(numcones_micrad{f}, tot_cones);
    drawnow;

    
%     pause;
end
figure(12); hold off;


totcones = cell2mat(num_cones_in_annulus')';
degsubs = cell2mat(numcones_degrad')';
mmsubs = cell2mat(numcones_micrad')';

combined = nan(size(totcones,1), size(totcones,2)*3);
combined(:,1:3:end) = degsubs;
combined(:,2:3:end) = mmsubs;
combined(:,3:3:end) = totcones;

dattable = table();
indarr = 1:3:size(combined,2)+1;
for j=1:length(indarr)-1
    beginv = indarr(j);
    endv = indarr(j+1)-1;
    dattable.(fNames{j}) = combined(:,beginv:endv);

end
writetable(dattable, 'totalcones_table.csv');
% figure(13); hold off;