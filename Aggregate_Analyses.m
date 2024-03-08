clear
close all;

% Add our support libraries to the path.
basePath = which('Aggregate_Analyses.m');
[basePath ] = fileparts(basePath);
path(path, fullfile(basePath,'lib'));

% Find the global dataset.
[agg_fName, globalpath ]=uigetfile(fullfile(pwd,'*.mat'),'Select the aggregate file.');

% Individuals should be above this.
individual_path = getparent(globalpath,2,'full');
[fNames] = read_folder_contents(individual_path, 'mat');

% Load the global dataset.
load(fullfile(globalpath, agg_fName));

% Find the micron/pixel equivalent to deg/pix
[scalinginfo,~,lut] = determine_scaling(individual_path,fNames,[],'microns (mm density)');
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
value_map_variance = imresize((value_map_variance), downsampled_size, 'lanczos3');
value_map_stddev = imresize(sqrt(value_map_variance), downsampled_size, 'lanczos3');

combined_sum_map = imresize(combined_sum_map, downsampled_size, 'nearest');

value_map_variance = value_map_variance.*threshold_mask;
value_map_stddev = value_map_stddev.*threshold_mask;
maskeddensity = avg_density.*threshold_mask;
maskedconf = avg_confidence.*threshold_mask;

value_map_variance(value_map_variance==0)= NaN;
value_map_stddev(value_map_stddev==0) = NaN;
maskedconf(maskedconf==0) = NaN;
maskeddensity(maskeddensity==0) = NaN;

maxstriplen = min([downsampled_size-global_fovea_coords global_fovea_coords]); %Find out the largest strip we can make in pixels.

%% Obtain global strip averages of density.

position = (0:maxstriplen)*downsamp_scaling;
strip_radius = round(1/(2*downsamp_scaling)); 
strip_length = length(position)-1;

figure(10); clf; 
% Superior
sup_strip = flipud(maskeddensity(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius));
avg_sup_strip = mean(sup_strip,2, 'omitnan');
plot(position,avg_sup_strip); hold on;
% Inferior
inf_strip = maskeddensity(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius);
avg_inf_strip = mean(inf_strip,2, 'omitnan');
plot(position,avg_inf_strip);
% Nasal
nasal_strip = maskeddensity(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,global_fovea_coords(1):global_fovea_coords(1)+strip_length);
avg_nasal_strip = mean(nasal_strip, 1, 'omitnan');
plot(position, avg_nasal_strip);
% Temporal
temporal_strip = fliplr(maskeddensity(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius, global_fovea_coords(1)-strip_length:global_fovea_coords(1)));
avg_temp_strip = mean(temporal_strip, 1, 'omitnan');
plot(position,avg_temp_strip);
legend('Superior','Inferior','Nasal','Temporal')


if strcmp(unit, 'microns (mm density)')    
    xlabel('Radial distance (microns)')
    ylabel('Density (cells/mm^2)')
elseif strcmp(unit, 'degrees')
    xlabel('Radial distance (degrees)')
    ylabel('Density (cells/degrees^2)')
end
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

%% Create and extract total cone annuli from average map.   
[colmesh, rowmesh] = meshgrid(1:size(maskeddensity,2), 1:size(maskeddensity,1));

horzmask = abs(colmesh-global_fovea_coords(1))>=abs(rowmesh-global_fovea_coords(2));
vertmask = abs(colmesh-global_fovea_coords(1))<abs(rowmesh-global_fovea_coords(2));

distmesh = sqrt((colmesh-global_fovea_coords(1)).^2 + ...
                (rowmesh-global_fovea_coords(2)).^2);

maxdist = max(distmesh(:));
radius = 5; %In pixels

r=1;    
radii = 0:radius:maxdist;

micrad = nan(size(radii));
rad = nan(size(radii));
num_cia = nan(size(radii));
num_cia_horz = nan(size(radii));
num_cia_vert = nan(size(radii));
num_cia_sep = nan(size(radii));

parfor r=1:length(radii)-1

%     while inner_rad + radius < maxdist

    aoi = double(distmesh > radii(r) & distmesh <= radii(r+1));
    aoi(aoi == 0) = NaN;
    
    horzaoi = aoi.*horzmask;
    horzaoi(horzaoi == 0) = NaN;
    vertaoi = aoi.*vertmask;
    vertaoi(vertaoi == 0) = NaN;

    annular_area = sum(aoi(:).*(downsamp_scaling*downsamp_scaling), 'omitnan');
    horz_annular_area = sum(horzaoi(:).*(downsamp_scaling*downsamp_scaling), 'omitnan');
    vert_annular_area = sum(vertaoi(:).*(downsamp_scaling*downsamp_scaling), 'omitnan');

    aoi_density = maskeddensity.*aoi;
    horz_aoi_density = maskeddensity.*horzaoi;
    vert_aoi_density = maskeddensity.*vertaoi;    
    
    rad(r) = downsamp_scaling*(radii(r+1)+radii(r))/2;
    num_cia(r) = mean(aoi_density(:), 'omitnan').*annular_area;

    num_cia_horz(r) = mean(horz_aoi_density(:), 'omitnan').*horz_annular_area;
    num_cia_vert(r) = mean(vert_aoi_density(:), 'omitnan').*vert_annular_area;
    num_cia_sep(r) = num_cia_horz(r) + num_cia_vert(r);
    
%         r=r+1;
%         inner_rad = inner_rad + radius;
end


tot_cones = cumsum(num_cia,'omitnan');
tot_cones(isnan(num_cia)) = NaN;

tot_cones_horz= cumsum(num_cia_horz,'omitnan');
tot_cones_horz(isnan(num_cia_horz)) = NaN;
tot_cones_vert= cumsum(num_cia_vert,'omitnan');
tot_cones_vert(isnan(num_cia_vert)) = NaN;

tot_cones_sep = cumsum(num_cia_sep,'omitnan');
tot_cones_sep(isnan(num_cia_sep)) = NaN;


figure(12); 
subplot(2,1,1);hold on;
plot(rad, tot_cones);
subplot(2,1,2);hold on;
plot(rad, tot_cones_sep,rad,tot_cones_horz, rad, tot_cones_vert);
drawnow; hold off;


%% STDDEV
figure(100); clf; hold on;

sup_stddev = sqrt(mean(flipud(value_map_variance(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius)), 2, 'omitnan'));
plot(position, sup_stddev);
inf_stddev = sqrt(mean(value_map_variance(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius), 2, 'omitnan'));
plot(position, inf_stddev); 
nasal_stddev = sqrt(mean(value_map_variance(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,global_fovea_coords(1):global_fovea_coords(1)+strip_length), 1, 'omitnan'));
plot(position,nasal_stddev); 
temp_stddev = sqrt(mean(fliplr(value_map_variance(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius, global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1, 'omitnan'));
plot(position,temp_stddev); 

legend('Superior','Inferior','Nasal','Temporal')
if strcmp(unit, 'microns (mm density)')    
    xlabel('Radial distance (microns)')
    ylabel('Density (cells/mm^2)')
elseif strcmp(unit, 'degrees')
    xlabel('Radial distance (degrees)')
    ylabel('Density (cells/degrees^2)')
end

hold off;

saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_stddev.svg']) );
saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_stddev.png']) );


%% Obtain global strip averages of confidence.

figure(11); clf; 
% Superior
avg_sup_strip_conf = mean(flipud(maskedconf(global_fovea_coords(2)-strip_length:global_fovea_coords(2),global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius)),2, 'omitnan');
plot(position,avg_sup_strip_conf); hold on;
% Inferior
avg_inf_strip_conf = mean(maskedconf(global_fovea_coords(2):global_fovea_coords(2)+strip_length, global_fovea_coords(1)-strip_radius:global_fovea_coords(1)+strip_radius),2, 'omitnan');
plot(position,avg_inf_strip_conf);
% Nasal
avg_nasal_strip_conf = mean(maskedconf(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius,global_fovea_coords(1):global_fovea_coords(1)+strip_length), 1, 'omitnan');
plot(position, avg_nasal_strip_conf);
% Temporal
avg_temp_strip_conf = mean(fliplr(maskedconf(global_fovea_coords(2)-strip_radius:global_fovea_coords(2)+strip_radius, global_fovea_coords(1)-strip_length:global_fovea_coords(1))), 1, 'omitnan');
plot(position,avg_temp_strip_conf);
legend('Superior','Inferior','Nasal','Temporal')

if strcmp(unit, 'microns (mm density)')    
    xlabel('Radial distance (microns)')
elseif strcmp(unit, 'degrees')
    xlabel('Radial distance (degrees)')
end
ylabel('Average confidence (AU)')

hold off;

saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgconf.svg']) );
saveas(gcf, fullfile(globalpath,[num2str(length(fNames)) 'subjects_directional_avgconf.png']) );



%% Analyze the total number of cones as a function of eccentricity by fitting our curves

mindensity = nan(size(fNames));
maxdensity = nan(size(fNames));
numcones_micrad = cell(size(fNames));
numcones_degrad = cell(size(fNames));
num_cones_in_annulus = cell(size(fNames));
num_cones_in_annulus_sep = cell(size(fNames));

figure(12); clf; hold on;
figure(9); clf;
if strcmp(unit, 'microns (mm density)')    
    xlabel('Radial distance (mm)')
    ylabel('Density (cells/mm^2)')
elseif strcmp(unit, 'degrees')
    xlabel('Radial distance (degrees)')
    ylabel('Density (cells/degrees^2)')
end

%%
for f=1:length(fNames)
    disp(['Loading, downsampling and prepping: ' fNames{f}])
    load( fullfile(individual_path, fNames{f}), 'density_map_comb', 'blendederr_comb','fovea_coords', 'imsize');
    
    indiv_size = [montage_rect{f}(3,2)-montage_rect{f}(1,2)+1 montage_rect{f}(3,1)-montage_rect{f}(1,1)+1];
    indiv_shifted_density = nan(downsampled_size);
    indiv_shifted_confidence = nan(downsampled_size);
   

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

    if strcmp(unit, 'microns (mm density)')
        farrangemask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 800/scaling; 
    elseif strcmp(unit, 'degrees')
        farrangemask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 2.5/scaling;
    end
    
    confcombfar = blendederr_comb.*~farrangemask;
    confcombfar(confcombfar==0)=NaN;
    lowconffar = quantile(confcombfar(~isnan(confcombfar)), [0.1]);

    confcombclose = blendederr_comb.*farrangemask;
    confcombclose(confcombclose==0)=NaN;
    lowconfclose = quantile(confcombclose(~isnan(confcombclose)), [0.01]);

    if strcmp(unit, 'microns (mm density)')
        closerangemask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 400/scaling;
    elseif strcmp(unit, 'degrees')
        closerangemask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 1.25/scaling;
    end
    denscombclose = density_map_comb.*closerangemask;
    denscombclose(denscombclose==0)=NaN;
    densclose = quantile(denscombclose(~isnan(denscombclose)), [0.5]);

    % If we're looking at data with these filenames, drop their long
    % distance data as its quality is questionable.
    if contains(fNames{f},'11101') || contains(fNames{f},'11092')
        if strcmp(unit, 'microns (mm density)')
            good_data_mask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 2100/scaling; 
        elseif strcmp(unit, 'degrees')
            good_data_mask = sqrt((X-fovea_coords(1)).^2 + (Y-fovea_coords(2)).^2) <= 7/scaling;
        end

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
    plot(fovea_coords(1)*downsample_factor, fovea_coords(2)*downsample_factor, 'r*')
    hold off;
    
    % Put the dataset it the position we had it before.
    indiv_shifted_density( rowrange, colrange) =  density_map_comb;
    indiv_shifted_confidence( rowrange, colrange) = blendederr_comb;

    indiv_shifted_density(indiv_shifted_density==0) = NaN;
    indiv_shifted_confidence(indiv_shifted_confidence==0) = NaN;

    [colmesh, rowmesh] = meshgrid(1:size(indiv_shifted_density,2), 1:size(indiv_shifted_density,1));
    colmesh = colmesh-global_fovea_coords(1);
    rowmesh = rowmesh-global_fovea_coords(2);
    horzmask = abs(colmesh)>=abs(rowmesh);
    vertmask = abs(colmesh)<abs(rowmesh);

    distmesh = sqrt(colmesh.^2 + rowmesh.^2);
    maxdist = max(distmesh(:));

    if strcmp(unit, 'microns (mm density)')
        colmesh = (colmesh.*downsamp_scaling)/1000;
        rowmesh = (rowmesh.*downsamp_scaling)/1000;
    elseif strcmp(unit, 'degrees')
        colmesh = colmesh.*downsamp_scaling;
        rowmesh = rowmesh.*downsamp_scaling;
    end

    radius = 5; %In pixels

    % Extract all horizontal-dominant and vertical-dominant data for
    % fitting.
    horzpolar = imcart2pseudopolar(indiv_shifted_density.*horzmask, 1, 0.5,global_fovea_coords,'cubic',0,maxstriplen);
    horzpolar(horzpolar==0) = NaN;
    horzpolar = movmean(mean(horzpolar,'omitnan'), radius);
    vertpolar = imcart2pseudopolar(indiv_shifted_density.*vertmask, 1, 0.5,global_fovea_coords,'cubic',0,maxstriplen);
    vertpolar(vertpolar==0) = NaN;
    vertpolar = movmean(mean(vertpolar,'omitnan'), radius);
    
    if strcmp(unit, 'microns (mm density)')
        position = ((0:maxstriplen)*(downsamp_scaling))/1000;
    
        plot(position,horzpolar,position,vertpolar)
        hold on;
    
        [pks,locs]=findpeaks(horzpolar,'MinPeakProminence',1000);
        locs = locs(pks>55000);
        pks = pks(pks>55000);

        if ~isempty(locs)
            horzcutoff = locs(end);
            horzpk = pks(end);
        else
            [horzpk, horzcutoff]=max(horzpolar);
        end
        
    
        [pks,locs]=findpeaks(vertpolar,'MinPeakProminence',1000);
        locs = locs(pks>55000);
        pks = pks(pks>55000);
        
        if ~isempty(locs)
            vertcutoff = locs(end);
            vertpk = pks(end);
        else
            [vertpk, vertcutoff]=max(horzpolar);
        end

        plot(position(horzcutoff),horzpk,'b*')
        plot(position(vertcutoff),vertpk,'r*')
        hold off;
        drawnow;
    elseif strcmp(unit, 'degrees')
        position = (0:maxstriplen)*(downsamp_scaling);

        plot(position,horzpolar,position,vertpolar)
        hold on;
    
        [pks,locs]=findpeaks(horzpolar,'MinPeakProminence',1000);
        locs = locs(pks>4000);
        pks = pks(pks>4000);

        if ~isempty(locs)
            horzcutoff = locs(end);
            horzpk = pks(end);
        else
            [horzpk, horzcutoff]=max(horzpolar);
        end
        
    
        [pks,locs]=findpeaks(vertpolar,'MinPeakProminence',1000);
        locs = locs(pks>4000);
        pks = pks(pks>4000);
        
        if ~isempty(locs)
            vertcutoff = locs(end);
            vertpk = pks(end);
        else
            [vertpk, vertcutoff]=max(horzpolar);
        end

        plot(position(horzcutoff),horzpk,'b*')
        plot(position(vertcutoff),vertpk,'r*')
        hold off;
        drawnow;
    end
    
    horzpolar(1:horzcutoff) = NaN;
    vertpolar(1:vertcutoff) = NaN;
    all_mm = [fliplr(-position) position];
    allpolar = [fliplr(vertpolar) horzpolar]/max([fliplr(vertpolar) horzpolar]);

    horz_position = position(horzcutoff+1:end);
    horzpolar = horzpolar(horzcutoff+1:end);
    vert_position = position(vertcutoff+1:end);
    vertpolar = vertpolar(vertcutoff+1:end);

    allmax = max(max(horzpolar), max(vertpolar));
    horzpolar=horzpolar/allmax;
    vertpolar=vertpolar/allmax;


    coeffs = fitLinkedCauchy(horz_position, horzpolar,vert_position,vertpolar);
    
    fitted_density = coeffs(1)./((1+(colmesh./coeffs(3)).^2+(rowmesh./coeffs(8)).^2)) + ...
                    coeffs(2)./((1+(colmesh./coeffs(4)).^2+(rowmesh./coeffs(9)).^2)) + (coeffs(5) + coeffs(10))/2;

    fitted_density  = allmax.*fitted_density ; % 'Un-normalize' the values

    % figure(122);
    % imagesc(fitted_density);


    % Create and extract total cone annuli.    
    r=5;    
    radii = 0:r:maxstriplen;
    
    micrad = nan(size(radii));
    rad = nan(size(radii));
    num_cia = nan(size(radii));
    num_cia_horz = nan(size(radii));
    num_cia_vert = nan(size(radii));
    num_cia_sep = nan(size(radii));
    
% figure(9000); hold on;


    parfor r=1:length(radii)-1
        aoi = double(distmesh > radii(r) & distmesh <= radii(r+1));
        aoi(aoi == 0) = NaN;
        
        horzaoi = aoi.*horzmask;
        horzaoi(horzaoi == 0) = NaN;
        vertaoi = aoi.*vertmask;
        vertaoi(vertaoi == 0) = NaN;
    
        annular_area = sum(aoi(:).*(downsamp_scaling*downsamp_scaling), 'omitnan');
        horz_annular_area = sum(horzaoi(:).*(downsamp_scaling*downsamp_scaling), 'omitnan');
        vert_annular_area = sum(vertaoi(:).*(downsamp_scaling*downsamp_scaling), 'omitnan');

        aoi_density = fitted_density.*aoi;
        horz_aoi_density = fitted_density.*horzaoi;
        vert_aoi_density = fitted_density.*vertaoi;

%         r*downsamp_scaling          
%         mean(aoi_density(:), 'omitnan').*annular_area
%         sum(aoi_density(:).*(downsamp_scaling*downsamp_scaling), 'omitnan')
        % plot(r, mean(aoi_density(:), 'omitnan'),'r*');
        % plot(r, annular_area,'k.');
        rad(r) = downsamp_scaling*(radii(r+1)+radii(r))/2;
        if strcmp(unit, 'microns (mm density)')
            num_cia(r) = mean(aoi_density(:), 'omitnan').*(annular_area/1000^2);
            num_cia_horz(r) = mean(horz_aoi_density(:), 'omitnan').*(horz_annular_area/1000^2);
            num_cia_vert(r) = mean(vert_aoi_density(:), 'omitnan').*(vert_annular_area/1000^2);
            num_cia_sep(r) = num_cia_horz(r) + num_cia_vert(r);
    
        
        elseif strcmp(unit, 'degrees')
            num_cia(r) = mean(aoi_density(:), 'omitnan').*annular_area;
            num_cia_horz(r) = mean(horz_aoi_density(:), 'omitnan').*horz_annular_area;
            num_cia_vert(r) = mean(vert_aoi_density(:), 'omitnan').*vert_annular_area;
            num_cia_sep(r) = num_cia_horz(r) + num_cia_vert(r);
        end
    end
    
    maxdensity(f) = max(fitted_density(:),[],'omitnan');
    mindensity(f) = min(fitted_density(:),[],'omitnan');

    numcones_rad{f} = rad;
    num_cones_in_annulus{f} = num_cia;
    num_cones_in_annulus_sep{f} = num_cia_sep;
    
    tot_cones = cumsum(num_cones_in_annulus{f},'omitnan');
    tot_cones(isnan(num_cia)) = NaN;
    
    num_cones_in_annulus{f} = tot_cones;

    tot_cones = cumsum(num_cones_in_annulus_sep{f},'omitnan');
    tot_cones(isnan(num_cia)) = NaN;

    num_cones_in_annulus_sep{f} = tot_cones;
    
    figure(12); hold on;
    plot(numcones_rad{f}, num_cones_in_annulus{f});
    drawnow;

    figure(9);
    subplot(1,2,1);  hold on;
    plot( mindensity(f), max(num_cones_in_annulus{f},[],'omitnan'), 'k*');
    xlabel('Min Density')
    subplot(1,2,2); hold on;
    plot( maxdensity(f), max(num_cones_in_annulus{f},[],'omitnan'), 'k*');
    xlabel('Max Density')
    drawnow;
     pause(1);
end
figure(12); hold off;


%%


totcones = cell2mat(num_cones_in_annulus)';
subs = cell2mat(numcones_rad')';

figure(99);
subplot(1,2,1); 
plot( mindensity, totcones(end-1,:),'b.');
xlabel('Min Density')
subplot(1,2,2);
plot( maxdensity, totcones(end-1,:), 'r.');
xlabel('Max Density')
drawnow;
dlmwrite(fullfile(globalpath, [char(datetime('today')) '_minmax_density_v_tot.csv']), [totcones(end-1,:)' mindensity maxdensity], ',');


combined = nan(size(totcones,1), size(totcones,2)*3);
combined(:,2:3:end) = subs;
combined(:,3:3:end) = totcones;

dattable = table();
indarr = 1:3:size(combined,2)+3;
for j=1:length(indarr)-1
    beginv = indarr(j);
    endv = indarr(j+1)-1;
    dattable.(fNames{j}) = combined(:,beginv:endv);

end

writetable(dattable, fullfile(globalpath, [char(datetime('today')) '_totalcones_table.csv']));
