%% Robert F Cooper
% 1/13/21
%
% This script's sole job is to run multiple montages in series; to save time, 
% instead of operating on individual files though, it operates on a series
% of folders.


thisfolder = pwd;

[foveafile, foveafold] = uigetfile(fullfile(pwd,'*.xlsx'), 'Select the folders containing the foveas you wish to assess.');

[~,filelist,~] = xlsread(fullfile(foveafold,foveafile));

[lutfname, lutfolder] = uigetfile(fullfile(pwd,'*.csv'),'Select scaling LUT, OR cancel if you want to input the scale directly.');

unit='microns (mm density)';

%%
restartf=1;
endf=length(filelist);

%%
for f=restartf:endf
    restartf=f;

    something=false;
    if ~isempty(filelist{f,1})
        
        splitlist=split(filelist{f,1},'_');
        filedir=fullfile(lutfolder, [splitlist{2} '_' splitlist{4}],'confocal');
        disp(['***** Running analysis on: ' [splitlist{2} '_' splitlist{4}] ' *****']);

        [ scaling, ~, lut ]=determine_scaling(filedir, filelist(f,1), fullfile(lutfolder,lutfname) ,unit );

        foveaim = imread(fullfile(filedir,filelist{f,1}));
        
        % Look 
        all_scales = nan(length(filelist),1);

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
        
        imsize = round(size(foveaim).*rel_scale);

        if contains(filelist{1},'11102_OD')  % This factor was added to account for trial lenses.
            imsize = round( imsize(1:2) * 1.13 );
        elseif contains(filelist{1},'11051_OD')  % This factor was added to account for trial lenses.
            imsize = round( imsize(1:2) * 1.1369 );
        end
        
        foveaim = double(imresize(foveaim, imsize(1:2)));
        
        [spac, spacmap, confmap, summap, imbox] = fit_fourier_spacing(foveaim, 96, false);
        
        imsize = size(foveaim);
        densim = zeros(imsize(1:2));
        
        weightedmap = spacmap./confmap;
        
        density_map = sqrt(3)./ (2*(weightedmap*scaling).^2);
        density_map(isinf(density_map))=0;
        
        densim(imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3)) = density_map;
        
        result_path = fullfile(filedir,'Results_fovea_map');
        mkdir(result_path)

       
        lower01 = quantile(density_map(~isnan(density_map)),0.01);
        
        scaled_density = density_map-lower01;
        upper99 = quantile(scaled_density(~isnan(scaled_density)),0.99);
        
        scaled_density = 255.*(scaled_density./ upper99 );
        scaled_density(scaled_density>255) = 255;
        
        imwrite(uint8(scaled_density),parula(256), fullfile(result_path, [splitlist{2} '_' splitlist{4} '_Fovea.tif']));
        safesave(fullfile(result_path, [splitlist{2} '_' splitlist{4} '_Fovea.mat']),...
            weightedmap,density_map,densim, spacmap, confmap, summap, imbox);
        
    end
end

function []= safesave(newsave, x1, x2, x3, x4, x5, x6, x7)
    disp(['Saving to:' newsave]);
    weightedmap = x1;
    density_map = x2;
    densim = x3;
    spacmap = x4;
    confmap = x5;
    summap = x6;
    imbox = x7;

    save(newsave, 'weightedmap','density_map','densim', 'spacmap', 'confmap', 'summap', 'imbox', '-v7.3')
end