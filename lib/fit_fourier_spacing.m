function [avg_pixel_spac, interped_spac_map, interped_conf_map, sum_map, imbox ] = fit_fourier_spacing(test_image, roi_size, ref_loc, rdivider)
% FUNCTION [avg_pixel_spac, interped_spac_map, interped_conf_map, sum_map, imbox ] = fit_fourier_spacing(test_image, roi_size, ref_loc, rdivider)
% 
% #### Inputs:
% 
% - **test_image**: The image that will be analyzed. The only requirement is that it is a 2d, grayscale (1 channel) image.
% - **roi_size**: The side length (in pixels) of a sliding roi window- The roi will march along the image you've provided at a rate of 1/rdivider the size of the ROI, creating a "map" of spacing of the image.
% - **ref_loc**: If supplied alongside a function definition for roi_size, then each calculated roi size will be referenced to a this coordinate.
% - **rdivider**: The roi_size divider of our sliding window relative to our roi_size. If undefined, it will default to 8, corresponding to 1/8 of the roi_size.
% 
% #### Outputs:
% 
% - **avg_pixel_spac**: The average spacing of the image.
% - **interped_spac_map**: The spacing map of the input image (in pixel spacing).
% - **interped_conf_map**: The confidence map of the input image.
% - **sum_map**: The map corresponding to the amount of ROI overlap across the output map.
% - **imbox**: The bounding region of valid (nonzero, NaN, or Inf) pixels.
%
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
%

if ~exist('test_image','var') || isempty(test_image)
    [filename, pathname] = uigetfile('*.tif', 'Pick an image to assess');

    test_image = imread( fullfile(pathname, filename) );
end

if size(test_image,3) >1
    test_image = test_image(:,:,1);
end
im_size = size(test_image);

if ~exist('roi_size','var')
    roi_size = 128;
end

if ~exist('ref_loc','var')
    ref_loc = 0;
end

if ~exist('rdivider','var')
    rdivider = 8;    
end


imcomps = bwconncomp( imclose(test_image>0,ones(5)) );
imbox = regionprops(imcomps, 'BoundingBox');

boxsizes = zeros(size(imbox,1),1);
for i=1:size(imbox,1)
    boxsizes(i)= imbox(i).BoundingBox(3)*imbox(i).BoundingBox(4);
end   
[~, maxsizeind]=max(boxsizes);
imbox = floor(imbox(maxsizeind).BoundingBox);

imbox(imbox<=0) = 1;
width_diff = im_size(2)-(imbox(1)+imbox(3));
if width_diff  < 0 
    imbox(3) = imbox(3)+width_diff;
end
height_diff = im_size(1)-(imbox(2)+imbox(4));
if height_diff  < 0 
    imbox(4) = imbox(4)+height_diff;
end


if ~isa(roi_size, 'function_handle') && any( im_size(1:2) <= roi_size)
    % Our roi size should always be divisible by 2 (for simplicity).
    if rem(max(roi_size),2) ~= 0
        roi_size = max(roi_size)-1;
    end
    
    % If our image is oblong and smaller than a default roi_size, then pad
    % it to be as large as the largest edge.
    padsize = ceil((max(roi_size)-im_size)/2);
    roi = {padarray(test_image,padsize,'both')};
    
else
       
    if ~isa(roi_size, 'function_handle')
        roi_step = floor(roi_size/rdivider);

        roi = cell(round((size(test_image)-roi_size)/roi_step));
        
        % Our roi size should always be divisible by 2 (for simplicity).
        if rem(roi_size,2) ~= 0
            roi_size = roi_size-1;
        end
        
        for i=imbox(2):roi_step:imbox(2)+imbox(4)-roi_size
            for j=imbox(1):roi_step:imbox(1)+imbox(3)-roi_size

                numzeros = sum(sum(test_image(i:i+roi_size-1, j:j+roi_size-1)<=10));

                if numzeros < (roi_size*roi_size)*0.05
                    roi{round(i/roi_step)+1,round(j/roi_step)+1} = test_image(i:i+roi_size-1, j:j+roi_size-1);
                else
                    roi{round(i/roi_step)+1,round(j/roi_step)+1} =[];
                end
            end
        end
    else %If roi_size is a function, that means we'll expect it to vary.
        
        j=round(imbox(1));
        i=round(imbox(2)); %ref_loc takes on the foveal coords if the size is a function handle.

        rsize = round(roi_size( sqrt( (i-ref_loc(2)).^2+(j-ref_loc(1)).^2))); %Update the estimated size of our roi on every loop.
        
        roi_step = floor(rsize/rdivider);
        roi = cell(round((size(test_image)-rsize)/roi_step));

        if (rsize*0.8) < imbox(3) && (rsize*0.8) < imbox(4)           
            rsize = min(imbox(3:4))-1;
        else
            warning(['Size of valid region is too small for requested rsize (' num2str(rsize) 'x' num2str(rsize) ' vs ' num2str(imbox(3)) 'x' num2str(imbox(4)) ').'])
            sum_map=[];
            interped_spac_map=[];
            interped_conf_map=[];
            avg_pixel_spac=NaN;
            imbox=[];
            return;
        end
        
        roi_step = 0; % First loop we don't want to move, so make the step 0
        while (i+rsize-1)<(imbox(2)+imbox(4))
            while (j+rsize-1)<(imbox(1)+imbox(3))
                
                j=round(j+roi_step);
                rsize = round(roi_size( sqrt( (i-ref_loc(2)).^2+(j-ref_loc(1)).^2))); %Update the estimated size of our roi on every loop.
                
                if rsize > imbox(3) || rsize > imbox(4) %If our size is larger than one of the image sides, crop it.
                    rsize = min(imbox(3:4))-1;
                end
                if (i+rsize-1) > im_size(1)
                    rsize = im_size(1)-i+1;
                end
                if (j+rsize-1) > im_size(2)
                    rsize = im_size(2)-j+1;
                end

                roi_step = floor(rsize/rdivider);                

                numzeros = sum(sum(test_image(i:(i+rsize-1), j:(j+rsize-1))<=10));
                if numzeros < (rsize*rsize)*0.05
                    roi{round(i/roi_step)+1,round(j/roi_step)+1} = test_image(i:i+rsize-1, j:j+rsize-1);
                else
                    roi{round(i/roi_step)+1,round(j/roi_step)+1} =[];
                end

            end
            
            j=round(imbox(1)); % reset the column index to 0
            i=round(i+roi_step); %update our row position
        end
    end
end

clear test_image

numind = size(roi,1)*size(roi,2);
par_pixel_spac = nan(size(roi));
par_confidence = nan(size(roi));


goodr = find(~cellfun(@isempty,roi))';
length(goodr)
% If we're not running inside of a pool already, make one to speed this up.
if ~isempty(getCurrentTask())
    for i=1:length(goodr)
        
        r=goodr(i);
    
        power_spect = fftshift(fft2( double(roi{r}) ));
        power_spect = log10(abs(power_spect).^2);
        rhostart=1; % Exclude the DC term from our radial average
        
    %         figure(100); imagesc(roi{r}); axis image;
    %         power_spect_export = power_spect-min(power_spect(:));
    %         power_spect_export = power_spect_export./max(power_spect_export(:));
    %         power_spect_export = power_spect_export.*255;
    % 
    %         imwrite(uint8(power_spect_export),['pwr_spect ' num2str(r) '.tif']);
    
        rhosampling = .5;
        thetasampling = 1;
        
        [polarroi, power_spect_radius] = imcart2pseudopolar(power_spect, rhosampling, thetasampling,[],'makima' , rhostart);
        polarroi = circshift(polarroi,-90/thetasampling,1);
        
        %figure(101); imagesc(polarroi); axis image;
        
        % This really isn't needed because of Hermitian symmetry...
        %upper_n_lower = [thetasampling:45 136:225 316:360]/thetasampling;
        left_n_right = [46:135 226:315]/thetasampling;
        %upper_n_lower_fourierProfile = mean(polarroi(upper_n_lower,:)); %
        % Currently ignored, but could be used
        left_n_right_fourierProfile = mean(polarroi(left_n_right,:));
        
    
        [par_pixel_spac(i), ~, par_confidence(i)] = fourierFit(left_n_right_fourierProfile,[], false);
        par_pixel_spac(i) = 1/ (par_pixel_spac(i) / ((power_spect_radius*2)/rhosampling));
        
        
    end
else
    parfor i=1:length(goodr)
        
        r=goodr(i);
    
        power_spect = fftshift(fft2( double(roi{r}) ));
        power_spect = log10(abs(power_spect).^2);
        rhostart=1; % Exclude the DC term from our radial average
        
        rhosampling = .5;
        thetasampling = 1;
        
        [polarroi, power_spect_radius] = imcart2pseudopolar(power_spect, rhosampling, thetasampling,[],'makima' , rhostart);
        polarroi = circshift(polarroi,-90/thetasampling,1);
        
        %figure(101); imagesc(polarroi); axis image;
        
        % This really isn't needed because of Hermitian symmetry...
        %upper_n_lower = [thetasampling:45 136:225 316:360]/thetasampling;
        left_n_right = [46:135 226:315]/thetasampling;
        %upper_n_lower_fourierProfile = mean(polarroi(upper_n_lower,:)); %
        % Currently ignored, but could be used
        left_n_right_fourierProfile = mean(polarroi(left_n_right,:));
        
    
        [par_pixel_spac(i), ~, par_confidence(i)] = fourierFit(left_n_right_fourierProfile,[], false);
        par_pixel_spac(i) = 1/ (par_pixel_spac(i) / ((power_spect_radius*2)/rhosampling));
        
        
    end
end

pixel_spac = nan(size(roi));
confidence = nan(size(roi));

pixel_spac(goodr) = par_pixel_spac(1:length(goodr));
confidence(goodr) = par_confidence(1:length(goodr));

avg_pixel_spac = mean(pixel_spac(~isnan(pixel_spac)) );
interped_spac_map = avg_pixel_spac;
interped_conf_map = confidence;


%% If we've sampled over the region, then create the heat map
if length(roi) > 1
    
    interped_spac_map=zeros(im_size);
    interped_conf_map=zeros(im_size);
    sum_map=zeros(im_size);
    
    if ~isa(roi_size, 'function_handle')
        
        roi_coverage=roi_size; 
        
        for i=imbox(2):roi_step:imbox(2)+imbox(4)-roi_size
            for j=imbox(1):roi_step:imbox(1)+imbox(3)-roi_size

                if ~isnan( pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1) )

                    thiserr = confidence(round(i/roi_step)+1,round(j/roi_step)+1)^2;

                    interped_conf_map(i:i+roi_coverage-1, j:j+roi_coverage-1) = interped_conf_map(i:i+roi_coverage-1, j:j+roi_coverage-1) +(thiserr); 
                    thisspac = pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1);


                    interped_spac_map(i:i+roi_coverage-1, j:j+roi_coverage-1) = interped_spac_map(i:i+roi_coverage-1, j:j+roi_coverage-1) + ((thiserr*thisspac));

                    sum_map(i:i+roi_coverage-1, j:j+roi_coverage-1) = sum_map(i:i+roi_coverage-1, j:j+roi_coverage-1) + 1;

                else

                end
            end
        end

    else %If roi_size is a function, that means we'll expect it to vary.
       
        j=round(imbox(1));
        i=round(imbox(2)); %ref_loc takes on the foveal coords if the size is a function handle.
        rsize = round(roi_size( sqrt( (i-ref_loc(2)).^2+(j-ref_loc(1)).^2))); %Update the estimated size of our roi on every loop.

        if (rsize*0.8) < imbox(3) && (rsize*0.8) < imbox(4)           
            rsize = min(imbox(3:4))-1;
        else
            warning(['Size of valid region is too small for requested rsize (' num2str(rsize) 'x' num2str(rsize) ' vs ' num2str(imbox(3)) 'x' num2str(imbox(4)) ').'])
            sum_map=[];
            interped_spac_map=[];
            interped_conf_map=[];
            imbox=[];
            avg_pixel_spac=NaN;
            return;
        end

        roi_step = 0;
             
        while (i+rsize-1)<(imbox(2)+imbox(4))
            while (j+rsize-1)<(imbox(1)+imbox(3))
                                         
                j=round(j+roi_step);
                rsize = round(roi_size( sqrt( (i-ref_loc(2)).^2+(j-ref_loc(1)).^2))); %Update the estimated size of our roi on every loop.
                
                if rsize > imbox(3) || rsize > imbox(4) %If our size is larger than one of the image sides, crop it.
                    rsize = min(imbox(3:4))-1;
                end
                if (i+rsize-1) > im_size(1)
                    rsize = im_size(1)-i+1;
                end
                if (j+rsize-1) > im_size(2)
                    rsize = im_size(2)-j+1;
                end
                roi_step = floor(rsize/rdivider);


                if ~isnan( pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1) )


                    thiserr = confidence(round(i/roi_step)+1,round(j/roi_step)+1)^2;                    
                    
                    interped_conf_map(i:i+rsize-1, j:j+rsize-1) = interped_conf_map(i:i+rsize-1, j:j+rsize-1) +thiserr;    
                    thisspac = pixel_spac(round(i/roi_step)+1,round(j/roi_step)+1);


                    interped_spac_map(i:i+rsize-1, j:j+rsize-1) = interped_spac_map(i:i+rsize-1, j:j+rsize-1) + thiserr*thisspac;

                    sum_map(i:i+rsize-1, j:j+rsize-1) = sum_map(i:i+rsize-1, j:j+rsize-1) + 1;

                end
                
            end
            j=round(imbox(1)); % reset the column index to "0"
            i=round(i+roi_step); %update our row position
        end
    end

    interped_spac_map = interped_spac_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    interped_conf_map = interped_conf_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    sum_map = sum_map( imbox(2):imbox(2)+imbox(4), imbox(1):imbox(1)+imbox(3) );
    
    if nargout <= 1

       figure(1);clf; imagesc(interped_spac_map./interped_conf_map); axis image;    
        figure(2);clf; imagesc((interped_conf_map./sum_map)); colormap hot;
        figure(3);clf; imagesc(sum_map); axis image; colormap gray;
    end
end
  
end