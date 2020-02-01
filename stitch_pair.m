function [im_stitched, im_stitched_ref, im1_mod, im2_mod, im1_ref, im2_ref, im1_D, im2_D] = stitch_pair(im1,im2,varargin)
% Function to stitch two images together.

% An approximate alignment is first found by matching SURF features. Next, a
% local displacement field is calculated using imregdemons to align the
% second image onto the first image. This displacement is then interpolated
% in the intersecting region such that the displacement field is maximal near the
% boundaries that touch the second image and vanishing near the boundaries
% that touch the first image. Using this interpolated displacement field, a
% modified version of image 2 is created. Finally, a local displacement field is found
% to align the first image with the modified version of image 2.

% Input
% im1 = double array representing the first image
% im2 = double array representing the second image

% Optional call
% stitch_pair(im1,im2,im1_ref,im2_ref,...)
% im1_ref = imref2D object referencing the first image
% im2_ref = imref2D object referencing the second image

% Output
% im_stitched = stitched image.
% im_stitched_ref = imref2D object referencing the stitched image
% im1_mod = modified version of the first image that was stitched
% im2_mod = modified version of the second image that was stitched
% im1_ref = imref2D object referencing the first image
% im2_ref = imref2D object referencing the second image
% im1_D = displacement field used to align the first image
% im2_D = displacement field used to align the second image
%% Parse inputs
if ~isempty(varargin) && isa(varargin{1},'imref2d')
    im1_ref = varargin{1};
    varargin(1) = [];
else
    im1_ref = imref2d(size(im1));
end
%% Parse optional parameters
p = inputParser;
addParameter(p, 'PlotStitchedImage', nargout==0); % Plot the stitched image.
addParameter(p, 'Plots', false); % Plot intermediate steps of the stitching process.
addParameter(p, 'AdjustContrast', true); % Compute displacement fields on contrast-adjused images. The final stitched image preserves the intensity range of the input images.
addParameter(p, 'FieldSmoothing', 1.5); % 'AccumulatedFieldSmoothing' parameter used in imregdemons.
addParameter(p, 'Waitbar', true); % Display progress waitbars.
parse(p, varargin{:});
options = p.Results;
%% Adjust the intensity to saturate at least 1% for both images.
if options.AdjustContrast
    saturation_level = 0.01;
    im1_threshold = stretchlim(im1,saturation_level);
    im2_threshold = stretchlim(im2,saturation_level);
    threshold = max(im1_threshold,im2_threshold);
    
    im1_adj = imadjust(im1,threshold);
    im2_adj = imadjust(im2,threshold);
else
    im1_adj = im1;
    im2_adj = im2;
end

if options.Plots
    figure;imshowpair(im1_adj,im2_adj,'montage');tightax;title('Input Images');
    title('Image 1 and Image 2');tightax;
end
%% Estimate the initial translation by matching SURF features.
im1_points = detectSURFFeatures(im1_adj);    
[im1_features, im1_points] = extractFeatures(im1_adj, im1_points);

im2_points = detectSURFFeatures(im2_adj);    
[im2_features, im2_points] = extractFeatures(im2_adj, im2_points);

% Match features.
indexPairs = matchFeatures(im1_features,im2_features, 'Unique', true);

im1_matched_points = im1_points(indexPairs(:,1), :);
im2_matched_points = im2_points(indexPairs(:,2), :);

% Estimate the transformation.
im2_translation_estimate = estimateGeometricTransform(im2_matched_points,im1_matched_points,'similarity', 'Confidence', 99.9, 'MaxNumTrials', 3000);

% Round the translation to pixel size.
im2_translation_estimate.T = round(im2_translation_estimate.T);
%% Computer a better translation estimate with registration (minimizes Mattes Mutual Information).
[optimizer, metric] = imregconfig('multimodal');
%optimizer.InitialRadius = 1e-4;
%optimizer.MaximumStepLength = 1e-3;

% imregtform can sometimes not converge. In this case, imregtform will return the initial transformation.
warning('off')
im2_translation = imregtform(im2_adj, im1_adj, 'translation', optimizer, metric,'InitialTransformation',im2_translation_estimate);
warning('on')

% Round the translation to a pixel.
im2_translation.T = round(im2_translation.T);
%% Translate im2.
% Align the second image with the first image.
im2_ref = imref2d(size(im2_adj));
im2_ref.XWorldLimits = im1_ref.XWorldLimits(1) + [0 diff(im2_ref.XWorldLimits)];
im2_ref.YWorldLimits = im1_ref.YWorldLimits(1) + [0 diff(im2_ref.YWorldLimits)];
[im2_translated,im2_ref] = imwarp(im2_adj,im2_ref,im2_translation);

if options.Plots
    figure; imshowpair(im1_adj,im1_ref,im2_translated,im2_ref,'scaling','joint');
    title('Image 1(Green) and translated Image 2(Magenta)'); tightax;
end
%% Find the region intersecting im1 and im2
% Find the position of the intersecting region in world coordinates for each image.
Intersection_lim_world = [max(im2_ref.XWorldLimits(1),im1_ref.XWorldLimits(1)) min(im2_ref.XWorldLimits(2),im1_ref.XWorldLimits(2)) max(im2_ref.YWorldLimits(1),im1_ref.YWorldLimits(1)) min(im2_ref.YWorldLimits(2),im1_ref.YWorldLimits(2))];

[intersect1_lim(1:2), intersect1_lim(3:4)] = worldToIntrinsic(im1_ref,Intersection_lim_world(1:2),Intersection_lim_world(3:4));
intersect1_lim = round(intersect1_lim + 0.5*[1 -1 1 -1]);

[intersect2_lim(1:2), intersect2_lim(3:4)] = worldToIntrinsic(im2_ref,Intersection_lim_world(1:2),Intersection_lim_world(3:4));
intersect2_lim = round(intersect2_lim + 0.5*[1 -1 1 -1]);

%% Crop the input images to retain only the intersecting regions.
%Intersection1_rect = [Intersection1_lim_intrinsic([1 3]) Intersection1_lim_intrinsic([2 4])-Intersection1_lim_intrinsic([1 3])];
%[im1_intersect,rect2] = imcrop(im1_a,Intersection1_rect);
im1_intersect = im1_adj(intersect1_lim(3):intersect1_lim(4),intersect1_lim(1):intersect1_lim(2));

%Intersection2_rect = [Intersection2_lim_intrinsic([1 3]) Intersection2_lim_intrinsic([2 4])-Intersection2_lim_intrinsic([1 3])];
%[im2_intersect,rect2] = imcrop(im2_translated,Intersection2_rect);
im2_intersect = im2_translated(intersect2_lim(3):intersect2_lim(4),intersect2_lim(1):intersect2_lim(2));
%% Find empty regions in im1 that intersect with non-empty regions of im2.
% Im1 may have empty patches(==0) often coming from previous stitching. 
% If some of these empty patches intersect with non-empty regions of im2, 
% the calculated displacement field washes away the non-empty pixels of im2 to match the empty region of im1.
% This is unwanted as non-empty regions of im2 correspond to relevant pixel intensities that should be kept.
% To resolve this issue, the displacement field that moves im2 on im1 is
% set to zero in intersecting regions where im1 has zero pixel values and im2 has non-zero pixel values.

% Find the positions of zero pixels in im1.
im1_iszero = im1_adj==0;

% Remove isolated zero pixels.
im1_iszero = bwmorph(bwmorph(im1_iszero,'erode',2),'dilate',2);

% Find the empty patches.
im1_emptypatches = bwconncomp(im1_iszero);

% Remove extremely small patches.
N_pixels_min = ceil(numel(im1_adj)*0.00001);
im1_emptypatches.PixelIdxList = im1_emptypatches.PixelIdxList(cellfun(@(x) numel(x)>N_pixels_min,im1_emptypatches.PixelIdxList));
im1_emptypatches.NumObjects = numel(im1_emptypatches.PixelIdxList);

% Define a mask that identifies the pixels belonging to empty patches.
im1_has_emptypatches = im1_emptypatches.NumObjects > 0;
if im1_has_emptypatches
    im1_emptypatches_logind = false(size(im1_adj));
    im1_emptypatches_props = regionprops(im1_emptypatches);
    for i=1:im1_emptypatches.NumObjects
        BoundingBox(1:2) = ceil(im1_emptypatches_props(i).BoundingBox(1:2));
        BoundingBox(3:4) = BoundingBox(1:2) + im1_emptypatches_props(i).BoundingBox(3:4)-1;
        
        im1_emptypatches_logind(BoundingBox(2):BoundingBox(4),BoundingBox(1):BoundingBox(3)) = true;
        im1_emptypatches_props(i).BoundingBox2 = BoundingBox;
    end
    
    % Restrict the empty patches logical array to the intersecting region
    % in im1.
    im1_emptypatches_intersect = im1_emptypatches_logind(intersect1_lim(3):intersect1_lim(4),intersect1_lim(1):intersect1_lim(2));
end
%% Pad the intersecting region to avoid drifting displacement fields that produce zero pixels at the boundary.
% To pad the arrays pixels from im1 or im2 are used. Pixels for im1 are
% used where the intersecting boundary touches im1 (and vice versa for
% im2).

% First, determine the type of the intersecting boundary. 
% type = 1 indicates that the boundary is contiguous to im1.
% type = 2 indicates that the boundary is contiguous to im2.
top_boundary_type = 2 - double(im1_ref.YWorldLimits(1) < im2_ref.YWorldLimits(1));
bottom_boundary_type = 2 - double(im1_ref.YWorldLimits(2) > im2_ref.YWorldLimits(2));
left_boundary_type = 2 - double(im1_ref.XWorldLimits(1) < im2_ref.XWorldLimits(1));
right_boundary_type = 2 - double(im1_ref.XWorldLimits(2) > im2_ref.XWorldLimits(2));

% Pad the intersection arrays.
pad_size = 15;
%pad_size = 20;

padding_params.pad_size = pad_size;
padding_params.left_pad_type = left_boundary_type;
padding_params.bottom_pad_type = bottom_boundary_type;
padding_params.right_pad_type = right_boundary_type;
padding_params.top_pad_type = top_boundary_type;
padding_params.ims = {im1_adj im2_translated};
padding_params.ims_intersect_lim = {intersect1_lim intersect2_lim};

[im1_intersect_padded,im2_intersect_padded] = pad_intersections(im1_intersect,im2_intersect,padding_params,'Method','neighbours');

% Pad the empty patches logical array of im1.
if im1_has_emptypatches
    im1_emptypatches_intersect_padded = padarray(im1_emptypatches_intersect,pad_size*ones(1,2),'replicate');
end
%% Register im2 on im1.
% Register the padded arrays. Padding helps to remove washing artifacts
% that are often seen at the boundaries.

[im2_D_full, im2_intersect_padded_displaced] = ...
    imregdemons(im2_intersect_padded,im1_intersect_padded,[500 400 300 200 200],...
    'PyramidLevels',5,'AccumulatedFieldSmoothing',options.FieldSmoothing,'DisplayWaitbar',options.Waitbar);

% [im2_D_full, im2_intersect_padded_displaced] = ...
%     imregdemons(im2_intersect_padded,im1_intersect_padded,[200 200 200 100 100 50],...
%     'PyramidLevels',6,'AccumulatedFieldSmoothing',options.FieldSmoothing,'DisplayWaitbar',options.Waitbar);

% [im2_D_full, im2_intersect_padded_displaced] = ...
%     imregdemons(im2_intersect_padded,im1_intersect_padded,[400 300 200 100 50],...
%     'PyramidLevels',5,'AccumulatedFieldSmoothing',options.FieldSmoothing,'DisplayWaitbar',options.Waitbar);

% Register on binary image (not good). Too much speckle noise.
% im1_intersect_bin = im1_intersect > 110;
% im2_intersect_bin = im2_intersect > 110;
% [im2_D_ref_bin] = imregdemons(im2_intersect_bin,im1_intersect_bin,[500 400 300 200 200],'PyramidLevels',5,'AccumulatedFieldSmoothing',1);
% im2_intersect_ref_bin = imwarp(im2_intersect,im2_D_ref_bin);

if options.Plots
    figure;imshowpair(im2_intersect_padded,im2_intersect_padded_displaced)
    title('Image 2 (Green) and Image 2 displaced (Magenta) (Padded Images)');tightax
    
    figure;imshowpair(im1_intersect_padded,im2_intersect_padded_displaced);
    title('Image 1(Green) and fully-displaced Image2(Magenta) (Padded Images)');tightax;
end
%% Calculate the displacement field interpolator in the intersecting region.
% The interpolation is such that the displacement field is maximum near im1 and minimum near im2.

% Define a mask that identifies the boundaries type in the intersection
% region.
interpolator_size = size(im2_intersect_padded);
boundary_mask = zeros(interpolator_size);

boundary_mask(1,:) = top_boundary_type;
boundary_mask(end,:) = bottom_boundary_type;
boundary_mask(:,1) = left_boundary_type;
boundary_mask(:,end) = right_boundary_type;
if im1_has_emptypatches
    boundary_mask(im1_emptypatches_intersect_padded) = 2;
end

% Calculate the interpolator.
t_interp = calculate_interpolator(boundary_mask);

% Scale the interpolator so that it equals 1 near the im1 boundaries (type=1) and 0
% near the im2 (type=2) boundaries.
t_interp = 2 - t_interp;

% Plot the interpolator.
if 0 && options.Plots
    figure;imagesc(t_interp);colorbar;
    title('Image 2 Displacement field Interpolator (1=Image1, 0=Image2)');
    ylabel('row index');
    xlabel('col index');
end
%% Interpolate the displacement field that moves im2 on im1.
% Interpolate the displacement field of im2. Close to the borders
% contiguous to im1, t_interp ~= 1, i.e., the displacement on im2
% is almost fully applied. On the other hand, close to the borders
% contiguous to im2, t_interp ~= 0, i.e., the displacement on im2
% is not applied.
im2_D = repmat(t_interp,[1 1 2]).*im2_D_full;

% Move im2 with the interpolated displacement field.
im2_intersect_displaced = imwarp(im2_intersect_padded,im2_D);

% Plot the displaced second image.
if options.Plots
    figure;imshowpair(im2_intersect_padded,im2_intersect_displaced);
    title('Image 2 (Green) and partially-displaced Image2 (Magenta)');tightax
end
%% Compute the displacement field that moves im1 on the displaced im2.
[im1_D,im1_intersect_padded_displaced] = ...
    imregdemons(im1_intersect_padded,im2_intersect_displaced,[500 400 300 200 200],...
    'PyramidLevels',5,'AccumulatedFieldSmoothing',options.FieldSmoothing,'DisplayWaitbar',options.Waitbar);
%[im1_D,im1_intersect2] = imregdemons(im1_intersect,im_intersect_ref,[500 400 300 200 200],'PyramidLevels',5,'AccumulatedFieldSmoothing',0.7);
%[im2_D,im2_intersect2] = imregdemons(im2_intersect,im_intersect_ref,[500 400 300 200 200],'PyramidLevels',5,'AccumulatedFieldSmoothing',0.7);

if options.Plots
    figure;imshowpair(im1_intersect_padded_displaced,im2_intersect_displaced);
    title('Displaced Image 1 (Green) and partially-displaced Image2 (Magenta)');tightax;
end
%% Apply the computed transformations on the original images.
% Translate im2.
[im2_translated,im2_ref] = imwarp(im2,im2_translation);

% Get the intersections.
im1_intersect = im1(intersect1_lim(3):intersect1_lim(4),intersect1_lim(1):intersect1_lim(2));
im2_intersect = im2_translated(intersect2_lim(3):intersect2_lim(4),intersect2_lim(1):intersect2_lim(2));

% Pad the intersections.
padding_params.ims = {im1 im2_translated};
[im1_intersect_padded,im2_intersect_padded] = pad_intersections(im1_intersect,im2_intersect,padding_params,'Method','neighbours');

% Warp the images with their respective displacement field.
im1_intersect_displaced = imwarp(im1_intersect_padded,im1_D);
im2_intersect_displaced = imwarp(im2_intersect_padded,im2_D);

% Linearly blend the pixel intensities of two displaced images to get the final
% intersecting region.
% The same interpolator that was used on the image 2 displacement field is used again.
im_class_func = str2func(class(im1));
im_intersect_final = im_class_func(t_interp.*double(im1_intersect_displaced) + (1-t_interp).*double(im2_intersect_displaced));

% Remove pad.
im_intersect_final = im_intersect_final(pad_size+1:end-pad_size,pad_size+1:end-pad_size);
%% Overwrite the intersecting region in each image with the intersecting patch calculated above.
im1_mod = im1;
im1_mod(intersect1_lim(3):intersect1_lim(4),intersect1_lim(1):intersect1_lim(2)) = im_intersect_final;

im2_mod = im2_translated;
im2_mod(intersect2_lim(3):intersect2_lim(4),intersect2_lim(1):intersect2_lim(2)) = im_intersect_final;

% Plot the final result of the stitching.
if options.PlotStitchedImage || options.Plots
    figure;imshowpair(im1_mod,im1_ref,im2_mod,im2_ref,'scaling','joint');
    title('Stiched Image 1(Green) and Image 2(Magenta)')
end
%% Fuse both images to produce the final stiched image.
XWorldLimits = [im1_ref.XWorldLimits; im2_ref.XWorldLimits];
YWorldLimits = [im1_ref.YWorldLimits; im2_ref.YWorldLimits];
im_stitched_XWorldLimits = [min(XWorldLimits(:,1)) max(XWorldLimits(:,2))];
im_stitched_YWorldLimits = [min(YWorldLimits(:,1)) max(YWorldLimits(:,2))];
im_stitched_size = round([diff(im_stitched_YWorldLimits) diff(im_stitched_XWorldLimits)]);
im_stitched_origins = [repmat(im_stitched_XWorldLimits(1),1,2) repmat(im_stitched_YWorldLimits(1),1,2)];

im1_pos = round(([XWorldLimits(1,:) YWorldLimits(1,:)] - im_stitched_origins + [1 0 1 0]));
im2_pos = round(([XWorldLimits(2,:) YWorldLimits(2,:)] - im_stitched_origins + [1 0 1 0]));

im_stitched = zeros(im_stitched_size,class(im1_mod));
im_stitched(im1_pos(3):im1_pos(4),im1_pos(1):im1_pos(2)) = im1_mod;
im_stitched(im2_pos(3):im2_pos(4),im2_pos(1):im2_pos(2)) = im2_mod;

im_stitched_ref = imref2d(im_stitched_size);

% Using imfuse changes the data format of the input image.
%[im_stitched, im_stitched_ref] = imfuse(im1_mod,im1_ref,im2_mod,im2_ref,'blend','Scaling','joint');
end
%% Function to pad images in the intersecting region.
function [im1_intersect_padded, im2_intersect_padded] = pad_intersections(im1_intersect,im2_intersect,params,varargin)

% Parse optional parameters
p = inputParser;
addParameter(p, 'Method', 'neighbours');
parse(p, varargin{:});
options = p.Results;
pad_size = params.pad_size;

% Pad images.
switch options.Method
    case 'replicate'
        im1_intersect_padded = padarray(im1_intersect,pad_size*ones(1,2),'replicate');
        im2_intersect_padded = padarray(im2_intersect,pad_size*ones(1,2),'replicate');
    case 'zero'
        im1_intersect_padded = padarray(im1_intersect,pad_size*ones(1,2),0);
        im2_intersect_padded = padarray(im2_intersect,pad_size*ones(1,2),0);
    case 'neighbours'
        ims = params.ims;
        ims_size = cellfun(@size,ims,'Uni',0);
        ims_intersect_lim = params.ims_intersect_lim;
        
        % Top pad.
        pad_type = params.top_pad_type;
        in_lim = ims_intersect_lim{pad_type};
        top_center_band_height = in_lim(3) - max((in_lim(3) - pad_size),1);
        top_center_band = ims{pad_type}(in_lim(3) - (top_center_band_height:-1:1), in_lim(1):in_lim(2));
        top_pad = [zeros(pad_size) [zeros(pad_size - top_center_band_height,size(top_center_band,2)); top_center_band] zeros(pad_size)];
        
        % Bottom pad.
        pad_type = params.bottom_pad_type;
        in_lim = ims_intersect_lim{pad_type};
        bottom_center_band_height = min(in_lim(4) + pad_size,ims_size{pad_type}(1)) - in_lim(4);
        bottom_center_band = ims{pad_type}(in_lim(4) + (1:bottom_center_band_height), in_lim(1):in_lim(2));
        bottom_pad = [zeros(pad_size) [bottom_center_band;zeros(pad_size-bottom_center_band_height,size(bottom_center_band,2))] zeros(pad_size)];
        
        % Left pad.
        pad_type = params.left_pad_type;
        in_lim = ims_intersect_lim{pad_type};
        left_center_band_width = in_lim(1) - max((in_lim(1) - pad_size),1);
        left_center_band = ims{pad_type}(in_lim(3):in_lim(4), in_lim(1) - (left_center_band_width:-1:1));
        left_pad = [zeros(pad_size); [zeros(size(left_center_band,1),pad_size - left_center_band_width) left_center_band] ; zeros(pad_size)];
        
        % Right pad.
        pad_type = params.right_pad_type;
        in_lim = ims_intersect_lim{pad_type};
        right_center_band_width = min(in_lim(2) +  pad_size,ims_size{pad_type}(2)) - in_lim(2);
        right_center_band = ims{pad_type}(in_lim(3):in_lim(4), in_lim(2) + (1:right_center_band_width));
        right_pad = [zeros(pad_size); [right_center_band zeros(size(right_center_band,1),pad_size-right_center_band_width)] ; zeros(pad_size)];
        
        % Create a pad template for both images.
        pad_template = zeros(size(im1_intersect) + 2*pad_size);
        pad_template(1:pad_size,:) = top_pad;
        pad_template(end-pad_size+1:end,:) = bottom_pad;
        pad_template(:,1:pad_size) = left_pad;
        pad_template(:,end-pad_size+1:end) = right_pad;
        
        % Pad arrays
        im1_intersect_padded = pad_template;
        im1_intersect_padded(pad_size+1:end-pad_size,pad_size+1:end-pad_size) = im1_intersect;
        im2_intersect_padded = pad_template;
        im2_intersect_padded(pad_size+1:end-pad_size,pad_size+1:end-pad_size) = im2_intersect;
    otherwise
        error('Padding method is undefined.');
end
end