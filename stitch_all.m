function [image_stitched, images_info] = stitch_all(images,varargin)
% Function to stitch several images using built-in methods of the image processing toolbox.
% The algorithm first aligns image by matching SURF features
% and fine-tunes the alignment by finding a local displacement field necessary for alignment. 
% The local field is found using the built-in function 'imregdemons'.

% Input
% images = cell array of images. The images are assumed to be matrices of
% various type (uint8,uint16,uint32,double)

% Output
% image_stitched = fused image containing all input images.
% images_info = structure containing information about stitching and
% preprocessing. The field 'im_ref' contains an imref2D
% object that can be used to locate each input image within the fused image.
%%
if ~iscell(images)
    error('The input images must be in a cell array')
end
N_images = numel(images);
%% Parse optional parameters
p = inputParser;
addParameter(p, 'PlotStitchedImage', nargout==0); % Plot the stitched image.
addParameter(p, 'Waitbar', nargout==0); % Display stitching progress waitbars.
addParameter(p, 'Plots', false); % Plot intermediate images detailing the stitching process.
parse(p, varargin{:});
options = p.Results;
%% Create contrast-adjusted images.
% This is sometimes necessary in order to find SURF features.
images_info = struct();
saturation_level = 0.01;
for i=1:N_images
    images_info(i).contrast_thresholds = stretchlim(images{i},saturation_level);
end
contrast_thresholds_all = cell2mat({images_info.contrast_thresholds});
contrast_threshold = [min(contrast_thresholds_all(:,1)) max(contrast_thresholds_all(:,2))];

images_adjusted = cell(N_images,1);
for i=1:N_images
    images_adjusted{i} = imadjust(images{i},contrast_threshold);
end
%% Detect SURF features in contrast-adjusted images.
surf_features = struct();
for i=1:N_images
    im_points = detectSURFFeatures(images_adjusted{i});
    [surf_features(i).features, surf_features(i).points] = extractFeatures(images_adjusted{i}, im_points);
end
%% Match images based on SURF features.
% For each image, find another image that has the most matched SURF features.
for i=1:N_images
    image_inds = setdiff(1:N_images,i);
    surf_features(i).MatchFeatures_ind = cell(1,N_images);
    
    for j=1:numel(image_inds)
        image_ind = image_inds(j);
        surf_features(i).MatchFeatures_ind{image_ind} = ...
            matchFeatures(surf_features(i).features,surf_features(image_ind).features, 'Unique', true);
    end    
end
%% Determine the order in which images will be asembled.
% Find the image that has the most matched features with other images.
for i=1:N_images
    surf_features(i).N_matched_features = cellfun(@numel,surf_features(i).MatchFeatures_ind);
    [~,surf_features(i).best_match_ind] = max(surf_features(i).N_matched_features);
end
N_matched_features_all = cell2mat({surf_features.N_matched_features}');
N_matched_features_sum = sum(N_matched_features_all,2);
[~,start_image_ind] = max(N_matched_features_sum);

% Initialize the assembling order.
assembling_order_ind = nan(1,N_images);
assembling_order_ind(1) = start_image_ind;
assembling_order_ind(2) = surf_features(start_image_ind).best_match_ind;

for j=3:N_images
    % Among the remaining images, find the image whose best match is
    % already stitched and which has the highest matched features with
    % images that are already matched.
    matched_images_ind = assembling_order_ind(1:j-1);
    unmatched_images_ind = setdiff(1:N_images,matched_images_ind);
    candidate_images_ind = find(ismember([surf_features.best_match_ind],matched_images_ind) & ~ismember(1:N_images,matched_images_ind));
    
    if isempty(candidate_images_ind)
        % In this case, find the image among the unmatched images that has
        % the most matched features with the matched images.
        [~,ind] = max(sum(N_matched_features_all(unmatched_images_ind,matched_images_ind),2));
        next_image_ind = unmatched_images_ind(ind);
    else
        [~,ind] = max(sum(N_matched_features_all(candidate_images_ind,matched_images_ind),2));
        next_image_ind = candidate_images_ind(ind);
    end
    
    assembling_order_ind(j) = next_image_ind;
end

% Check that all images are accounted for in the assembling order.
if numel(unique(assembling_order_ind)) ~= N_images || any(~ismember(assembling_order_ind,1:N_images))
    error('Some images are missing in the assembling order.');
end
%% Stich images
% Stitch the first two images using the first image as reference.
[image_stitched, image_stitched_ref, ~, ~, im1_ref, stitched_image_ref] = ...
    stitch_pair(images{assembling_order_ind(1)},images{assembling_order_ind(2)},'Waitbar',options.Waitbar,'Plots',options.Plots);
images_info(assembling_order_ind(2)).im_ref = stitched_image_ref;
images_info(assembling_order_ind(1)).im_ref = im1_ref;

% Sitch the rest of the images.
for j=3:N_images
    % Sitch the next image onto the previously-stitched image.
    image_to_stitch = images{assembling_order_ind(j)};
    [image_stitched, image_stitched_ref, ~, ~, ~, stitched_image_ref] = ...
        stitch_pair(image_stitched,image_to_stitch,image_stitched_ref,'Waitbar',options.Waitbar,'Plots',options.Plots);
    images_info(image_ind).im_ref = stitched_image_ref;
end
end