%% Read all images and store them in a cell array.
image_files = dir('./images/*.png');
N_images = numel(image_files);
images = cell(N_images,1);
for i=1:N_images
    Filename = fullfile(image_files(i).folder,image_files(i).name);
    images{i} = imread(Filename);
end
%% Stitch images.
image_stitched = stitch_all(images,'PlotStitchedImage',1);
stitched_image_filename = fullfile(image_files(1).folder,'Neuron-stitched.png');
imwrite(image_stitched,stitched_image_filename);