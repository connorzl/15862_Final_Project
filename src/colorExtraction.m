close all
clear all

%% read in image
img_name = '../data/peach.png';
img = im2double(imread(img_name));
[imh, imw, ~] = size(img);

canvasScale = 2;
numRows = imh * canvasScale;
numCols = imw * canvasScale;

img_large = imresize(img, canvasScale);

layers = load('edgeclip_layers.mat');
layer0 = layers.layer0;
layer1 = layers.layer1;
layer2 = layers.layer2;
layer3 = layers.layer3;

[layer0,mask0] = add_color_to_strokes(layer0,numRows,numCols,img_large);
[layer1,mask1] = add_color_to_strokes(layer1,numRows,numCols,img_large);
[layer2,mask2] = add_color_to_strokes(layer2,numRows,numCols,img_large);
[layer3,mask3] = add_color_to_strokes(layer3,numRows,numCols,img_large);

save('color_layers.mat','layer0','layer1','layer2','layer3');

function [layer,mask] = add_color_to_strokes(layer,numRows,numCols,img_large)
mask = zeros(numRows,numCols,3);
for i = 1:size(layer,1)
    curr_stroke = layer(i);
    if size(curr_stroke.stroke_pixels,1) == 0
        continue
    end
    curr_color = [0, 0, 0];
    
    for j = 1:size(curr_stroke.stroke_pixels,1)
        r = curr_stroke.stroke_pixels(j,1);
        c = curr_stroke.stroke_pixels(j,2);
        pixel_color = squeeze(img_large(r,c,:));
        curr_color = curr_color + pixel_color';
    end
    if mod(i,50) == 0
        disp(i);
    end
    
    curr_stroke.color = curr_color / size(curr_stroke.stroke_pixels,1);
    for j = 1:size(curr_stroke.stroke_pixels,1)
        r = curr_stroke.stroke_pixels(j,1);
        c = curr_stroke.stroke_pixels(j,2);
        mask(r,c,:) = curr_stroke.color';
    end
    layer(i) = curr_stroke;
end
end