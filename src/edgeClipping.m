close all
clear all

%% read in image
img = im2double(imread('../data/peach.png'));
[imh, imw, ~] = size(img);

canvasScale = 2;
numRows = imh * canvasScale;
numCols = imw * canvasScale;

img_large = imresize(img, canvasScale);
img_grayscale = rgb2gray(img_large);
wb = 36;

layers = load('orientation_layers.mat');
layer0 = layers.layer0;
layer1 = layers.layer1;
layer2 = layers.layer2;
layer3 = layers.layer3;

%% blur layers
sigma_0 = round(0.3*wb);
sigma_1 = round(0.15*wb);
sigma_2 = round(0.07*wb);

img_0 = imgaussfilt(img_grayscale, sigma_0);
img_1 = imgaussfilt(img_grayscale, sigma_1);
img_2 = imgaussfilt(img_grayscale, sigma_2);
img_3 = img_grayscale;
%{
figure;
imshow(img_0);
figure;
imshow(img_1);
figure;
imshow(img_2);
figure;
imshow(img_3);
%}
edge_img_0 = edge(img_0,'Canny',0.4);
edge_img_1 = edge(img_1,'Canny',0.3);
edge_img_2 = edge(img_2,'Canny',0.2);
edge_img_3 = edge(img_3,'Canny',0.1);
%{
figure;
imshow(edge_img_0);
figure;
imshow(edge_img_1);
figure;
imshow(edge_img_2);
figure;
imshow(edge_img_3);
%}

for s = 1:size(layer0,1)
    curr_stroke = layer0(s);
    r = curr_stroke.r;
    c = curr_stroke.c;
    ang = curr_stroke.ang;
    
    dX = 0;
    dY = 0;
    while true
        % grow stroke both directions
        
        % for loop in a radius along stroke width to check if hit edge
        
        % set l accordingly
        
    end
end