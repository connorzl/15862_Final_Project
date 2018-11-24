%% read in image
img = im2double(imread('../data/peach.png'));
[imh, imw, ~] = size(img);

canvasScale = 2;
numRows = imh * canvasScale;
numCols = imw * canvasScale;

% brush width
wb = 36;

layers = load('layers.mat');
baseLayer = layers.baseLayer;
layer1 = layers.layer1;
layer2 = layers.layer2;
layer3 = layers.layer3;

% scale image to canvas size and convert to grayscale
img_large = imresize(img, canvasScale);
img_grayscale = rgb2gray(img_large);

% find strong strokes for base layer
thresh = 0.1;
[canvas, baseLayer] = findStrongStrokes(baseLayer, numRows, numCols, wb,...
    thresh, img_grayscale);
%figure;
%imshow(canvas);

% find strong strokes for other layers
[canvas, layer1] = findStrongStrokes(layer1, numRows, numCols, wb/2,...
    thresh, img_grayscale);
figure;
imshow(canvas);

%% find strong strokes
function [canvas, layer] = findStrongStrokes(layer, numRows, numCols, width,...
    thresh, img_grayscale)

kernelSize = [width width];
kernel = fspecial('gaussian',kernelSize);
img_blur = imfilter(img_grayscale,kernel,'same');
[Gx,Gy] = imgradientxy(img_blur,'sobel');

numStrong = 0;
for i = 1:size(layer)
    S = layer(i);
    v = [Gx(S.r,S.c), Gy(S.r,S.c)];
    if norm(v) >= thresh
        numStrong = numStrong + 1;
        S.strong = 1;
        layer(i) = S;
    end
end

canvas = zeros(numRows, numCols);
for i = 1:size(layer)
    S = layer(i);
    if S.strong
        canvas(S.r,S.c) = 1;
    end
end

disp(numStrong);
end
