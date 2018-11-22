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

%% find strong strokes
img_large = imresize(img, canvasScale);
img_grayscale = rgb2gray(img_large);

width = wb;
kernelSize = [width width];
kernel = fspecial('gaussian',kernelSize);
img_blur = imfilter(img_grayscale,kernel,'same');

[Gx,Gy] = imgradientxy(img_blur,'sobel');
thresh = 1;

numStrong = 0;
for i = 1:size(baseLayer)
    S = baseLayer(i);
    v = [Gx(S.c), Gy(S.r)];
    if norm(v) >= thresh
        numStrong = numStrong + 1;
        S.strong = 1;
        baseLayer(i) = S;
    end
end

baseCanvas = zeros(numRows, numCols);
for i = 1:size(baseLayer)
    S = baseLayer(i);
    if S.strong
        baseCanvas(S.r,S.c) = 1;
    end
end
imshow(baseCanvas);
