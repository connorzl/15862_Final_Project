close all
clear all

%% read in image
img_name = '../data/peach.png';
img = im2double(imread(img_name));
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
thresh = 0.25;
useGradientShop = false;

%%
% find strong strokes for base layer
[canvas, baseLayer, Gx0, Gy0] = findStrongStrokes(baseLayer, numRows, numCols, wb,...
    thresh, img_grayscale, useGradientShop);
figure;
imshow(canvas);

% find strong strokes for other layers
[canvas, layer1, Gx1, Gy1] = findStrongStrokes(layer1, numRows, numCols, wb/2,...
    thresh, img_grayscale, useGradientShop);
hold on
imshow(canvas);

[canvas, layer2, Gx2, Gy2] = findStrongStrokes(layer2, numRows, numCols, wb/3,...
    thresh, img_grayscale, useGradientShop);
imshow(canvas);

[canvas, layer3, Gx3, Gy3] = findStrongStrokes(layer3, numRows, numCols, wb/6,...
    thresh, img_grayscale, useGradientShop);
imshow(canvas);

%% compute radial basis function
P = [];
P_g = [];
[P,P_g] = getStrongStrokes(P,P_g,baseLayer,Gx0,Gy0,numRows,numCols);
[P,P_g] = getStrongStrokes(P,P_g,layer1,Gx1,Gy1,numRows,numCols);
[P,P_g] = getStrongStrokes(P,P_g,layer2,Gx2,Gy2,numRows,numCols);
[P,P_g] = getStrongStrokes(P,P_g,layer3,Gx3,Gy3,numRows,numCols);

[layer0, temp_x, temp_y] = find_orientation(numRows,numCols,baseLayer,P,P_g);
[layer1, temp_x, temp_y] = find_orientation(numRows,numCols,layer1,P,P_g);
[layer2, temp_x, temp_y] = find_orientation(numRows,numCols,layer2,P,P_g);
[layer3, temp_x, temp_y] = find_orientation(numRows,numCols,layer3,P,P_g);
save('orientation_layers.mat','layer0','layer1','layer2','layer3');

%% visualize strong and all gradients
n = 1;
xs = zeros(numRows,numCols);
ys = zeros(numRows,numCols);
for i = 1:size(P_g,1)
    r = P(i,1);
    c = P(i,2);
    xs(r,c) = P_g(i,1);
    ys(r,c) = P_g(i,2);
end
figure;
quiver(xs(1:n:end,1:n:end), ys(1:n:end,1:n:end),100, '.');
axis image;
axis ij;
hold on;

n = 1;

quiver(temp_x(1:n:end,1:n:end), temp_y(1:n:end,1:n:end),100, '.');
axis image;
axis ij;

%% extract gradients for strong strokes
function [P,T] = getStrongStrokes(P,T,layer,Gx,Gy,numRows,numCols)
for i = 1:size(layer)
    S = layer(i);
    if S.strong
        coords = [S.r S.c];
        if S.r > 5 && S.r < numRows-10 && S.c > 5 && S.c < numCols-5
            P = [P; coords];
            T = [T; Gx(S.r,S.c) Gy(S.r,S.c)];
        end
    end
end
end
%% find strong strokes based on gradient threshold
function [canvas, layer, Gx, Gy] = findStrongStrokes(layer, numRows, numCols, width,...
    thresh, img_grayscale, useGradientShop)

kernelSize = [width width];
kernel = fspecial('gaussian',kernelSize);
img_blur = imfilter(img_grayscale,kernel,'replicate','same');
[Gx,Gy] = imgradientxy(img_blur,'sobel');

if useGradientShop
    gradients = load('long_edge_gradients.mat');
    Gx = gradients.s_x;
    Gy = gradients.s_y;
else
    kernelSize = [width width];
    kernel = fspecial('gaussian',kernelSize);
    img_blur = imfilter(img_grayscale,kernel,'replicate','same');
    [Gx,Gy] = imgradientxy(img_blur,'sobel');
end

numStrong = 0;
for i = 1:size(layer)
    S = layer(i);
    v = [Gx(S.r,S.c), Gy(S.r,S.c)];
    if norm(v) >= thresh
        numStrong = numStrong + 1;
        S.strong = 1;
        S.ang = mod(atan2(v(2),v(1)) + pi/2, 2*pi);
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

disp(sprintf('numStrong: %d\n', numStrong));
end

% estimate gradients based on nearest strong gradients
function [layer, temp_x, temp_y] = find_orientation(numRows,numCols,layer,P,P_g)
temp_x = zeros(numRows, numCols);
temp_y = zeros(numRows, numCols);

nearest_num = 20;
for s=1:size(layer,1)
    curr_stroke = layer(s);
    
    if curr_stroke.strong
        continue
    end
    
    r = curr_stroke.r;
    c = curr_stroke.c;
    inv_dists = zeros(size(P,1),2);
    for i=1:size(P,1)
        P_r = P(i,1);
        P_c = P(i,2);
        dist = sqrt((r-P_r)^2 + (c-P_c)^2);
        inv_dists(i,:) = [1/dist, i];
    end
    sorted_dists = sortrows(inv_dists,1,'descend');
    sorted_dists = sorted_dists(1:nearest_num,:);
    total_dists = sum(sorted_dists(:,1));
    
    x_grad = 0;
    y_grad = 0;
    for i=1:size(sorted_dists,1)
        weight = sorted_dists(i,1) / total_dists;
        x_grad = x_grad + weight * P_g(sorted_dists(i,2),1);
        y_grad = y_grad + weight * P_g(sorted_dists(i,2),2);
    end
    temp_x(r,c) = x_grad;
    temp_y(r,c) = y_grad;
    curr_stroke.ang = mod(atan2(y_grad,x_grad) + pi/2, 2*pi);
    if isnan(curr_stroke.ang)
        disp(sorted_dists);
    end
    layer(s) = curr_stroke;
end
end
