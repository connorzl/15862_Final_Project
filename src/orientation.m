close all

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
thresh = 0.15;

% find strong strokes for base layer
[canvas, baseLayer, Gx0, Gy0] = findStrongStrokes(baseLayer, numRows, numCols, wb,...
    thresh, img_grayscale);
figure;
imshow(canvas);

% find strong strokes for other layers
[canvas, layer1, Gx1, Gy1] = findStrongStrokes(layer1, numRows, numCols, wb/2,...
    thresh, img_grayscale);
hold on
imshow(canvas);

[canvas, layer2, Gx2, Gy2] = findStrongStrokes(layer2, numRows, numCols, wb/4,...
    thresh, img_grayscale);
imshow(canvas);

[canvas, layer3, Gx3, Gy3] = findStrongStrokes(layer3, numRows, numCols, round(wb/8),...
    thresh, img_grayscale);
imshow(canvas);

%% compute radial basis function
P = [];
Tx = [];
[P,Tx] = getStrongStrokes(P,Tx,baseLayer,Gx0,numRows,numCols);
[P,Tx] = getStrongStrokes(P,Tx,layer1,Gx1,numRows,numCols);
[P,Tx] = getStrongStrokes(P,Tx,layer2,Gx2,numRows,numCols);
[P,Tx] = getStrongStrokes(P,Tx,layer3,Gx3,numRows,numCols);

eg = 0.001;    % sum-squared error goal
sc = 20;    % spread constant
net_x = newrb(P,Tx,eg,sc);

P = [];
Ty = [];
[P,Ty] = getStrongStrokes(P,Ty,baseLayer,Gy0,numRows,numCols);
[P,Ty] = getStrongStrokes(P,Ty,layer1,Gy1,numRows,numCols);
[P,Ty] = getStrongStrokes(P,Ty,layer2,Gy2,numRows,numCols);
[P,Ty] = getStrongStrokes(P,Ty,layer3,Gy3,numRows,numCols);

net_y = newrb(P,Ty,eg,sc);

%% visualize strong gradients
n = 1;
xs = zeros(numRows,numCols);
ys = zeros(numRows,numCols);
for i = 1:size(Tx,2)
    y = P(1,i);
    x = P(2,i);
    xs(y,x) = Tx(i);
    ys(y,x) = Ty(i);
end
quiver(xs(1:n:end,1:n:end), ys(1:n:end,1:n:end),100, '.');
axis image;
axis ij;
%%
input = zeros(2,numRows*numCols);
for i = 1:numRows
    for j = 1:numCols
        coord = [i; j];
        input(:,(i-1)*numCols + j) = coord;
    end
end
disp('done building input');
res_x = net_x(input);
res_y = net_y(input);
disp('done evaluating');

xs = zeros(numRows,numCols);
ys = zeros(numRows,numCols);
for i = 1:numRows*numCols
    r = ceil(i / numCols);
    c = mod(i-1, numCols) + 1;
    xs(r,c) = res_x(i);
    ys(r,c) = res_y(i);
end

%% visualize
figure;
n = 3;

quiver(xs(1:n:end,1:n:end), ys(1:n:end,1:n:end),10, '.');
axis image;
axis ij;

%% extract gradients for strong strokes

function [P,T] = getStrongStrokes(P,T,layer,G,numRows,numCols)
for i = 1:size(layer)
    S = layer(i);
    if S.strong
        coords = [S.r; S.c];
        if S.r > 5 && S.r < numRows-10 && S.c > 5 && S.c < numCols-5
            P = [P coords];
            T = [T G(S.r,S.c)];
        end
    end
end
end
%% find strong strokes
function [canvas, layer, Gx, Gy] = findStrongStrokes(layer, numRows, numCols, width,...
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
