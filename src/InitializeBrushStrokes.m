%% read in image
img_name = '../data/peach.png';
img = im2double(imread(img_name));
[imh, imw, ~] = size(img);

canvasScale = 2;
numRows = imh * canvasScale;
numCols = imw * canvasScale;

% regeneration width
wr = 18;

% brush width
wb = 36;

% default brush stroke
N = ceil(numCols / wr * numRows / wr);
S = struct('r',0,'c',0,'ang',0,'w',0,'l1',0,'l2',0,...
    'color',[0,0,0],'opacity',0,'strong',0,'stroke_pixels',[]);

%% compute base layer
baseLayer = repmat(S, N, 1);

baseLayerMask = zeros(numRows, numCols);
baseLayerMask(wr+1:numRows-wr-1, wr+1:numCols-wr-1) = 1;

baseLayer = computeBrushStrokes(baseLayerMask,baseLayer,numRows,numCols,wr);
visualizeLayer(baseLayer,numRows,numCols,baseLayerMask);

%% initialize layer parameters
img_large = imresize(img, canvasScale);
img_grayscale = rgb2gray(img_large);
sigma_1 = 0.3 * wb;
sigma_2 = 0.2 * wb;
sigma_3 = 0.15 * wb;

%% initialize layer 1
layer1 = repmat(S, N, 1);
img_1 = edge(img_grayscale,'Canny',[],sigma_1);
layer1Mask = computeLayerMask(wr, find(img_1), numRows, numCols);

layer1 = computeBrushStrokes(layer1Mask,layer1,numRows,numCols,wr);
visualizeLayer(layer1,numRows,numCols,layer1Mask);

%% initialize layer 2
layer2 = repmat(S, N, 1);
img_2 = edge(img_grayscale,'Canny',[],sigma_2);
layer2Mask = computeLayerMask(0.5*wr, find(img_2), numRows, numCols);

layer2 = computeBrushStrokes(layer2Mask,layer2,numRows,numCols,wr);
visualizeLayer(layer2,numRows,numCols,layer2Mask);

%% initialize layer 3
layer3 = repmat(S, N, 1);
img_3 = edge(img_grayscale,'Canny',[],sigma_3);
layer3Mask = computeLayerMask(round(0.25*wr), find(img_3), numRows, numCols);

layer3 = computeBrushStrokes(layer3Mask,layer3,numRows,numCols,wr);
visualizeLayer(layer3,numRows,numCols,layer3Mask);

%% save to file
save('layers.mat','baseLayer','layer1','layer2','layer3');
save('layer_masks.mat','baseLayerMask','layer1Mask','layer2Mask','layer3Mask');

%% helper functions
function layerMask = computeLayerMask(rad, inds, numRows, numCols)
layerMask = zeros(numRows, numCols);
for i = 1:size(inds,1)
    index = inds(i);
    c = ceil(index / numRows);
    r = mod(index-1, numRows) + 1;
    
    for x = c-rad:c+rad
        for y = r-rad:r+rad
            if x < 1 || x > numCols || y < 1 || y > numRows
                continue
            end
            
            v = [x-c, y-r];
            if norm(v) <= rad
                layerMask(y,x) = 1;
            end
        end
    end
end
end

function [strokes] = computeBrushStrokes(layerMask,strokes,numRows,numCols,wr)
numStrokes = 0;
while any(layerMask(:))
    validCenters = find(layerMask);
    index = validCenters(randi(size(validCenters,1),1));
    
    c = ceil(index / numRows);
    r = mod(index-1, numRows) + 1;
    
    numStrokes = numStrokes + 1;
    S = strokes(numStrokes);
    S.r = r;
    S.c = c;
    strokes(numStrokes) = S;
    
    % remove possible centers
    for i = r-wr:r+wr
        for j = c-wr:c+wr
            if i < 1 || i > numRows || j < 1 || j > numCols
                continue
            end
            v = [r-i, c-j];
            if norm(v) <= wr
                layerMask(i,j) = 0;
            end
        end
    end
end

strokes = strokes(1:numStrokes);
end

% visualize results
function [] = visualizeLayer(layer,numRows,numCols,layerMask)
baseCanvas = zeros(numRows, numCols);
for i=1:size(layer)
    S = layer(i);
    baseCanvas(S.r, S.c) = 1;
end
figure;
imshow(baseCanvas);
figure;
imshow(layerMask);
end