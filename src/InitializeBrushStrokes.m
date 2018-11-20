%% read in image
img = imread('../data/peach.png');
[imh, imw, ~] = size(img);

canvasScale = 2;
numLayers = 4;

numRows = imh * canvasScale;
numCols = imw * canvasScale;

% regeneration width
wr = 15;

N = ceil(numCols / wr * numRows / wr);
S = struct('r',0,'c',0,'ang',0,'w',0,'l1',0,'l2',0,...
    'color',[0,0,0],'opacity',0);
baseLayer = repmat(S, N, 1);

possibleCenters = zeros(numRows, numCols);
possibleCenters(wr:numRows-wr, wr:numCols-wr) = 1;

layer1 = [];
layer2 = [];
layer3 = [];

%% compute strokes for base layer
numStrokes = 0;
while existsGap(baseLayer, numStrokes, wr)
    index = randi(numRows * numCols, 1);
    r = ceil(index / numCols);
    c = mod(index-1, numCols) + 1;
    
    validGap = false;
    for i = 1:numStrokes
        S = baseLayer(i);
        
        v = [S.r - r, S.c - c];
        if norm(v) >= 2*wr
            validGap = true;
            break
        end
    end
    
    if validGap || numStrokes == 0
        numStrokes = numStrokes + 1;
        S = baseLayer(numStrokes);
        S.r = r;
        S.c = c;
        baseLayer(numStrokes) = S;
        
        % 
    end
end

%baseLayer = baseLayer(1:numStrokes);

function validGap = existsGap(layer, numStrokes, wr)
validGap = false;
for i = 1:numStrokes
    S1 = layer(i);
    for j = i+1:numStrokes
        S2 = layer(j);
        v = [S1.r - S2.r, S1.c - S2.c];
        if norm(v) >= 2*wr
            validGap = true;
            break
        end
    end
end
disp(validGap);
end

