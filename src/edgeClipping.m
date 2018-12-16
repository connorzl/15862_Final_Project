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
img_grayscale = rgb2gray(img_large);
wb = 36;

layers = load('orientation_layers.mat');
layer0 = layers.layer0;
layer1 = layers.layer1;
layer2 = layers.layer2;
layer3 = layers.layer3;

layer_masks = load('layer_masks.mat');
mask0 = layer_masks.baseLayerMask;
mask1 = layer_masks.layer1Mask;
mask2 = layer_masks.layer2Mask;
mask3 = layer_masks.layer3Mask;

%% blur layers
sigma_0 = round(0.2*wb);
img_0 = imgaussfilt(img_grayscale, sigma_0);
edge_img_0 = ~(edge(img_0,'Canny',0.4));

[layer0,layerMask0] = add_strokes_to_layer(numRows,numCols,layer0,edge_img_0,wb,true);
[layer1,layerMask1] = add_strokes_to_layer(numRows,numCols,layer1,mask1,wb/2,false);
[layer2,layerMask2] = add_strokes_to_layer(numRows,numCols,layer2,mask2,wb/3,false);
[layer3,layerMask3] = add_strokes_to_layer(numRows,numCols,layer3,mask3,round(wb/4),false);

save('edgeclip_layers.mat','layer0','layerMask0','layerMask1','layerMask2','layerMask3','layer1','layer2','layer3');

function [layer, mask] = add_strokes_to_layer(numRows,numCols,layer,edge_img,wb,is_base_layer)
mask = zeros(numRows, numCols);
inc = 0;
for s = 1:size(layer,1)
    if mod(inc,50) == 0
        disp(inc);
    end
    inc = inc + 1;
    curr_stroke = layer(s);
    ang = curr_stroke.ang;
    
    dX = 1;
    dY = tan(ang);
    if dY > tan(pi/2-0.01)
        dX = 0;
        dY = 1;
    end
    
    [mask,stroke_length1,mask_pixels1] = grow_stroke(dX,dY,wb,numRows,numCols,edge_img,curr_stroke,mask,is_base_layer);
    curr_stroke.l1 = stroke_length1;
    
    dX = -dX;
    dY = -dY;
    [mask,stroke_length2,mask_pixels2] = grow_stroke(dX,dY,wb,numRows,numCols,edge_img,curr_stroke,mask,is_base_layer);
    curr_stroke.l2 = stroke_length2;
    
    curr_stroke.stroke_pixels = [mask_pixels1; mask_pixels2];
    layer(s) = curr_stroke;
end
end

function [mask,stroke_length,mask_pixels] = grow_stroke(dX,dY,wb,numRows,numCols,layer_mask,curr_stroke,mask,is_base_layer)
no_edge_hit = true;
stroke_length = 0;
y = curr_stroke.r;
x = curr_stroke.c;
local_mask = zeros(numRows,numCols);
while no_edge_hit
            
    circle = [];
    for new_y=round(y-wb/2):round(y+wb/2)
        for new_x=round(x-wb/2):round(x+wb/2)
            if new_y < 1 || new_y > numRows || new_x < 1 || new_x > numCols
                continue
            end
            if sqrt((new_x-x)^2 + (new_y-y)^2) <= wb/2
                circle = [circle; new_y new_x];
            end
        end
    end
    
    if no_edge_hit
        if sqrt((x-curr_stroke.c)^2 + (y-curr_stroke.r)^2) > 2 * wb
            break
        else
            stroke_length = sqrt((x-curr_stroke.c)^2 + (y-curr_stroke.r)^2);
        end
        for i=1:size(circle,1)
            local_mask(circle(i,1), circle(i,2)) = 1;
            mask(circle(i,1), circle(i,2)) = 1;
        end
    end
    
    if stroke_length > 2 * wb
        break
    end
    
    x = x + dX;
    y = y + dY;
    if (x < 1 || x > numCols || y < 1 || y > numRows)
        break
    end

    if ~layer_mask(round(y), round(x))
        no_edge_hit = false;
        break
    end
end

mask_idxs = find(local_mask);
mask_pixels = zeros(size(mask_idxs,1),2);
for j=1:size(mask_idxs,1)
    idx = mask_idxs(j);
    idx_c = ceil(idx / numRows);
    idx_r = mod(idx-1, numRows) + 1;
    mask_pixels(j,:) = [idx_r,idx_c];
end
end
