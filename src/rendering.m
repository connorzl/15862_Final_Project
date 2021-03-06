%%
close all

%% read in image
img_name = '../data/peach.png';
img = im2double(imread(img_name));
[imh, imw, ~] = size(img);

wb = 36;

canvasScale = 2;
numRows = imh * canvasScale;
numCols = imw * canvasScale;

canvas = ones(numRows,numCols,3);
canvas_alphas = ones(numRows,numCols);

layers = load('color_layers.mat');
layer0 = layers.layer0;
layer1 = layers.layer1;
layer2 = layers.layer2;
layer3 = layers.layer3;

%%
% load texture, pad to be a square, get color for this texture
text_img = im2double(rgb2gray(imread('../data/imp_brushstrokes.jpg')));
textures = text_img(49:146,4:302);
textures = textures(18:79,13:285);
alphas = text_img(49:146,305:603);
alphas = alphas(18:79,13:285);

useAlpha = true;

[canvas0, canvas_alphas0] = renderLayer(layer0,wb,textures,alphas,canvas,canvas_alphas,numRows,numCols,useAlpha);
disp("done canvas0");
[canvas1, canvas_alphas1] = renderLayer(layer1,wb/2,textures,alphas,canvas0,canvas_alphas0,numRows,numCols,useAlpha);
disp("done canvas1");
[canvas2, canvas_alphas2] = renderLayer(layer2,wb/3,textures,alphas,canvas1,canvas_alphas1,numRows,numCols,useAlpha);
disp("done canvas2");
[canvas3, canvas_alphas3] = renderLayer(layer3,wb/6,textures,alphas,canvas2,canvas_alphas2,numRows,numCols,useAlpha);
disp("done canvas3");

function [canvas,canvas_alphas] = ...
    renderLayer(layer,wb,textures,alphas,canvas,canvas_alphas,numRows,numCols,useAlpha)
for s=1:size(layer,1)
    if mod(s,50) == 0
        disp(s);
    end
    
    stroke = layer(s);
    if stroke.l1 + stroke.l2 == 0
        continue
    end
    
    color_texture = get_color_texture(textures, stroke.color);
    % Scale texture and alpha
    stroke_texture = imresize(color_texture,[wb, stroke.l1 + stroke.l2]);
    stroke_alpha = imresize(alphas, [wb, stroke.l1 + stroke.l2]);
    
    % center of the stroke
    row = round(size(stroke_texture,1)/2);
    if (stroke.ang > pi/2 && stroke.ang < 3*pi/2)
        col = round(stroke.l1 / (stroke.l1 + stroke.l2) * size(stroke_texture,2));
    else
        col = round(stroke.l2 / (stroke.l1 + stroke.l2) * size(stroke_texture,2));
    end
    
    % pad texture and alpha
    centered_texture = center_im(stroke_texture,col);
    centered_alpha = center_im(stroke_alpha,col);
    
    % rotate texture and alpha
    angle = stroke.ang * 360 / (2 * pi);
    angle_texture = imrotate(centered_texture, -angle);
    angle_alpha = imrotate(centered_alpha, -angle);
    
    % translate to canvas
    row_start = floor(stroke.r - size(angle_texture,1)/2);
    col_start = floor(stroke.c - size(angle_texture,2)/2);
    
    for i = row_start:row_start+size(angle_texture,1)-1
        for j = col_start:col_start+size(angle_texture,2)-1
            if i < 1 || i > numRows || j < 1 || j > numCols
                continue
            end
            
            if angle_alpha(i-row_start+1,j-col_start+1) > 0
                texture_pixel = angle_texture(i-row_start+1,j-col_start+1,:);
                
                if useAlpha
                    texture_alpha = angle_alpha(i-row_start+1,j-col_start+1);
                    A_prime = squeeze(canvas_alphas(i,j) * canvas(i,j,:));
                    B_prime = squeeze(texture_alpha * texture_pixel);
                    
                    canvas(i,j,:) = B_prime + (1 - texture_alpha) * A_prime;
                    canvas_alphas(i,j) = texture_alpha + (1-texture_alpha) * canvas_alphas(i,j);
                else
                    canvas(i,j,:) = texture_pixel;
                end
            end
        end
    end
end
end

function centered_im = center_im(im,col)
[imh,imw,~] = size(im);
centered_im = im;
left_pad = 0;
top_pad = 0;
if imh < imw
    rad = max(col,imw-col);
    if col > imw-col
        top_pad = rad-imh/2;
        bot_pad = rad-(imh/2);
        right_pad = rad-(imw-col);
        centered_im = padarray(centered_im,[0;right_pad],0,'post');
    else
        top_pad = rad-imh/2;
        bot_pad = rad-imh/2;
        left_pad = rad-col;
        centered_im = padarray(centered_im,[0;left_pad],0,'pre');
    end
    centered_im = padarray(centered_im,top_pad,'pre');
    centered_im = padarray(centered_im,bot_pad,'post');
else
    rad = imh/2;
    left_pad = rad-col;
    right_pad = rad-(imw-col);
    if left_pad > 0
        centered_im = padarray(centered_im,[0;left_pad],0,'pre');
    end
    
    if right_pad > 0
        centered_im = padarray(centered_im,[0;right_pad],0,'post');
    end
end
end

function color_texture = get_color_texture(textures, color)
color_texture = zeros(size(textures,1),size(textures,2),3);
color_texture(:,:,1) = textures * color(1);
color_texture(:,:,2) = textures * color(2);
color_texture(:,:,3) = textures * color(3);
end