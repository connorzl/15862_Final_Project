
layers = load('color_layers.mat');
layer0 = layers.layer0;
layer1 = layers.layer1;
layer2 = layers.layer2;
layer3 = layers.layer3;

wb = 36;

% load texture, pad to be a square, get color for this texture
text_img = im2double(rgb2gray(imread('./data/imp_brushstrokes.jpg')));
textures = text_img(49:146,4:302);
textures = textures(18:79,13:285);
alphas = text_img(49:146,305:603);
alphas = alphas(18:79,13:285);

for s=1:size(layer0,1)
    stroke = layer0(s)
    
    if stroke.l1 + stroke.l2 == 0
        continue
    end
    
    color_texture = get_color_texture(textures, stroke.color);
    % Scale texture and alpha
    stroke_texture = imresize(color_texture,[wb, stroke.l1 + stroke.l2]);
    stroke_alpha = imresize(alphas, [wb, stroke.l1 + stroke.l2]);
    row = round(size(stroke_texture,1)/2);
    col = round(stroke.l2 / (stroke.l1 + stroke.l2) * size(stroke_texture,2));
    centered_texture = center_im(stroke_texture,row,col);
    centered_alpha = center_im(stroke_alpha,row,col);
    angle_texture = imrotate(centered_texture, stroke.ang);
    angle_alpha = imrotate(centered_alpha, stroke.ang);
end

% scale texture to have length length1 + length2

% find texture origin based on ratio of length1 and length2

% rotate texture by angle

% translate texture to the brush stroke center and alpha composite

function centered_im = center_im(im,row,col)
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
    centered_im = padarray(centered_im,[0;left_pad],0,'pre');
    centered_im = padarray(centered_im,[0;right_pad],0,'post');
end

figure;
imshow(centered_im);
hold on;
plot(col+left_pad,row+top_pad,'b*');
pause;

    
end


function color_texture = get_color_texture(textures, color)
color_texture = zeros(size(textures,1),size(textures,2),3);
color_texture(:,:,1) = textures * color(1);
color_texture(:,:,2) = textures * color(2);
color_texture(:,:,3) = textures * color(3);
end