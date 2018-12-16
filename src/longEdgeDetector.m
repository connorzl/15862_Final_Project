close all
clear all

%%
img_name = '../data/peach.png';

img = im2double(imread(img_name));
[imh, imw, ~] = size(img);

canvasScale = 2;
numRows = imh * canvasScale;
numCols = imw * canvasScale;

img_large = imresize(img, canvasScale);
img_grayscale = rgb2gray(img_large);

[Gmag,Gdir] = imgradient(img_grayscale);
Gdir = Gdir + 360;
Gdir = mod(Gdir,360);

% compute p_m
p_m = zeros(numRows,numCols);
mu = zeros(numRows,numCols);
sigma = zeros(numRows,numCols);

% compute mu
for i = 1:numRows
    for j = 1:numCols
        validPixels = 0;
        sum = 0;
        for l = i-2:i+2
            for k = j-2:j+2
                if l < 1 || l > numRows || k < 1 || k > numCols
                    continue
                end
                validPixels = validPixels + 1;
                sum = sum + Gmag(l,k);
            end
        end
        mu(i,j) = sum / validPixels;
    end
end

% compute sigma
for i = 1:numRows
    for j = 1:numCols
        validPixels = 0;
        sqDiffs = 0;
        for l = i-2:i+2
            for k = j-2:j+2
                if l < 1 || l > numRows || k < 1 || k > numCols
                    continue
                end
                validPixels = validPixels + 1;
                sqDiffs = sqDiffs + (Gmag(l,k) - mu(i,j))^2;
            end
        end
        
        sigma(i,j) = sqrt(sqDiffs / (validPixels-1));
    end
end

% compute pm
for i = 1:numRows
    for j = 1:numCols
        p_m(i,j) = ( Gmag(i,j) - mu(i,j) ) /  (sigma(i,j) + 0.000001);
    end
end

%% compute m0 and m1
numIterations = 60;
m0 = zeros(numRows,numCols,numIterations);
m1 = zeros(numRows,numCols,numIterations);
var = pi / 5;
disp("Computing m0 and m1");
for t = 2:numIterations
    for i = 1:numRows
        for j = 1:numCols
            % m0
            p_dir = Gdir(i,j);   
            x = sqrt(2) * cos(p_dir * 180 / pi);
            y = sqrt(2) * sin(p_dir * 180 / pi);
            [q_xs,q_ys,weights] = bilinearInterpolation(x,y);
            for k = 1:4
                q_x = q_xs(k);
                q_y = q_ys(k);
                if q_x < 1 || q_x > numCols || q_y < 1 || q_y > numRows
                    continue
                end
                q_dir = Gdir(q_y,q_x);
                
                w_alpha = weights(k);
                w_theta = exp( - (p_dir - q_dir)^2 / (2 * var^2));
                q_m = p_m(q_y,q_x);
                
                m0(i,j,t) = w_alpha * w_theta * (q_m + m0(q_y,q_x,t-1));
            end
            
            % m1, only difference is projection angle
            p_dir = mod(Gdir(i,j) + 180,360);    
            x = sqrt(2) * cos(p_dir * 180 / pi);
            y = sqrt(2) * sin(p_dir * 180 / pi);
            [q_xs,q_ys,weights] = bilinearInterpolation(x,y);
            for k = 1:4
                q_x = q_xs(k);
                q_y = q_ys(k);
                if q_x < 1 || q_x > numCols || q_y < 1 || q_y > numRows
                    continue
                end
                q_dir = Gdir(q_y,q_x);
                
                w_alpha = weights(k);
                w_theta = exp( - (p_dir - q_dir)^2 / (2 * var^2));
                q_m = p_m(q_y,q_x);
                
                m1(i,j,t) = w_alpha * w_theta * (q_m + m1(q_y,q_x,t-1));
            end  
        end
    end
    disp(t);
end

%% now compute local gradient saliency
e_l = m0(:,:,60) + m1(:,:,60) + p_m;
e_o = Gdir;

[u_x,u_y] = imgradientxy(img_grayscale);
s_x = cos(e_o).^2 .* e_l .* u_x;
s_y = sin(e_o).^2 .* e_l .* u_y;

figure;
n = 2;
quiver(s_x(1:n:end,1:n:end), s_y(1:n:end,1:n:end),100, '.');
axis image;
axis ij;

figure;
n = 2;
quiver(u_x(1:n:end,1:n:end), u_y(1:n:end,1:n:end),100, '.');
axis image;
axis ij;

save('long_edge_gradients.mat', 's_x', 's_y');

%% bilinear interpolation function
function [q_x, q_y, weights] = bilinearInterpolation(x,y)

x_nw = floor(x);
y_nw = floor(y);

x_ne = ceil(x);
y_ne = floor(y);

x_sw = floor(x);
y_sw = ceil(y);

x_se = ceil(x);
y_se = ceil(y);

q_x = [x_nw x_ne x_sw x_se];
q_y = [y_nw y_ne y_sw y_se];
q = [q_x;q_y];

low_x = 1 - (x - floor(x));
high_x = 1 - low_x;

low_y = 1 - (y - floor(y));
high_y = 1 - low_y;

weights = [low_x*low_y high_x*low_y low_x*high_y high_x*high_y];
end