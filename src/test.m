%% 11, 18

for r=1:size(textures,1)
    hit = false;
    for c=1:size(textures,2)
        if textures(r,c) < 252/255
            hit = true;
            disp(r);
            return;
        end
    end
end

%% 285
for r=size(textures,1):-1:1
    hit = false;
    for c=1:size(textures,2)
        if textures(r,c) < 252/255
            hit = true;
            disp(r);
            return;
        end
    end
end

%% 17
for c=1:size(alphas,2)
    for r=1:size(alphas,1)
        if alphas(r,c) > 0
            hit = true;
            disp(c);
            return;
        end
    end
end

%% 288
for c=size(alphas,2):-1:1
    for r=1:size(alphas,1)
        if alphas(r,c) > 0
            hit = true;
            disp(c);
            return;
        end
    end
end
