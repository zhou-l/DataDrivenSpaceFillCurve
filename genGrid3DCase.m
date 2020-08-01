
data = zeros(4*4*4,3);
C = zeros(4*4*4,3);
cnt = 1;
for z = 1:4
    for y = 1:4
        for x = 1:4
            data(cnt,:) =[x,y,z];
            if mod(z,2) == 0
                C(cnt,:) = [0.5,0.5,0.5];
            else
                C(cnt,:) = [0,0,1];
            end
            cnt = cnt + 1;
        end
    end
end
S = 100;
figure, scatter3(data(:,1),data(:,2),data(:,3), S, C, 'filled');