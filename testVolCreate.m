function V = testVolCreate(dimX, dimY, dimZ, nSphere, R)
if nargin < 3
    dimZ = 1;
    nSphere = 1;
end

if nargin < 4
    nSphere = 1;
end

if nargin < 5
    R = 4;
end

V = zeros(dimX, dimY, dimZ);

minR = 3;
maxR = R;
% R = 10;
Rlist = zeros(nSphere, 1);
for i = 1:length(Rlist)
    Rlist(i) = round(minR + rand()*(maxR - minR));
end
Ctrs = zeros(nSphere, 3);


% if dimZ == 1
%     Ctrs(1,:) = [10,10,1];
    rctrs = rand(nSphere,3);
    Ctrs(:,:,:) = rctrs(:,:,:) .* [dimX, dimY, dimZ];
    if dimZ == 1
        Ctrs(:,3) = 1;
    end
%     if nSphere > 1
%         Ctrs(2,:) = [20,25,1];
%     end
%     
%     if nSphere > 2
%         Ctrs(3,:) = [25,12,1];
%     end
%     
%     if nSphere > 3
%         Ctrs(4,:) = [7,27,1];
%     end
%     
%     if nSphere > 4
%         Ctrs(5,:) = [25, 20, 1];
%     end
% end

for z = 1:dimZ
    for y = 1:dimY
        for x = 1:dimX
            for i = 1:nSphere
                dist = norm([y,x,z] - Ctrs(i,:));
%                if dist <= R
                if dist <= Rlist(i)
                   if i == 2
                        V(y,x,z) = max(-255/R * dist+ 255, 0);
                   else
                   V(y,x,z) = max(-255/R * dist+ 255,0);
                   end
               end
            end
        end
    end
end
% imagesc(V);
% colorbar;
% fn = sprintf('testSph%i_%i_%i_%i.txt', nSphere, dimX, dimY, dimZ);
% dlmwrite(fn, V);
 