function V = testVolCreate3D(dimX, dimY, dimZ, nSphere)
if nargin < 3
    dimZ = 1;
    nSphere = 1;
end

if nargin < 4
    nSphere = 1;
end

V = zeros(dimX, dimY, dimZ);

R = 5;
Ctrs = zeros(nSphere, 3);


Ctrs(1,:) = [8,8,8];
if nSphere > 1
    Ctrs(2,:) = [20,25,8];
end

if nSphere > 2
    Ctrs(3,:) = [25,12,15];
end

if nSphere > 3
    Ctrs(4,:) = [7,27,10];
end

if nSphere > 4
    Ctrs(5,:) = [25, 20, 17];
end


for z = 1:dimZ
    for y = 1:dimY
        for x = 1:dimX
            for i = 1:nSphere
                dist = norm([y,x,z] - Ctrs(i,:));
               if dist <= R 
                   if i == 2
                        V(y,x,z) = max(-255/R * dist+ 255,0);
                   else
                   V(y,x,z) = max(-255/(R) * dist+ 255,0);
%                    if -255/(R) * dist+ 255 < 100 
%                        V(y,x,z) = 0;
%                    end
                   end
               end
            end
        end
    end
end

% fn = sprintf('testSph%i_%i_%i_%i.txt', nSphere, dimX, dimY, dimZ);
% dlmwrite(fn, V);
 