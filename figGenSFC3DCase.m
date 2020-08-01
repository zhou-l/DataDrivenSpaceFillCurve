function [myLT, myVO, HLT, HVO, lineLT, lineVO, dimX, dimY, dimZ] = figGenSFC3DCase(filename)
close all;
global V;
if nargin < 1 
    useSphereVol = true;
elseif isempty(filename)
    useSphereVol = true;
else
    useSphereVol = false;
end

% useSphereVol = true;
dimY = 4;
dimX = 4;
dimZ = 4;
V = zeros(dimY,dimX,dimZ);
V(1:2,1:2,1:2) = 1;

if useSphereVol
    dimX = 16;
    dimY = 16;
    dimZ = 16;
    nSphere = 1;
    V = testVolCreate3D(dimX, dimY, dimZ, nSphere);
    V = V ./ 255;
    V = imdiffusefilt(V);
    V = V .* 255;

%     volumeViewer(V);
%     filepath = [];
%     name = 'testSphere';
%     figure;
%     isoVal = 128;
%     patch(isocaps(V,isoVal),...
%    'FaceColor','interp','EdgeColor','none');
% patch(isocaps(V,isoVal),...
%    'FaceColor','interp','EdgeColor','none');
% p1 = patch(isosurface(V,isoVal),...
%    'FaceColor','blue','EdgeColor','none');
% isonormals(V,p1);
% view(3); 
% axis vis3d tight
% camlight left
% colormap('jet');
% lighting gouraud


else
    % TODO: load from file!
    % filename = './data/11-12_LoResCleaver_30PerEndoCyl_Labeled_Solutions010.nhdr';
    %filename = './data/heartCrop.nhdr';
    headerInfo = nhdr_nrrd_read(filename, true);
    [filepath, name, ext] = fileparts(filename);
    ltfilename = strcat(filepath, '/LT', name, '.csv');
    vofilename = strcat(filepath, '/VO', name, '.csv');
    % headerInfo = nhdr_nrrd_read('./data/heartCrop.nhdr', true);
    V = headerInfo.data;
    dimX = size(V,2);
    dimY = size(V,1);
    dimZ = size(V,3);
    for z = 1:dimZ
        for y = 1:dimY
            for x = 1:dimX
               if isnan(V(x,y,z))
                   V(x,y,z) = -6;
               end
            end
        end
    end
    volumeViewer(V);
end


% slice(V,64,64,64);
%% Ours
% Compute context-based SFC
[cLT, cVisitOrder] = SFCmine3D(V);
myLT = cLT;
myVO = cVisitOrder;
col = 1:length(myLT);
dlmwrite(ltfilename, cLT, ',');
dlmwrite(vofilename, cVisitOrder - [1 1 1], ',');
%% scan line
[LT, visitOrder] = linearizeScanline3D(V);
figure, plot(1:length(LT),LT(:,1));
title('Scanline');
lineLT = LT;
lineVO = visitOrder;
figure;
% p = plot3(visitOrder(:,2), visitOrder(:,1), visitOrder(:,3) ,'b-');
% p.LineWidth = 2;
 lineColorCoded(visitOrder(:,2)',visitOrder(:,1)', visitOrder(:,3)',col);
%% hilbert
order = log2(dimX);
% order = 3;
[hX,hY,hZ] = hilbert3(order);
maxHX = max(hX);
maxHY = max(hY);
maxHZ = max(hZ);
hX = floor((2^order-1) .* 0.5 .* (hX + maxHX) ./ maxHX + 1);
hY = floor((2^order-1) .* 0.5 .* (hY + maxHY) ./ maxHY + 1);
hZ = floor((2^order-1) .* 0.5 .* (hZ + maxHZ) ./ maxHZ + 1);
hX = hX';
hY = hY';
hZ = hZ';
LT = zeros(length(hX),1);
for i = 1:length(hX)
    LT(i,:) = V(hY(i,:),hX(i,:),hZ(i,:));
    if isnan(LT(i,:))
        LT(i,:) = -6;
    end
end
figure, plot(1:length(LT),LT(:,1));


title('Hilbert');
% hX = hX;
% hY = dimY - (hY) + 1;
% hZ = hZ ;
M = int16([hX,hY,hZ]);
figure;
lineColorCoded(hY', hX', hZ', col);

HLT = LT;
HVO = M;

ltfilename = strcat(filepath, '/LTHilbert', name, '.csv');
vofilename = strcat(filepath, '/VOHilbert', name, '.csv');

dlmwrite(vofilename, M - int16([1 1 1]), ',');
dlmwrite(ltfilename,LT, ',');


