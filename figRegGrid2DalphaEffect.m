%% Test the alpha parameter 
% return values:
% listACvals -- has numAlphas rows of autocorrelations of data values
% listACdist -- has numAlphas rows of autocorrelations of Eucledian
% distance
function [listACvals, listACdist]= figRegGrid2DalphaEffect(filename, alphaList)
close all;
V = imread(filename);

if length(size(V))==3
    V = double(rgb2gray(V));
else
    V = double(V(:,:,1));
end

if nargin < 2
    alphas = [0, 0.01, 0.2, 0.4, 0.6, 0.8];
else
    alphas = alphaList;
end

  V = V ./ 255;
  
 V = padImgToPow2(V);
 dV = size(V);
 dimY = dV(1);
 dimX = dV(2);
global useLocWeight;
global zalpha;
useLocWeight = true;
maxLags = 64;
figure;
[Gmag,Gdir] = imgradient(V);
[WpXR, WpXL, WpYU, WpYD] = buildSmallCircsDualGraph(V);
[GpXR, GpXL, GpYU, GpYD] = buildSmallCircsDualGraph(Gmag);
numAlphas = 6;
listACvals = zeros(numAlphas, maxLags);
listACdist = zeros(numAlphas, maxLags);
% alphas = [0, 0.001, 0.01, 0.05, 0.1, 0.2];
% alphas = [0, 0.01, 0.2, 0.4, 0.6, 0.8];
numAlphas = length(alphas);
for i = 1:numAlphas
    zalpha = alphas(i);%(i-1)*0.2;
    [T,mstSet] = findMinSpanTree2(WpXR, WpXL, WpYU, WpYD, GpXR, GpXL, GpYU, GpYD);
    [LT, visitOrder] = linearizeHamCycleMerge(T, V, [dimY, 1]);
%     myLT = LT;
%     myVO = visitOrder;
% 
     subplot(3,6, i + 6);
     plot(1:length(LT),LT(:,1));
     if i == 1
         h2 = text(10, 0,'Linearization');
        set(h2, 'rotation', 90);
     end
     
    if i == 1
        h3 = text(0, 0,'Autocorr of Dist');
        set(h3, 'rotation', 90);
    end

%     %     plot(1:length(LT),LT(:,2));
%     hold off;
%     title('Locality Hamilton');

    hamLT = LT;
    hamOrder = visitOrder;
    avgAutoCorrOurs = compAvgAutoCorr(LT, maxLags);
    listACvals(i,:) = avgAutoCorrOurs;
    
    dist = vecnorm(visitOrder(:,:), 2, 2);
    avgAutoCorrDistOurs = compAvgAutoCorr(dist, maxLags);
    listACdist(i,:) = avgAutoCorrDistOurs;
    
    subplot(3,6,i+12);  
    plot(1:length(avgAutoCorrDistOurs), avgAutoCorrDistOurs);
    ylim([0.97, 1.0]);


    % disp(visitOrder)
%     hold on;
    travOrder = zeros(dimY,dimX);
    for j = 1:length(hamOrder)
        travOrder(hamOrder(j,1),hamOrder(j,2)) = j;
    end
    subplot(3, 6, i), imagesc(travOrder);
    strTitle = sprintf('\\alpha = %.3f', zalpha);
    title(strTitle);
    line(hamOrder(:,2), hamOrder(:,1));
    if i == 1
        h1 = text(-1, 0,'SFC');
        set(h1, 'rotation', 90);
    end

%     hold off;
end

figure, imshow(V);