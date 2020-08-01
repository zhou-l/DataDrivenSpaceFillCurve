
function [LT, VO, travOrder] = SFCmine(I)
    global V;
    global useLocWeight;
    V = I;
    figure, imagesc(V);
    dimX = size(V,2);
    dimY = size(V,1);
    [WpXR, WpXL, WpYU, WpYD] = buildSmallCircsDualGraph(V);
    % figure, imagesc(WpXR);
    % figure, imagesc(WpXL);
    % figure, imagesc(WpYD);
    % figure, imagesc(WpYU);
     useLocWeight = true; 
    [Gmag,Gdir] = imgradient(V);
    [GpXR, GpXL, GpYU, GpYD] = buildSmallCircsDualGraph(Gmag);
    [T,mstSet] = findMinSpanTree2(WpXR, WpXL, WpYU, WpYD, GpXR, GpXL, GpYU, GpYD);
%     [LT,visitOrder] = linearizeHamCycle(T, V, [dimY,1]);
    [LT,VO] = linearizeHamCycleMerge(T, V, [dimY,1]);
    travOrder = zeros(size(V,1), size(V,2), 1);
    figure;
    hold on;
    for i = 1:length(VO)
        travOrder(VO(i,1),VO(i,2)) = i;
    end
    hold off;
    figure, imagesc(travOrder);
end