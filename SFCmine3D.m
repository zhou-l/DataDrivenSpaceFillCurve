
function [LT, visitOrder] = SFCmine3D(I)
    global V;
    global useLocWeight;
    V = I;
    dimX = size(V,2);
    dimY = size(V,1);
    dimZ = size(V,3);
    
%     figure, montage(I,'DisplayRange',[])
    [WpXR, WpXL, WpYU, WpYD, WpZF, WpZB] = buildSmallCircsDualGraph3D(V);
    % figure, imagesc(WpXR);
    % figure, imagesc(WpXL);
    % figure, imagesc(WpYD);
    % figure, imagesc(WpYU);
     useLocWeight =  true; 
    [Gx, Gy, Gz] = imgradientxyz(V);
    Gmag = sqrt(Gx.*Gx + Gy.*Gy + Gz.*Gz);
    [GpXR, GpXL, GpYU, GpYD, GpZF, GpZB] = buildSmallCircsDualGraph3D(Gmag);
    disp('dual graph built!');

    [T,mstSet] = findMinSpanTree3D2(WpXR, WpXL, WpYU, WpYD, WpZB, WpZF, GpXR, GpXL, GpYU, GpYD, GpZB, GpZF);
    disp('mst tree built.');
    [LT,visitOrder] = linearizeHamCycleMerge3D(T, V, [dimY,1,1]);
    disp('hamCyc built.');
    figure, plot(1:length(LT),LT(:,1));  %plot(1:length(LT),LT(:,2));
    title('Ours');

     figure;
%      p = plot3(visitOrder(:,1),visitOrder(:,2),visitOrder(:,3),'b-');
%  p.LineWidth = 2;
     x = visitOrder(:,2)';
     y = visitOrder(:,1)';
     z = visitOrder(:,3)';
     col = 1:length(visitOrder);
     lineColorCoded(x,y,z,col);

    
    
     
     testHam = true;
     if testHam
         occ = zeros(dimY,dimX,dimZ);
          for i = 1:length(visitOrder)
              occ(visitOrder(i,1),visitOrder(i,2),visitOrder(i,3)) = occ(visitOrder(i,1),visitOrder(i,2),visitOrder(i,3)) + 1;
          end
         badP = find(occ(:,:,:) ~= 1);
         if ~isempty(badP)
             disp(badP);
         else
             disp('all voxels are visited exactly once!');
         end
     end

 
     
    travOrder = zeros(size(V,1), size(V,2), size(V,3));
%     figure;
%     hold on;
    for i = 1:length(visitOrder)
        travOrder(visitOrder(i,1),visitOrder(i,2),visitOrder(i,3)) = i;
    end
%     hold off;
%     figure, montage(travOrder,'DisplayRange',[])
%       figure; volumeViewer(travOrder);
end