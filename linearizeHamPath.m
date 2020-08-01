%Traversing the image I with a Hamilton path with given entry point and
%exit point without a precomputed min spanning tree

function [LT, visitOrder, exitPix] = linearizeHamPath(I, entryPix, exitEdge, isPlot)
      % simple graph test cases work up to 3*3
%       I = [1 1 1 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
%       entryPix = [3,1];
%       exitPix = [1,3];
%       [Status, visitOrder, LT] = hamPathGrid(I, entryPix, exitPix);
%       return;
%     I = zeros(8);
if nargin < 4
    zisPlot = false;
else
%% for fig genneration
%     I = testVolCreate(8, 8, 1, 1);
    zisPlot = isPlot;
end
% force plot
% zisPlot = true;
if zisPlot
 figure;imagesc(I);
end
% Use nearest neighbor strategy
     if size(I,1) * size(I,2) >= 36
         [Status, LT, visitOrder, exitPix] = hamPathPartitionGrid(I, entryPix, exitEdge);
         if Status == 0
             disp('HamPath Status Error!');
         end
         if zisPlot
            figure;imagesc(I);line(visitOrder(:,2), visitOrder(:,1), 'lineWidth', 2, 'color', 'white');
         end
         return;
     end
     
     exitPixs = findCompatibleExitPixs(size(I), entryPix, exitEdge);
     if  ~isempty(exitPixs)
%          for i = 1:size(exitPixs,1)
%              i = size(exitPixs,1);
%              exitPix = exitPixs(i,:);
%               [Status, visitOrder, LT] = hamPathGrid(I, entryPix, exitPix);
            [Status, visitOrder, LT, exitPix] = hamPathGridRepNN(I, entryPix, exitPixs); 
%             visitOrder
         if Status == 0
             disp('HamPath Status Error!');
         end
%          end
     else
         Status = 0;
         visitOrder = [];
         LT = [];
         exitPix = [];
     end
      if zisPlot
         figure;
         imagesc(I); line(visitOrder(:,2), visitOrder(:,1));
     end
end




