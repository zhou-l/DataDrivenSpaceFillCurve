%Traversing the image I with a Hamilton path with given entry point and
%exit point without a precomputed min spanning tree

function [LT, visitOrder, exitPix] = linearizeHamPath3D(V, entryPix, exitFace)
      % simple graph test cases work up to 3*3
%       I = [1 1 1 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
%       entryPix = [3,1];
%       exitPix = [1,3];
%       [Status, visitOrder, LT] = hamPathGrid(I, entryPix, exitPix);
%       return;
% Use nearest neighbor strategy

LT = [];
visitOrder = [];
exitPix = [];
global w;
global h;
global d;
     w = size(V,1);
     h = size(V,3);
     d = size(V,2);
     visitOrder = [];
     LT = [];
             exitPix = [];
     if w * h * d > 64
         disp('Cannot handle this case yet.');
         return;
     end
     
     needPartition = false;
     exitPixs = findCompatibleExitPixs3D(size(V), entryPix, exitFace);
     if  ~isempty(exitPixs)
         %% convert the volume connectivity to adjacency matrix
            Madj = createAdjMatrixVol(V, size(V));
            entryId = ((entryPix(3)-1)*h + entryPix(2)-1)*w + entryPix(1)-1+1;
            minCost = realmax;
            minpath_id = -1;
            hamPath = [];
            for i = size(exitPixs,1):-1:1
                if exitPixs(i,1) == entryPix(1) && exitPixs(i,2) == entryPix(2) && exitPixs(i,3) == entryPix(3)
                    continue;
                end
               %%Generate Hampath!!!
               exitId = ((exitPixs(i,3)-1)*h + exitPixs(i,2)-1)*w + exitPixs(i,1)-1+1;
               hamPath = hamiltonian(Madj, entryId, exitId);
               if isempty(hamPath) % hampath may be invalid!
                   continue;
               end
               [visitOrder, LT] = reconstrPath(hamPath, V, size(V));
                LTG = gradient(LT);
                gradMag = sum(abs(LTG));
%                disp(gradMag)
              
               if gradMag < minCost
                    exitPix = exitPixs(i,:);
                    minCost = gradMag;
                    minpath_id = i;
                    % shall we quit to save time???
                    break;
               end
            end
             % compute visitOrder
         disp('chosen path number:');
         disp(minpath_id);
         if minpath_id == -1 
             disp('HamPath configuration error!');
             disp('Cannot generate exitPixs!');
             needPartition = true;
         end
%          if isempty(hamPath)% Status == 0
%              disp('HamPath Status Error!');
%          end
     else
         disp('Cannot generate exitPixs!');
         needPartition = true;
         visitOrder = [];
         LT = [];
         exitPix = [];
     end
     
     % in case we cannot find a path in vol of 4*4*4, we have to do
     % partitioning.
     if needPartition && w == 4 && h == 4 && d == 4
       [Status, LT, visitOrder, exitPix] = hamPathPartitionGrid3D(V, entryPix, exitFace);
        if Status == 0
            disp('HamPath Status Error!');
        end   
        return;
     end
end

