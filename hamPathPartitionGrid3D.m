function [Status, LT, visitOrder, exitPix] = hamPathPartitionGrid3D(V, entryPix, exitFace)
    imSize = size(V);
    Status = 0;
    LT = [];
    visitOrder = [];
    exitPix = [];
    if mod(imSize(1)*imSize(2)*imSize(3),2) ~= 0
        warning('Cannot handle odd sized grid now!');
        return;
    end
    if imSize(1) == 4  && imSize(2) == 4 && imSize(3) == 4
        [grids, splitDirs] = partitionGrid3D(imSize, entryPix, exitFace);
        if splitDirs == 1 % xz plane
             if entryPix(2) > 2
                interimExitFace = 1; % go forward, then forward
            else
                interimExitFace = 6; % shouldn't happen really
             end
        elseif splitDirs == 2 % xy plane 
             if entryPix(3) > 2
                interimExitFace = 5; % go down, then forward
            else
                interimExitFace = 3; % go up, then forward
             end
        else % yz plane
            if entryPix(1) > 2 
                interimExitFace = 2; % go left, then forward
            else
                interimExitFace = 4; % go right, then forward
             end
        end

        % first block
        V1 = V(grids(1,1):grids(1,4),grids(1,2):grids(1,5),grids(1,3):grids(1,6));
        [exitPixs,isEntryInExitFace, entryPix1] = findCompatibleExitPixs3D(size(V1), entryPix, interimExitFace);
        Madj = createAdjMatrixVol(V1, size(V1));
        h = size(V1,3); d = size(V1,2); w = size(V1,1);
        entryId = ((entryPix1(3)-1)*d + entryPix1(2)-1)*w + entryPix1(1)-1+1;
        minCost = realmax;
        minpath_id = -1;
        LT1 = [];
        for i = size(exitPixs,1):-1:1
            if exitPixs(i,1) == entryPix1(1) && exitPixs(i,2) == entryPix1(2) && exitPixs(i,3) == entryPix1(3)
                continue;
            end
           exitId = ((exitPixs(i,3)-1)*d + exitPixs(i,2)-1)*w + exitPixs(i,1)-1+1;
           hamPath = hamiltonian(Madj, entryId, exitId);
           if isempty(hamPath) % hampath may be invalid!
               continue;
           end
           [svisitOrder, sLT] = reconstrPath(hamPath, V1, size(V1));
            LTG = gradient(sLT);
            gradMag = sum(abs(LTG));
%            disp(gradMag)

           if gradMag < minCost
                inExitPix = exitPixs(i,:);
                minCost = gradMag;
                minpath_id = i;
                visitOrder1 = svisitOrder;
                LT1 = sLT;
                Status = 1;
           end
        end
         % compute visitOrder

    %             [Status, visitOrder, LT, exitPix] = hamPathRepNN(V, entryPix, exitPixs); 
         disp('chosen path number:');
         disp(minpath_id);
         
        if Status == 0
            disp('Hampath cannot be created!');
            disp(size(V1));
            disp(entryPix1);
            return;
        end
        gridStart1 = grids(1,1:3) - [1,1,1];
        gridStart2 = grids(2,1:3) - [1,1,1];
        % second block
        gridDiff = (gridStart2 - gridStart1);
        gridDiff = gridDiff ./ norm(gridDiff);
        % get the transformed entry point
        inPix2 = inExitPix + gridStart1 + gridDiff;
        V2 = V(grids(2,1):grids(2,4),grids(2,2):grids(2,5),grids(2,3):grids(2,6));
        
        h = size(V2,3); d = size(V2,2); w = size(V2,1);
        [exitPixs,isEntryInExitFace, entryPix2] = findCompatibleExitPixs3D(size(V2), inPix2, exitFace);
        Madj = createAdjMatrixVol(V2, size(V2));
        entryId = ((entryPix2(3)-1)*d + entryPix2(2)-1)*w + entryPix2(1)-1+1;
        minCost = realmax;
        minpath_id = -1;
        LT2 =[];
        Status = 0;
        visitOrder2 = [];
        for i = size(exitPixs,1):-1:1
            if exitPixs(i,1) == entryPix2(1) && exitPixs(i,2) == entryPix2(2) && exitPixs(i,3) == entryPix2(3)
                continue;
            end
           exitId = ((exitPixs(i,3)-1)*d + exitPixs(i,2)-1)*w + exitPixs(i,1)-1+1;
           hamPath = hamiltonian(Madj, entryId, exitId);
           if isempty(hamPath) % hampath may be invalid!
               continue;
           end
           [svisitOrder, sLT] = reconstrPath(hamPath, V2, size(V2));
            LTG = gradient(sLT);
            gradMag = sum(abs(LTG));
%            disp(gradMag)

           if gradMag < minCost
                exitPix = exitPixs(i,:);
                minCost = gradMag;
                minpath_id = i;
                visitOrder2 = svisitOrder;
                LT2 = sLT;
                Status = 1;
           end
        end
         % compute visitOrder
        if Status == 0
            disp('Hampath cannot be created!');
            disp(size(V2));
            disp(entryPix2);
            return;
        end
%             [Status, visitOrder, LT, exitPix] = hamPathRepNN(V, entryPix, exitPixs); 
%      disp('chosen path number:');
%      disp(minpath_id);
%         [LT2, visitOrder2, exitPix] = linearizeHamPath3D(V(grids(2,1):grids(2,4),grids(2,2):grids(2,5),grids(2,3):grids(2,6)),...
%             entryPix2, exitFace);

        if isempty(LT2) || isempty(visitOrder2)
            disp('Something went wrong during Hampath generation.');
            disp(size(V)); 
            disp(entryPix);
             return;
        end
       disp('chosen path number:');
        disp(minpath_id);
        LT = [LT1;LT2];
        visitOrder = [gridStart1 + visitOrder1; gridStart2 + visitOrder2];
%         
%         expixs2 = findCompatibleExitPixs(gridDim, entryPix2, exitFace);
%         [Status, hamPath2, LT2, exitPix] = hamPathGridRepNN(V(grids(2,1):grids(2,3),grids(2,2):grids(2,4),:), entryPix2, expixs2);
%         % connect two subpaths
%         if Status == 0
%             warning(Status);
%         end
%         LT = [LT1(:,:); LT2(:,:)];
%    
%         visitOrder = [gridStart1 + hamPath1(1:end,:); gridStart2 + hamPath2(1:end,:)];
%         if norm(visitOrder(32,:) - visitOrder(33,:))>1
%             warning('error!');
%         end
%         if splitDirs == 1 % vertical split
%             visitOrder = [hamPath1(1:end,:); grids(2,1:2) -[1,1] + hamPath2(1:end,:)];
%         else
%             visitOrder = [hamPath1(1:end,:); [4 0] + hamPath2(1:end,:)];
%         end
        
    else
        % a more general case
        disp('cannot handle yet!');
    end
end 


function [grids, splitDirs] = partitionGrid3D(imSize, entryPix, exitFace)
 
    splitDirs = [];
%     exitPixs = zeros(imSize(1)/2, 2);% one entry point can have two 
    exitPixs = findCompatibleExitPixs3D(imSize, entryPix, exitFace);
    if imSize(1) == 4 && imSize(2) == 4 && imSize(3) == 4
        numGrids = 2;
        grids = zeros(numGrids, 6); % record min and max corners of all subgrids (3D)
        for  i = 1:length(exitPixs)
            dist = entryPix - exitPixs(i,:);
            if abs(dist(1)) < 2 && abs(dist(2)) < 2 && abs(dist(3)) < 2 % shorter than have the length?
                continue;
            else
                % try vertical cut
                plane1 = [1 2 1; 0 1 0]; % xz plane
                plane2 = [1 2 2; 0 0 1]; % xy plane
                plane3 = [2 1 1; 1 0 0]; % yz plane
                l_dir = entryPix - exitPixs(i,:) ./ norm(entryPix - exitPixs(i,:));
                l_p = entryPix;
                % intersection test plane 1
                [pt_i,rc] = line_plane_intersection(l_dir, l_p, plane1(4:6), plane1(1:3));
                if rc == 1 %only inersection at 1 point counts!
                    if entryPix(2) > 2
                        grids = [1 3 1 4 4 4; 1 1 1 4 2 4];
                    else
                        grids = [1 1 1 4 2 4; 1 3 1 4 4 4];
                    end
                    splitDirs = 1; %split at xz plane
                    return;
                end 
                % intersection test plane 2
                [pt_i,rc] = line_plane_intersection(l_dir, l_p, plane2(4:6), plane2(1:3));
                if rc == 1 %only inersection at 1 point counts!
                    if entryPix(1) > 2
                        grids = [1 1 3 4 4 4; 1 1 1 4 4 2];
                    else
                        grids = [1 1 1 4 4 2; 1 1 3 4 4 4];
                    end
                    splitDirs = 2; %split at xy plane
                    return;
                end 
                % intersection test plane 3
                [pt_i,rc] = line_plane_intersection(l_dir, l_p, plane3(4:6), plane3(1:3));
                if rc == 1 %only inersection at 1 point counts!
                    if entryPix(3) > 2
                        grids = [3 1 1 4 4 4; 1 1 1 2 4 4];
                    else
                        grids = [1 1 1 2 4 4; 3 1 1 4 4 4];
                    end
                    splitDirs = 3; %split at yz plane
                    return;
                end      
               
            end
        end
    else
        
    end
     
end