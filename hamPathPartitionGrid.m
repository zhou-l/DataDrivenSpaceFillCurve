function [Status, LT, visitOrder, exitPix] = hamPathPartitionGrid(I, entryPix, exitEdge)
    imSize = size(I);
    Status = 0;
    LT = [];
    visitOrder = [];
    exitPix = [];
    if mod(imSize(1)*imSize(2),2) ~= 0
        warning('Cannot handle odd sized grid now!');
        return;
    end
    if imSize(1) == 8  && imSize(2) == 8
        % split 8*8 tile: the most encountered case!
        [grids, splitDirs] = partitionGrid(imSize, entryPix, exitEdge);
        if splitDirs == 1
            gridDim = [8 4];
        else
            gridDim = [4 8];
        end
%         switch exitEdge
%         case 1 %bottom
%             if splitDirs == 1 % vertical
%                 if entryPix(2)<=4
%                     intermExitEdge = 2; % right
%                 else
%                     intermExitEdge = 4; % left 
%                 end
%             else
%                 if entryPix(1) <=4
%                     intermExitEdge = 1; % down
%                 else
%                     intermExitEdge = 3; % up
%                 end
%             end
%         case 2 %right
%             if splitDirs == 1 
%                 if entryPix(2) <= 4
%                     intermExitEdge = 2; % right
%                 else
%                     intermExitEdge = 4; % left
%                 end
%             else
% 
%                 if entryPix(1) <=4
%                     intermExitEdge = 1; % down
%                 else
%                     intermExitEdge = 3; % up
%                 end
%             end
%         case 3 %top
%             if splitDirs == 1
%                 if entryPix(2) <=4
%                     intermExitEdge = 2; % right
%                 else
%                     intermExitEdge = 4; % left
%                 end
% %              entryPix = entryPix - [8 4];
%             else
%                 if entryPix(1) <=4
%                     intermExitEdge = 1; % down
%                 else
%                     intermExitEdge = 3; % up
%                 end
%             end
%         case 4 %left
%             if splitDirs == 1
%                 if entryPix(2) <=4
%                     intermExitEdge = 2; % right
%                 else
%                     intermExitEdge = 4; % left
%                 end
% %              entryPix = entryPix - [8 4];
%             else
%                 if entryPix(1) <=4
%                     intermExitEdge = 1; % down
%                 else
%                     intermExitEdge = 3; % up
%                 end
%             end
%         end

        if splitDirs == 1
            if entryPix(2) <=4
                intermExitEdge = 2; % right
            else
                intermExitEdge = 4; % left
            end
%              entryPix = entryPix - [8 4];
        else
            if entryPix(1) <=4
                intermExitEdge = 1; % down
            else
                intermExitEdge = 3; % up
            end
        end
        % first block
        expixs = findCompatibleExitPixs(gridDim, entryPix, intermExitEdge);
        [Status, hamPath1, LT1, intExitPix] = hamPathGridRepNN(I(grids(1,1):grids(1,3), grids(1,2):grids(1,4),:), entryPix, expixs);
        if Status == 0
            warning(Status);
        end
        gridStart1 = grids(1,1:2) - [1,1];
        gridStart2 = grids(2,1:2) - [1,1];
        % second block
        gridDiff = (gridStart2 - gridStart1);
        gridDiff = gridDiff ./ norm(gridDiff);
        % get the transformed entry point
        entryPix2 = intExitPix + gridStart1 + gridDiff;

        expixs2 = findCompatibleExitPixs(gridDim, entryPix2, exitEdge);
        [Status, hamPath2, LT2, exitPix] = hamPathGridRepNN(I(grids(2,1):grids(2,3),grids(2,2):grids(2,4),:), entryPix2, expixs2);
        % connect two subpaths
        if Status == 0
            warning(Status);
        end
        LT = [LT1(:,:); LT2(:,:)];
   
        visitOrder = [gridStart1 + hamPath1(1:end,:); gridStart2 + hamPath2(1:end,:)];
        if norm(visitOrder(32,:) - visitOrder(33,:))>1
            warning('error!');
        end
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


function [grids, splitDirs] = partitionGrid(imSize, entryPix, exitEdge)
    colormat = zeros(imSize(1),imSize(2));
    for i = 1:imSize(1)
        for j = 1:imSize(2)
            s = mod(i + j,2);
            colormat(i,j) = s;
        end
    end
    numGrids = 1;
    splitDirs = [];
%     exitPixs = zeros(imSize(1)/2, 2);% one entry point can have two 
    entryColor = colormat(entryPix(1),entryPix(2));
    cnt = 1;
       switch(exitEdge)
        case 1 % bottom
            exitPixs = zeros(imSize(2)/2, 2);% one entry point can have two 
            exitPixs(:,1) = imSize(1); 

            for i = 1:size(colormat,2)
                if colormat(exitPixs(cnt,1),i)~= entryColor
                    exitPixs(cnt,2) = i;
                    cnt = cnt + 1;
                    if cnt >= size(exitPixs,1)+1
                        break;
                    end
                end
            end
        case 2 % right
            exitPixs = zeros(imSize(1)/2, 2);% one entry point can have two 
            exitPixs(:,2) = imSize(2);
            for i = 1:size(colormat,1)
                if colormat(i,exitPixs(cnt,2))~= entryColor
                    exitPixs(cnt,1) = i;
                    cnt = cnt + 1;
                    if cnt >= size(exitPixs,1)+1
                        break;
                    end
                end
            end
        case 3 % top
             exitPixs = zeros(imSize(2)/2, 2);% one entry point can have two 
           exitPixs(:,1) = 1;
           for i = 1:size(colormat,2)
                if colormat(exitPixs(cnt,1),i)~= entryColor
                    exitPixs(cnt,2) = i;
                    cnt = cnt + 1;
                    if cnt >= size(exitPixs,1)+1
                        break;
                    end
                end
            end
        case 4 % left
            exitPixs = zeros(imSize(1)/2, 2);% one entry point can have two 
            exitPixs(:,2) = 1;
            for i = 1:size(colormat,1)
             if colormat(i,exitPixs(cnt,2))~= entryColor
                    exitPixs(cnt,1) = i;
                    cnt = cnt + 1;
                    if cnt >= size(exitPixs,1)+1
                        break;
                    end
                end
            end   
        end   
    if imSize(1) == 8 && imSize(2) == 8
        numGrids = 2;
        grids = zeros(numGrids, 4); % record top left and bottom right corners of all subgrids
        for  i = 1:length(exitPixs)
            dist = entryPix - exitPixs(i,:);
            if abs(dist(1)) < 4 && abs(dist(2)) < 4
                continue;
            else
                % try vertical cut
                p1 = [1,4]; p2 = [8,4];
                if intersectLL2D(p1,p2,entryPix, exitPixs(i,:))
                    if entryPix(2) > 4
                        grids = [1 5 8 8; 1 1 8 4];
                    else
                        grids = [1 1 8 4; 1 5 8 8];
                    end
                    splitDirs = 1; % vertical split
                else
                    if entryPix(1) > 4
                       grids = [5 1 8 8; 1 1 4 8];
                    else
                       grids = [1 1 4 8; 5 1 8 8];
                    end
                    splitDirs = 2; % horizontal split
                end
                return;
            end
        end
    else
        
    end
     
end