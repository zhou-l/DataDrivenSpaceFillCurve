%coarse level image file name: clFilename
%next level image file name: nlFilename
function SFCMultiScaleMain(clFilename, nlFilename) 
    close all;
    cI = imread(clFilename);
    cI = double(rgb2gray(cI));
    maxSize = max(size(cI));
    nextPow2 = 2^ceil(log2(maxSize));
    cI = padarray(cI, [(nextPow2-size(cI,1))/2 (nextPow2-size(cI,2))/2], 'both');
    [clLT, clVisitOrder] = SFCmine(cI);
%      clLT = allNodesLinearFunc(clVisitOrder, cI);% visiting all nodes without aggregation

    figure, hold on;
    subplot(2,1,1), plot(1:length(clLT),clLT);title('Current Level SFC');
    
    nI = imread(nlFilename);
    nI = double(rgb2gray(nI));
    maxSize = max(size(nI));
    nextPow2 = 2^ceil(log2(maxSize));
    nI = padarray(nI, [(nextPow2-size(nI,1))/2 (nextPow2-size(nI,2))/2], 'both');
    [nlLT, nlVisitOrder] = refineLTnextLvl(clVisitOrder, cI, nI);
    subplot(2,1,2), plot(1:length(nlLT),nlLT);title('Next Level SFC');
    hold off;
    
    travOrder = zeros(size(nI,1), size(nI,2), 1);
    figure;
    for i = 1:length(nlVisitOrder)
        travOrder(nlVisitOrder(i,1),nlVisitOrder(i,2),:) = i;
    end
    imagesc(travOrder);
end

% To be consistent, we need input of a coarse level image clI, a fine level
% image nlI, and the visit order of the coarse level clVisitOrder
% next level must be 2*2 of the coarse level!
function [nlLT,nlVisitOrder] = refineLTnextLvl(clVisitOrder, clI, nlI)
    szNLI = size(nlI);
    szCLI = size(clI);
    clGraphSize = szCLI(1);
    nlGraphSize = szNLI(1);
    scale = (nlGraphSize/clGraphSize);
    scale2 = scale*scale;
    nlLT = zeros(szNLI(1)*szNLI(2),1);
    nlVisitOrder = zeros(size(nlLT,1)*size(nlLT,2),2);
   
    % generate linearized function by computing space filling curve for the subpicture of nlI for each pixel of clI 
    for i = 1:length(clVisitOrder)
        
            pos = clVisitOrder(i,:);
            %current pos
            y = pos(1); x=pos(2);
            %next pos
            if i < length(clVisitOrder)
                nextPos = clVisitOrder(i+1,:);
                outDir = nextPos - pos;
                if i > 1
                    inDir = pos - clVisitOrder(i-1,:);
                else
                    inDir = pos - [1,0];
                end
            else
                nextPos = clVisitOrder(i-1,:);
                outDir = pos - nextPos;
                inDir = pos - clVisitOrder(i-1,:);
            end
            % get exit edge 
            if norm(outDir) ~= 1
                warning('Error in coarse level SFC! Forced quit!');
                return; 
            end
            
            if outDir(1) == 1 && outDir(2) == 0 % down
                exitEdge = 1;
            elseif outDir(1) == -1 && outDir(2) == 0 %up
                exitEdge = 3;
            elseif outDir(2) == 1 && outDir(1) == 0 %right
                exitEdge = 2;
            elseif outDir(2) == -1 && outDir(1) == 0 %left
                exitEdge = 4;
            end
            
          %------------------------------------------
            % get the sub image that is occupied by [x,y] of the coarse
            % level image
            subI = nlI(scale*(y-1)+1:scale*y, scale*(x-1)+1:scale*x,:);
            % determine inNode position
% %             inNode = findBestInNode(subI, parentSubI, 1);
            if i == 1
                entryPix = [size(subI,1), 1];
            else
                if inDir(1) == 1 %down
                    entryPix(2) = exitPix(2);
                    entryPix(1) = 1;
                elseif inDir(1) == -1 %up
                  entryPix(2) = exitPix(2);
                  entryPix(1) = 2;
                end
                
                if inDir(2) == 1 % right
                    entryPix(1) = exitPix(1);
                    entryPix(2) = 1;
                elseif inDir(2) == -1 %left
                    entryPix(1) = exitPix(1);
                    entryPix(2) = 2;
                end
            end
            % need to use hamilton path
%            [subNlLT, subNlVisitOrder, exitPix] = linearizeHamPath2x2(subI,entryPix,exitEdge);
            % up to scale = 4
           [subNlLT, subNlVisitOrder, exitPix] = linearizeHamPath(subI,entryPix,exitEdge);
           nlLT(scale2*(i-1)+1:scale2*i,:) = subNlLT(:,:);
            nlVisitOrder(scale2*(i-1)+1:scale2*i,:) = [(y-1).* scale + subNlVisitOrder(:,1), (x-1) .* scale +  subNlVisitOrder(:,2)];
    end
end