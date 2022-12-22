
function [LT, visitOrder] = SFCOctTreePoint(Vlvls, OT, blockSizeFactor, minBlockSize)
LT = [];
visitOrder = [];
Idim = size(Vlvls{1});
global Isize;
Isize = Idim(1);
global Vs;
Vs = Vlvls;
global Tree;
Tree = OT;
global bsFactor;
bsFactor = blockSizeFactor;
global minBlkSize;
minBlkSize = minBlockSize;

%% the mode of octree: 1 -- Volume Octree; 2 -- Point Octree
global modeOCT;
modeOCT = OT.mode;
% get nodes 
Itl = [1,1,1];
Ibr = [Idim(1),Idim(2),Idim(3)];
outFace = 2; % should use a flexible out edge!!! Going out to the left face
entryPix = [1,1,1];

% SHOULD use encoding of the 3D hamCycle version!
% node configuration:
               %   p6   p7
               %  /|   /|
               % p1   p2|
               % | p5-|p8
               % p4---p3
% face configuration
               %      
               %  /|F3   /
               % |----|F4| 
               % | F1 | /z
               % |----|/
% F1: Front 
% F2: Left
% F3: Top
% F4: Right
% F5: Bottom
% F6: Back
              
lvl = length(Vlvls);
if size(Vlvls{end},1)<=4 % have to use hamPath
    ePix = getEquivPixAtLvl(entryPix, 1, lvl);
    % choose the best outEdge
    minCost = realmax;
     for toFace = 1:6
        [tLTv, tVO, tExitPix] = linearizeHamPath3D(Vlvls{end}, ePix, toFace);
        cost = norm(gradient(tLTv)',1); % amaximum absolute row sum
       
        if cost < minCost
            minCost = cost;
            topLTv = tLTv;
            topVisitOrder = tVO;
            outFace = toFace;
            exitPix = tExitPix;
        end
    end
else % could use hamCycle
    [topLTv, topVisitOrder] = SFCmine3D(Vlvls{end});
     topLTv = topLTv(:,1);
end
topLvlBlockSize = lvlToNodeSize(lvl);
topLT(:,1) = topLTv;
topLT(:,2) = repmat(topLvlBlockSize,length(topLTv),1);

inFace =2; % left
lastVal = -1;
lastPos = [0 0 0];
LT = [];
visitOrder = topVisitOrder;
idInsert = 1;

for i= 1:length(topVisitOrder)
    entryN = topVisitOrder(i,:);
    if i == length(topVisitOrder)
        exitN = [1,1,0]; % a virtual node for the last element
    else
        exitN = topVisitOrder(i+1,:);
    end
    
    [Itl,Ibr] = getNodeExtent(entryN, lvl);
    dir = exitN - entryN;
    if i ~= length(topVisitOrder) && norm(dir) > 1
        warning('Error: Distance between two nodes in path is greater than 1!');
        disp(i);
        return;
    end
    if i > 1
        switch outFace
            case 1 % front
                inFace = 6; % back
            case 2 % left
                inFace = 4; % right
            case 3 % top
                inFace = 5; % bottom
            case 4 % right
                inFace = 2; % left
            case 5 % bottom 
                inFace = 3; % top
            case 6 
        end
    end

    if dir(3) == -1% 5<-bottom 
       outFace = 5;
    elseif dir(3) == 1 % 3<-top
        outFace = 3;
    else
        if dir(1) == 1
            outFace = 4; %4<- right
        elseif dir(1) == -1
            outFace = 2; %2<- left
        else
            if dir(2) == 1
                outFace = 6; %6 <- back
            else
                outFace = 1; %1 <- front
            end
        end
    end
    
    LTzz = topLT(i,:);
    vOzz = topVisitOrder(i,:);
%     if Itl(1)==65 && Itl(2) == 65
    strProg = sprintf('Progress: top level block %i of %i', i, length(topVisitOrder));
    disp(strProg);
   [epix, LTzz, vOzz] = refine3D(inFace, outFace, Itl, Ibr, lastVal, lastPos);
%     end
    lastVal = LTzz(end,:);
    lastPos = vOzz(end,:);
    
    % set return parameters
    
%     idInsert = length(LTzz(:,1)) + idInsert - 1; %length(LT) - length(topLT) + i + 1;
    if i > 1
        LT = [LT; LTzz];
        visitOrder = [visitOrder; vOzz];
    else
        LT = LTzz;
%         LT = [LTzz;topLT(i+1:end,:)];
        visitOrder = vOzz;
    end
%    disp(LT);
end
end

function faceId = outDirToFace(outdir)

    if outdir(3) == -1% 5<-bottom 
       faceId = 5;
    elseif outdir(3) == 1 % 3<-top
        faceId = 3;
    else
        if outdir(1) == 1
            faceId = 4; %4<- right
        elseif outdir(1) == -1
            faceId = 2; %2<- left
        else
            if outdir(2) == 1
                faceId = 6; %6 <- back
            else
                faceId = 1; %1 <- front
            end
        end
    end
end

function faceId = inDirToFace(indir)
    if indir(3) == 1% 3<-top 
       faceId = 3;
    elseif indir(3) == -1 
        faceId = 5;% 5<-bottom
    else
        if indir(1) == 1
            faceId = 2; %2<- left
        elseif indir(1) == -1
            faceId = 4; %4<- right
        else
            if indir(2) == 1
                faceId = 1; %1 <- front
            else
                faceId = 6; %6 <- back
            end
        end
    end
end

function [exitPix, LTinout, VOinout, posInOrg] = refine3D(inFace, outFace, blockmin, blockmax, lastVal, lastNode)
global Vs;
global Isize;
% global LT;
    VOinout = [];
    LTinout = [];
    exitPix = [];
    posInOrg = [];

     blkSize = blockmax - blockmin + 1;
    [Lfine, Lcoarse] = getLvlsInBlock3D(blockmin, blockmax);
    %% Get the starting position of the current block
    if lastNode(1) == 0 && lastNode(2) == 0 && lastNode(3) == 0
        tblockmin = blockmin;
        tblockmax = blockmax;
    else
        lastL = log2(lastVal(2)) + 1;
        currL = log2(blockmax(1)-blockmin(1)+1)+1;
        %% need to think the correct version: The first node to be visited should have the same level as the last node 
         if lastL < currL && Lcoarse < lastL % the last node has a smaller size than the current one
            lastLen = lastVal(2);
            switch inFace
                case 1 % front: y is fixed
                    tblockmin(2) = blockmin(2);
                    tblockmin(1) = max(blockmin(1), lastNode(1));
                    tblockmin(3) = max(blockmin(3),lastNode(3));

                    tblockmax(2) = blockmin(2) + lastLen - 1;
                    tblockmax(1) = tblockmin(1) + lastLen - 1;
                    tblockmax(3) = tblockmin(3) + lastLen - 1;
      
                case 4 % right: x is fixed
                    tblockmin(1) = blockmax(1) - lastLen + 1;
                    tblockmin(2) = max(blockmin(2), lastNode(2));
                    tblockmin(3) = max(blockmin(3),lastNode(3));

                    tblockmax(1) = tblockmin(1) + lastLen - 1;
                    tblockmax(2) = tblockmin(2) + lastLen - 1;
                    tblockmax(3) = tblockmin(3) + lastLen - 1;

                case 6 % back: y is fixed
                    tblockmin(2) = blockmax(2) - lastLen + 1;
                    tblockmin(1) = max(blockmin(1), lastNode(1));
                    tblockmin(3) = max(blockmin(3),lastNode(3));

                    tblockmax(2) = blockmax(2);
                    tblockmax(1) = tblockmin(1) + lastLen - 1;
                    tblockmax(3) = tblockmin(3) + lastLen - 1;
                    
                case 2 % left: x is fixed
                     tblockmin(1) = blockmin(1);
                     tblockmin(2) = max(blockmin(2),lastNode(2));
                     tblockmin(3) = max(blockmin(3),lastNode(3));
                     
%                      tblockmax(3) = 
                     tblockmax(1) = blockmax(1);
                     tblockmax(2) = tblockmin(2)+lastLen - 1; % lastVal(2) := edge length
                     tblockmax(3) = tblockmin(3) + lastLen - 1;

                case 5 % bottom: z is fixed        
                    tblockmin(3) = blockmin(3);
                    tblockmin(2) = max(blockmin(2), lastNode(2));
                    tblockmin(1) = max(blockmin(1), lastNode(1));
                    
                    tblockmax(2) = tblockmin(2) + lastLen - 1;
                    tblockmax(1) = tblockmin(1) + lastLen - 1;
                    tblockmax(3) = tblockmin(3) + lastLen - 1;
                case 3 % top: z is fixed
                    tblockmin(3) = blockmax(3) - lastLen + 1;
                    tblockmin(2) = max(blockmin(2), lastNode(2));
                    tblockmin(1) = max(blockmin(1), lastNode(1));

                    tblockmax(2) = tblockmin(2) + lastLen - 1;
                    tblockmax(1) = tblockmin(1) + lastLen - 1;
                    tblockmax(3) = tblockmin(3) + lastLen - 1; 
            end
            if tblockmin(1) < blockmin(1) || tblockmin(2) < blockmin(2) || tblockmin(3) < blockmin(3) ||...
                tblockmax(1) > blockmax(1)|| tblockmax(2) > blockmax(2) || tblockmax(3) > blockmax(3)    
                tblockmin = blockmin;
                tblockmax = blockmax;
                disp('Potential Erroneous Case');
                disp(inFace);
            end
        else
            tblockmin = blockmin;
            tblockmax = blockmax;
          end
    end
    
   
    [inFaceNodes, outFaceNodes] = getInoutFaceNodes(tblockmin, tblockmax, inFace, outFace);
    
    doOneLevelFiner = false;
    forceRefine = false;
    posInOrg = false;
        % check levels in the block
%     if Lcoarse - Lfine >= 6
%         disp(Lcoarse-Lfine);
%         warning('danger! too many levels!');
%     end
    if Lfine == Lcoarse % only one level in the block
        % if blockSize == blockSize @ L
        if blkSize(1) == lvlToNodeSize(Lfine)
            % get value from Lcoarse == Lfine
             blockminInBlock = getEquivPixAtLvl(blockmin, 1, Lfine);
             blockrbInBlock = getEquivPixAtLvl(blockmax, 1, Lfine);
%             LT = Vs{Lfine}(blockmin(1):blockmax(1), blockmin(2):blockmax(2),:);
            %replace item
            LTinout(:,1) = Vs{Lfine}(blockminInBlock(1):blockrbInBlock(1),...
                blockminInBlock(2):blockrbInBlock(2),blockminInBlock(3):blockrbInBlock(3),:);
            LTinout(:,2) = repmat(blkSize(1),length(LTinout),1);
            VOinout = blockmin;
            posInOrg = true;
        else
            
            if blkSize(1) > 4
%                  forceRefine = true; % force the block to be smaller than 8*8
% should we limit the 3D block to be smaller than 4*4*4?
                doOneLevelFiner = true;
            else
                ePix = findBestEntry(inFaceNodes, Lfine, lastVal);
                blockminInBlock = getEquivPixAtLvl(blockmin, 1, Lfine);
                blockrbInBlock = getEquivPixAtLvl(blockmax, 1, Lfine);
                ePixInBlock = ePix - blockminInBlock + [1,1,1];
                if blockrbInBlock == blockminInBlock
                    visitOrderK = blockminInBlock;
                    LTK = Vs{Lfine}(blockminInBlock(1),blockminInBlock(2),blockminInBlock(3),:);
                    exitPix = ePixInBlock;
                else
                % hamPath 3D
                [LTK, visitOrderK, exitPix] = linearizeHamPath3D(Vs{Lfine}(blockminInBlock(1):blockrbInBlock(1),...
                    blockminInBlock(2):blockrbInBlock(2),blockminInBlock(3):blockrbInBlock(3)), ePixInBlock, outFace);
%                 [LTK, visitOrderK, exitPix] = linearizeHamPath(Vs{Lfine}(blockminInBlock(1):blockrbInBlock(1),...
%                     blockminInBlock(2):blockrbInBlock(2)),...
%                     ePixInBlock, outEdge);
                end
                if isempty(LTK)
                    disp('Something went wrong in hamPath generation for the block.');
                end
                LTinout(:,1) = LTK;
                LTinout(:,2) = repmat(lvlToNodeSize(Lfine),length(LTinout),1);
                visitOrderKOrg = getEquivPixAtLvl(visitOrderK - [1,1,1], Lcoarse, 1) + blockmin;
                VOinout = visitOrderKOrg;
                posInOrg = true;
            end
        end
    else
        doOneLevelFiner = true;
    end
    
    if doOneLevelFiner % more than one levels in the block
        % work with a level finer
        if forceRefine 
            if  blkSize(1) <= 4
                forceRefine = false;
                LL = Lcoarse;
            else
                LL = Lcoarse;
                solvableBlkSize = blkSize;
                % make sure that a finer level has block size equal to or smaller than 4*4
                while solvableBlkSize(1) > 4
                    LL = LL + 1;
                    solvableBlkSize(1) = solvableBlkSize(1)/2;
                end
            end
        else
            LL = Lcoarse;
        end
        
         ePix = findBestEntry(inFaceNodes, LL, lastVal);
         blockminInBlock = getEquivPixAtLvl(blockmin, 1, LL);
         blockrbInBlock = getEquivPixAtLvl(blockmax, 1, LL);
         ePixInBlock = ePix - blockminInBlock + [1,1,1];
%          if forceRefine % force one level finer to avoid large blocks for hamPath
%             [LTc, vOc, ePixc] = linearizeHamPath(Vs{LL}(blockminInBlock(1):blockrbInBlock(1),...
%              blockminInBlock(2):blockrbInBlock(2)), ePixInBlock, outEdge);
%          end
          if ~forceRefine
             szInBlock = blockrbInBlock-blockminInBlock + 1;
             if szInBlock(1) > 4
                  forceRefine = true;
                 LL = Lcoarse;
                solvableBlkSize = szInBlock;
                % make sure that a finer level has block size equal to or smaller than 4*4
                while solvableBlkSize(1) > 4
                    LL = LL + 1;
                    solvableBlkSize(1) = solvableBlkSize(1)/2;
                end
                  ePix = findBestEntry(inFaceNodes, LL, lastVal);
                 blockminInBlock = getEquivPixAtLvl(blockmin, 1, LL);
                 blockrbInBlock = getEquivPixAtLvl(blockmax, 1, LL);
                ePixInBlock = ePix - blockminInBlock + [1,1,1];
             end

          end
         [LTc, vOc, ePixc] = linearizeHamPath3D(Vs{LL}(blockminInBlock(1):blockrbInBlock(1),...
                    blockminInBlock(2):blockrbInBlock(2),blockminInBlock(3):blockrbInBlock(3)), ePixInBlock, outFace);
%             [LTc, vOc, ePixc] = linearizeHamPath(Vs{LL}(blockminInBlock(1):blockrbInBlock(1),...
%              blockminInBlock(2):blockrbInBlock(2)), ePixInBlock, outEdge);
         % convert back to the original (finest level)
         if(size(vOc,2) ~= 2)
             disp('caution!');
         end
         vOcOrg = getEquivPixAtLvl(vOc - [1,1,1], LL, 1) + blockmin;
         cblockSize = 2^ (LL - 1);
         cInFace = inFace; % adopt the inEdge of the coarser level 
         cOutFace = outFace; % adopt the outEdge of the coarser level
         
         LTc(:,1) = LTc;
         LTc(:,2) = repmat(cblockSize,length(LTc),1);
            
         LTz = [];
         VOz = [];
         for i = 1:length(vOcOrg)
               cblockmin = vOcOrg(i,:);
               
               shouldRef = needRefine(cblockmin, cblockSize);
               if forceRefine
                   shouldRef = true;
               end
               if shouldRef
                    entryN = vOc(i,:);
%                     updateEdges();
                    if i ~= length(vOc)
                        % try every possible direction?
                        %TODO: Choose the right direction!
                        exitN = vOc(i+1,:);
                        outdir = exitN - entryN;
                        if norm(outdir) > 1
                            warning('Error: Distance between two nodes in path is greater than 1!');
                            return;
                        end
                        % get edge id
                        cOutFace = outDirToFace(outdir);
                    else
                        cOutFace = outFace;
                    end
                
                    if i > 1
                        inDir = vOc(i,:) - vOc(i-1,:);
                        cInFace = inDirToFace(inDir);
                    end
                        
                     

                   % get in, out edges from Nc
                   cblockmax = cblockmin + [cblockSize-1,cblockSize-1,cblockSize-1];
%                    [NcInEdge, NcOutEdge] = getInoutEdgeNodes(cblockmin, cblockmax,...
%                        cInEdge, cOutEdge);
                   if isempty(VOz)
                       tlastNode = lastNode;
                       tlastVal = lastVal;
                   else
                       tlastNode = VOz(end,:);
                       tlastVal = LTz(end,:);
                   end
                   
                   [ePixNc, LTNc, vOzc, posConverted] = refine3D(cInFace, cOutFace,...
                       cblockmin, cblockmax, tlastVal, tlastNode);
%                    , LTc(i,:), vOc(i,:)
                   if posConverted
                       vOzcOrg = vOzc; 
                   else
                       vOzcOrg = vOzc;
%                        if vOzc(1)*2^(LL-1)<Isize && vOzc(2)*2^(LL-1)< Isize && vOzc(3)*2^(LL-1) < Isize
%                          vOzcOrg = getEquivPixAtLvl(vOzc - [1,1,1], LL, 1) + blockmin;
%                        else
%                           vOzcOrg = vOzc;
%                        end
                   end
                   % replace LTNc, vOc items: 
                    if i > 1
                       % keep the last i+1:end items
                       LTz = [LTz; LTNc];
%                        VOz = [VOz(1:idInsert-1,:); vOzc(1:end,:); vOc(i+1:end,:)];
                     
                       VOz = [VOz; vOzcOrg];
                      
                   else
%                        LTz = [LTNc(1:end,:); LTc(i+1:end,:)];
                        LTz = LTNc;
%                        VOz = [vOzc(1:end,:); vOc(i+1:end,:)];
                       VOz = vOzcOrg;
                    end
                    posInOrg = true;
               else
                   if i > length(LTc)
                       disp('index out of bound!');
                       return;
                   end
                   LTz = [LTz;LTc(i,:)];
                   VOz = [VOz;vOcOrg(i,:)];
               end
               
               if i == length(vOc)
                   if shouldRef
                       exitPix = ePixNc;
                   else
                       exitPix = ePixc;
                   end
               end
         end
         LTinout = LTz; 
         VOinout = VOz;
    end
    % need to convert BACK to locations in the finest level
    exitPix = VOinout(end,:);
end

function [Ptl, Pbr] = getNodeExtent(nodePos, lvl)
    blkSize = lvlToNodeSize(lvl);
    Ptl = (nodePos-1) .*blkSize + 1;
    Pbr = Ptl + [blkSize-1,blkSize-1,blkSize-1];
end

% get block levels
function [Lfine, Lcoarse] = getLvlsInBlock3D(blocktl, blockbr)
    global Vs;
    FS = Vs{1};
    global Tree;
    global bsFactor;
    global minBlkSize;
    % get sub image
%     subFS = FS(blocktl(1):blockbr(1),blocktl(2):blockbr(2),blocktl(3):blockbr(3));
    % find all nodes within the subimage
    % TODO:
    lvls = [];
    blks = [];
     blkSize = blockbr - blocktl + [1 1 1];
%     blkSize = blockbr - blocktl;
    octBlkSize =  blkSize(1);
    octBlktl = (blocktl - [1 1 1]);
    octBlkbr = octBlktl + octBlkSize;
    
    [binNos,minD,maxD] = Tree.rangeQuery([octBlktl octBlkbr]);
    
%     for i = 1:Tree.BinCount
%          boundsToCheck = Tree.BinBoundaries(i,:);
%           hasChildren = find(Tree.BinParents == i,1);
%         if isempty(hasChildren) 
%          if boundsToCheck(1) >= octBlktl(1) && boundsToCheck(2)>=octBlktl(2) && boundsToCheck(3)>=octBlktl(3)...
%              && boundsToCheck(4) <= octBlkbr(1) &&...
%                 boundsToCheck(5) <= octBlkbr(2) && boundsToCheck(6) <= octBlkbr(3)
%              % has to be leaf node!!!
%              lvls(end+1,:) = Tree.BinDepths(i);
%              blks(end+1,:) = boundsToCheck;
%          end
%         end
%      end
    %% NOTE: The OcTree library uses 0 for coarsest level !!!! Other parts of our method uses an inverted level notation!!!
    % check min/max nodes within the sub volume
%     [SBr, SBc, SBv] = find(subFS);
    if isempty(binNos) % In fact, this should not happen
%     if isempty(lvls)
        Lfine = int16(ceil(nodeSizeToLvl(octBlkSize(1))));%size(Vs,1);%1;
        Lcoarse = Lfine;%size(Vs,1);%1;
    else
%         maxD = max(lvls);
%         minD = min(lvls);
      Lfine =  convertOctLvlToSFCLvl(maxD, Tree.DepthMax);%    size(Vs,1) - (max(lvls)) + 1;
      Lcoarse = convertOctLvlToSFCLvl(minD,Tree.DepthMax); %size(Vs,1) -  (min(lvls)) + 1;    
     end
end
function [inFaceNodes, outFaceNodes] = getInoutFaceNodes(blockmin, blockmax, inFace, outFace)

    P4 = blockmin;
    P7 = blockmax;
    P1 = [blockmin(1),blockmin(2),blockmax(3)];
    P2 = [blockmax(1),blockmin(2),blockmax(3)];
    P3 = [blockmax(1),blockmin(2),blockmin(3)];
    P5 = [blockmin(1),blockmax(2),blockmin(3)];
    P6 = [blockmin(1),blockmax(2),blockmax(3)];
    P8 = [blockmax(1),blockmax(2),blockmin(3)];
    
    
    switch inFace
        case 1 %front P1,P2,P3,P4
            inFaceNodes(1,:) = P1;
            inFaceNodes(2,:) = P2;
            inFaceNodes(3,:) = P3;
            inFaceNodes(4,:) = P4;
        case 2 %left P1,P4,P5,P6
            inFaceNodes(1,:) = P1;
            inFaceNodes(2,:) = P4;
            inFaceNodes(3,:) = P5;
            inFaceNodes(4,:) = P6;
        case 3 %top P1,P2,P7,P6
            inFaceNodes(1,:) = P1;
            inFaceNodes(2,:) = P2;
            inFaceNodes(3,:) = P7;
            inFaceNodes(4,:) = P6;
        case 4 %right P2,P3,P8,P7
            inFaceNodes(1,:) = P2;
            inFaceNodes(2,:) = P3;
            inFaceNodes(3,:) = P8;
            inFaceNodes(4,:) = P7;
        case 5 %bottom P4,P3,P8,P5
            inFaceNodes(1,:) = P4;
            inFaceNodes(2,:) = P3;
            inFaceNodes(3,:) = P8;
            inFaceNodes(4,:) = P5;
        case 6 %bottom P6,P7,P8,P5
            inFaceNodes(1,:) = P6;
            inFaceNodes(2,:) = P7;
            inFaceNodes(3,:) = P8;
            inFaceNodes(4,:) = P5;
    end
      switch outFace
      case 1 %front P1,P2,P3,P4
            outFaceNodes(1,:) = P1;
            outFaceNodes(2,:) = P2;
            outFaceNodes(3,:) = P3;
            outFaceNodes(4,:) = P4;
        case 2 %left P1,P4,P5,P6
            outFaceNodes(1,:) = P1;
            outFaceNodes(2,:) = P4;
            outFaceNodes(3,:) = P5;
            outFaceNodes(4,:) = P6;
        case 3 %top P1,P2,P7,P6
            outFaceNodes(1,:) = P1;
            outFaceNodes(2,:) = P2;
            outFaceNodes(3,:) = P7;
            outFaceNodes(4,:) = P6;
        case 4 %right P2,P3,P8,P7
            outFaceNodes(1,:) = P2;
            outFaceNodes(2,:) = P3;
            outFaceNodes(3,:) = P8;
            outFaceNodes(4,:) = P7;
        case 5 %bottom P4,P3,P8,P5
            outFaceNodes(1,:) = P4;
            outFaceNodes(2,:) = P3;
            outFaceNodes(3,:) = P8;
            outFaceNodes(4,:) = P5;
        case 6 %bottom P6,P7,P8,P5
            outFaceNodes(1,:) = P6;
            outFaceNodes(2,:) = P7;
            outFaceNodes(3,:) = P8;
            outFaceNodes(4,:) = P5; 
    end 

end

function epix = getEquivPixAtLvl(entryPix, srcLvl, destLvl)
epix = ceil(entryPix ./ double(2^(destLvl-srcLvl)));
end

function blkSize = lvlToNodeSize(lvl)
blkSize = 2^(lvl-1);
end

%find the entry pixel to the block @ level with the lastVal
function epix = findBestEntry(inFaceNodes, level, lastVal)
    global Vs;
    epix = [];
    disp(level);
    if isempty(level)
        disp('Level error!');
        return;
    elseif level < 1 || level > size(Vs,1)
        disp('Level out of bound');
        return;
    end
    V = Vs{level};
%     SV = V(blocktl(1):blockrb(1), blocktl(2):blockrb(2));
    % find all pixels along the face inFace
    % get min max nodes from the face nodes
    nbl = min(inFaceNodes(:,:));
    ebl = max(inFaceNodes(:,:));

    nb = getEquivPixAtLvl(nbl, 1, level);% level);
    eb = getEquivPixAtLvl(ebl, 1, level);% level);
   if nb(1) >= eb(1) && nb(2) >= eb(2) && nb(3) >= eb(3)
       epix = eb;
       return;
   end
    minDif = realmax;
    for z = nb(3):eb(3)
        for y = nb(1):eb(1)
            for x = nb(2):eb(2)
                vdif = abs(V(y,x,z) - lastVal(1));
                if vdif < minDif
                     minDif = vdif;
                     epix = [y,x,z];
                end
            end
        end
    end
    
end


function lvl = nodeSizeToLvl(nodeSize)
   lvl = log2(nodeSize) + 1;
end

function shouldRef = needRefine(blocktl, blockSize)
 global modeOCT;
%  if modeOCT == 2 % point-based octree
%     [Lfine, Lcoarse] = getLvlsInBlock3D(blocktl, blocktl+ [blockSize,blockSize,blockSize]);
%  else
   [Lfine, Lcoarse] = getLvlsInBlock3D(blocktl, blocktl+ [blockSize,blockSize,blockSize] - [1,1,1]);
%  end
   Lblk = nodeSizeToLvl(blockSize);
   if Lcoarse == Lblk
       shouldRef = false;
   else
       shouldRef = true;
   end
end