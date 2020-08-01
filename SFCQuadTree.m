function [LT, visitOrder] = SFCQuadTree(Vlvls, QT, entryNode)
LT = [];
visitOrder = [];
% input: QT -- quad tree; Vlvls -- images of all levels
Idim = size(Vlvls{1});
global Isize;
Isize = Idim(1);

maxSize = max(Idim);
nextPow2 = 2^ceil(log2(maxSize));
    
global Vs;
Vs = Vlvls;
global S;
S = QT;
global QTr;
global QTc;
global QTv;
[QTr, QTc, QTv] = find(S);
% get nodes 
Itl = [1,1];
Ibr = [Idim(1),Idim(2)];

outEdge = 4; % should use a flexible out edge!!!
entryPix = [Idim(1),1];

lvl = length(Vlvls);
if size(Vlvls{end},1)<4 % have to use hamPath
    ePix = getEquivPixAtLvl(entryPix, 1, lvl);
    % choose the best outEdge
    minCost = realmax;
     for toEdge = 1:1%4
        [tLTv, tVO, tExitPix] = linearizeHamPath(Vlvls{end}, ePix, toEdge);
        cost = norm(gradient(tLTv)',1); % amaximum absolute row sum
        if cost < minCost
            minCost = cost;
            topLTv = tLTv;
            topVisitOrder = tVO;
            outEdge = toEdge;
            exitPix = tExitPix;
        end
    end
else % could use hamCycle
    [topLTv, topVisitOrder] = SFCmine(Vlvls{end});
     topLTv = topLTv(:,1);
end
figure, hold on;
fTop = flip(Vlvls{end},1);
imagesc(fTop);
line(topVisitOrder(:,2), 2 - topVisitOrder(:,1)+1);
% lineColorCoded(topVisitOrder(:,2), topVisitOrder(:,1), 0, 1:length(topVisitOrder));
hold off;

topLvlBlockSize = lvlToNodeSize(lvl);
topLT(:,1) = topLTv;
topLT(:,2) = repmat(topLvlBlockSize,length(topLTv),1);

inEdge =4; % left
lastVal = -1;
lastPos = [0 0];
LT = [];
visitOrder = topVisitOrder;
idInsert = 1;

for i= 1:length(topVisitOrder)
    
% figure;
% imagesc(Vlvls{1});

    entryN = topVisitOrder(i,:);
    if i == length(topVisitOrder)
        exitN = [1,0]; % a virtual node for the last element
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
        switch outEdge
            case 1
                inEdge = 3;
            case 2 
                inEdge = 4;
            case 3
                inEdge = 1;
            case 4
                inEdge = 2;
        end
    end

    if dir(1) == 1% 1<-bottom 
       outEdge = 1;
    elseif dir(1) == -1 % 3<-top
        outEdge = 3;
    else
        if dir(2) == 1
            outEdge = 2; %2<- right
        else
            outEdge = 4; %4<- left
        end
    end
    
    LTzz = topLT(i,:);
    vOzz = topVisitOrder(i,:);
%     if Itl(1)==65 && Itl(2) == 65
        [epix, LTzz, vOzz] = refine(inEdge, outEdge, Itl, Ibr, lastVal, lastPos);
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
    figure;
    imagesc(Vlvls{1});

    
% plot traverse order at block centers
   zdim = nextPow2;
   ncVisitOrder = visitOrder;
    while(zdim>=1)
        [vals,r,c] = qtgetblk(Vlvls{1},QT,zdim);
        
        if ~isempty(vals)
            for ii = 1:length(r)
                v = [r(ii),c(ii)];
                orderId1 = find(visitOrder(:,1) == v(1));
                orderId2 = find(visitOrder(:,2) == v(2));
                orderId = intersect(orderId1, orderId2);
                if ~isempty(orderId)
                    blkCtr = [r(ii)+zdim/2-0.5,c(ii)+zdim/2-0.5];
%                     orderBlocks(r(ii):r(ii)+dim-1,c(ii):c(i)+dim-1,:) = orderId;
                    ncVisitOrder(orderId,:) = blkCtr;
                end
            end
        end
        zdim = zdim / 2;
        
    end
     line(ncVisitOrder(:,2), ncVisitOrder(:,1));
     hold off;
end
end



%find the entry pixel to the block @ level with the lastVal
function epix = findBestEntry(inEdge, level, lastVal)
    global Vs;
    V = Vs{level};
%     SV = V(blocktl(1):blockrb(1), blocktl(2):blockrb(2));
    % find all pixels along the edge inEdge
    nbl = inEdge(1,:);
    ebl = inEdge(2,:);
    nb = getEquivPixAtLvl(nbl, 1, level);% level);
    eb = getEquivPixAtLvl(ebl, 1, level);% level);
   if nb(1) >= eb(1) && nb(2) >= eb(2)
       epix = eb;
       return;
   end
    minDif = realmax;
    for r = nb(1):eb(1)
        for c = nb(2):eb(2)
            vdif = abs(V(r,c) - lastVal(1));
            if vdif < minDif
                 minDif = vdif;
                 epix = [r,c];
            end
        end
        
    end
    
end

function [Ptl, Pbr] = getNodeExtent(nodePos, lvl)
    blkSize = lvlToNodeSize(lvl);
    Ptl = (nodePos-1) .*blkSize + 1;
    Pbr = Ptl + [blkSize-1,blkSize-1];
end

% get block levels
function [Lfine, Lcoarse, SBr, SBc, SBv] = getLvlsInBlock(blocktl, blockbr)
    global S;
    FS = full(S);
    % get sub image
    subFS = FS(blocktl(1):blockbr(1),blocktl(2):blockbr(2));
%     subS = sparse();
    [SBr, SBc, SBv] = find(subFS);
    Lfine = nodeSizeToLvl(min(SBv));
    Lcoarse = nodeSizeToLvl(max(SBv));
end

function [inEdgeNodes, outEdgeNodes] = getInoutEdgeNodes(blocktl, blockbr, inEdge, outEdge)

    switch inEdge
        case 1 %bottom
            inEdgeNodes(1,:) = [blockbr(1),blocktl(2)];
            inEdgeNodes(2,:) = [blockbr(1),blockbr(2)];
        case 2 %right
            inEdgeNodes(1,:) = [blocktl(1),blockbr(2)];
            inEdgeNodes(2,:) = [blockbr(1),blockbr(2)];
        case 3 %top
            inEdgeNodes(1,:) = [blocktl(1),blocktl(2)];
            inEdgeNodes(2,:) = [blocktl(1),blockbr(2)];
        case 4 %left
            inEdgeNodes(1,:) = [blocktl(1),blocktl(2)];
            inEdgeNodes(2,:) = [blockbr(1),blocktl(2)];
    end
    
     switch outEdge
        case 1 %bottom
            outEdgeNodes(1,:) = [blockbr(1),blocktl(2)];
            outEdgeNodes(2,:) = [blockbr(1),blockbr(2)];
        case 2 %right
            outEdgeNodes(1,:) = [blocktl(1),blockbr(2)];
            outEdgeNodes(2,:) = [blockbr(1),blockbr(2)];
        case 3 %top
            outEdgeNodes(1,:) = [blocktl(1),blocktl(2)];
            outEdgeNodes(2,:) = [blocktl(1),blockbr(2)];
        case 4 %left
            outEdgeNodes(1,:) = [blocktl(1),blocktl(2)];
            outEdgeNodes(2,:) = [blockbr(1),blocktl(2)];
    end
end


function shouldRef = needRefine(blocktl, blockSize)
   [Lfine, Lcoarse, SBr, SBc, SBv] = getLvlsInBlock(blocktl, blocktl+ [blockSize,blockSize] - [1,1]);
   Lblk = nodeSizeToLvl(blockSize);
   if Lcoarse == Lblk
       shouldRef = false;
   else
       shouldRef = true;
   end
end

%% NOTE: the edge direction to edge id functions are based on Matlab conventions: Y is increasing when going down!!!
function edgeId = outDirToEdge(outdir)
        if outdir(1) == 1 % going down
           edgeId = 1;% 1<-bottom
        elseif outdir(1) == -1 % going up
            edgeId =  3;% 3<-top
        else
            if outdir(2) == 1 % going right
                edgeId = 2; %2<- right
            else %going left
                edgeId = 4; %4<- left
            end
        end
end

function edgeId = inDirToEdge(indir)
        if indir(1) == 1 %going down
           edgeId = 3;% 3<-top edge
        elseif indir(1) == -1 %going down
            edgeId = 1;% 1<-bottom edge
        else
            if indir(2) == 1 % going right
                edgeId = 4; %4<- left
            else % dir(2) == -1 going left
                edgeId = 2; %2<- right
            end
        end
end

% recursive function for refinement
function [exitPix, LTinout, VOinout, posInOrg] = refine(inEdge, outEdge, blocktl, blockbr, lastVal, lastNode)
global Vs;
global Isize;
% global LT;
    VOinout = [];
    LTinout = [];
    exitPix = [];
    posInOrg = [];

     blkSize = blockbr - blocktl + 1;
    [Lfine, Lcoarse, SBr, SBc, SBv] = getLvlsInBlock(blocktl, blockbr);
    %%
    % Get the locations of the current block
    if lastNode(1) == 0 && lastNode(2) == 0
        tblocktl = blocktl;
        tblockbr = blockbr;
    else
        lastL = log2(lastVal(2)) + 1;
        currL = log2(blockbr(1)-blocktl(1)+1)+1;
         if lastL < currL % The level of current block is coarser than the previous one
            lastLen = lastVal(2);
            switch inEdge
                case 1 % bottom: y is fixed
                     tblocktl(2) = max(blocktl(2),lastNode(2));
                     tblocktl(1) = blockbr(1) - lastLen + 1;

                     tblockbr(1) = blockbr(1);
                     tblockbr(2) = tblocktl(2)+lastLen - 1; % lastVal(2) := edge length
                case 2 % right: x is fixed
                    tblocktl(2) = blockbr(2) - lastLen + 1;
                    tblocktl(1) = max(blocktl(1), lastNode(1));

                    tblockbr(2) = blockbr(2);
                    tblockbr(1) = tblocktl(1) + lastLen - 1;
                case 3 % top: y is fixed
                    tblocktl(1) = blocktl(1);
                    tblocktl(2) = max(blocktl(2), lastNode(2));

                    tblockbr(1) = blocktl(1) + lastLen - 1;
                    tblockbr(2) = tblocktl(2) + lastLen - 1;

                case 4 % left x is fixed
                    tblocktl(2) = blocktl(2);
                    tblocktl(1) = max(blocktl(1), lastNode(1));

                    tblockbr(2) = blocktl(2) + lastLen - 1;
                    tblockbr(1) = tblocktl(1) + lastLen - 1;
            end
        else
            tblocktl = blocktl;
            tblockbr = blockbr;
        end
    end
    
   
    [inEdgeNodes, outEdgeNodes] = getInoutEdgeNodes(tblocktl, tblockbr, inEdge, outEdge);
    
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
             blocktlInBlock = getEquivPixAtLvl(blocktl, 1, Lfine);
             blockrbInBlock = getEquivPixAtLvl(blockbr, 1, Lfine);
%             LT = Vs{Lfine}(blocktl(1):blockbr(1), blocktl(2):blockbr(2),:);
            %replace item
            LTinout(:,1) = Vs{Lfine}(blocktlInBlock(1):blockrbInBlock(1),...
                blocktlInBlock(2):blockrbInBlock(2),:);
            LTinout(:,2) = repmat(blkSize(1),length(LTinout),1);
            VOinout = blocktl;
            posInOrg = true;
        else
            
            if blkSize(1) > 8
%                  forceRefine = true; % force the block to be smaller than 8*8
                doOneLevelFiner = true;
            else
                ePix = findBestEntry(inEdgeNodes, Lfine, lastVal);
                blocktlInBlock = getEquivPixAtLvl(blocktl, 1, Lfine);
                blockrbInBlock = getEquivPixAtLvl(blockbr, 1, Lfine);
                ePixInBlock = ePix - blocktlInBlock + [1,1];
                [LTK, visitOrderK, exitPix] = linearizeHamPath(Vs{Lfine}(blocktlInBlock(1):blockrbInBlock(1),...
                    blocktlInBlock(2):blockrbInBlock(2)),...
                    ePixInBlock, outEdge);
                LTinout(:,1) = LTK;
                LTinout(:,2) = repmat(lvlToNodeSize(Lfine),length(LTinout),1);
                visitOrderKOrg = getEquivPixAtLvl(visitOrderK - [1,1], Lcoarse, 1) + blocktl;
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
            if  blkSize(1) <= 8
                forceRefine = false;
                LL = Lcoarse;
            else
                LL = Lcoarse;
                solvableBlkSize = blkSize;
                % make sure that a finer level has block size equal to or smaller than 4*4
                while solvableBlkSize(1) > 8
                    LL = LL + 1;
                    solvableBlkSize(1) = solvableBlkSize(1)/2;
                end
            end
        else
            LL = Lcoarse;
        end
        
         ePix = findBestEntry(inEdgeNodes, LL, lastVal);
         blocktlInBlock = getEquivPixAtLvl(blocktl, 1, LL);
         blockrbInBlock = getEquivPixAtLvl(blockbr, 1, LL);
         ePixInBlock = ePix - blocktlInBlock + [1,1];
%          if forceRefine % force one level finer to avoid large blocks for hamPath
%             [LTc, vOc, ePixc] = linearizeHamPath(Vs{LL}(blocktlInBlock(1):blockrbInBlock(1),...
%              blocktlInBlock(2):blockrbInBlock(2)), ePixInBlock, outEdge);
%          end
          if ~forceRefine
             szInBlock = blockrbInBlock-blocktlInBlock + 1;
             if szInBlock(1) > 8
                  forceRefine = true;
                 LL = Lcoarse;
                solvableBlkSize = szInBlock;
                % make sure that a finer level has block size equal to or smaller than 4*4
                while solvableBlkSize(1) > 8
                    LL = LL + 1;
                    solvableBlkSize(1) = solvableBlkSize(1)/2;
                end
                  ePix = findBestEntry(inEdgeNodes, LL, lastVal);
                 blocktlInBlock = getEquivPixAtLvl(blocktl, 1, LL);
                 blockrbInBlock = getEquivPixAtLvl(blockbr, 1, LL);
                ePixInBlock = ePix - blocktlInBlock + [1,1];
             end

         end
            [LTc, vOc, ePixc] = linearizeHamPath(Vs{LL}(blocktlInBlock(1):blockrbInBlock(1),...
             blocktlInBlock(2):blockrbInBlock(2)), ePixInBlock, outEdge);
         % convert back to the original (finest level)
         if(size(vOc,2) ~= 2)
             disp('caution!');
         end
         vOcOrg = getEquivPixAtLvl(vOc - [1,1], LL, 1) + blocktl;
         cblockSize = 2^ (LL - 1);
         cInEdge = inEdge; % adopt the inEdge of the coarser level 
         cOutEdge = outEdge; % adopt the outEdge of the coarser level
         
         LTc(:,1) = LTc;
         LTc(:,2) = repmat(cblockSize,length(LTc),1);
            
         LTz = [];
         VOz = [];
         for i = 1:length(vOcOrg)
               cblocktl = vOcOrg(i,:);
               
               shouldRef = needRefine(cblocktl, cblockSize);
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
                        cOutEdge = outDirToEdge(outdir);
                    else
                        cOutEdge = outEdge;
                    end
                
                    if i > 1
                        inDir = vOc(i,:) - vOc(i-1,:);
                        cInEdge = inDirToEdge(inDir);
                    end
                        
                     

                   % get in, out edges from Nc
                   cblockbr = cblocktl + [cblockSize-1,cblockSize-1];
%                    [NcInEdge, NcOutEdge] = getInoutEdgeNodes(cblocktl, cblockbr,...
%                        cInEdge, cOutEdge);
                   if isempty(VOz)
                       tlastNode = lastNode;
                       tlastVal = lastVal;
                   else
                       tlastNode = VOz(end,:);
                       tlastVal = LTz(end,:);
                   end
                   
                   [ePixNc, LTNc, vOzc, posConverted] = refine(cInEdge, cOutEdge,...
                       cblocktl, cblockbr, tlastVal, tlastNode);
%                    , LTc(i,:), vOc(i,:)
                   if posConverted
                       vOzcOrg = vOzc; 
                   else
                       if vOzc(1)*2^(LL-1)<Isize && vOzc(2)*2^(LL-1)< Isize
                         vOzcOrg = getEquivPixAtLvl(vOzc - [1,1], LL, 1) + blocktl;
                       else
                          vOzcOrg = vOzc;
                       end
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


function epix = getEquivPixAtLvl(entryPix, srcLvl, destLvl)
    epix = ceil(entryPix ./ (2^(destLvl-srcLvl)));
end


function Lvl = checkLevel(node, quadT)
    
    [Nr,Nc,Sv] = find(quadT);
    for i = 1:length(Sv)
        if isInsideBlk(node, [Nr(i),Nc(i)], Sv(i))
            Lvl = nodeSizeToLvl(Sv(i));
            break;
        end
    end
end

% check coarsest node in a sub image given the top left corner IPtl and
% bottom right IPbr
function lvl = checkCoarsestNodeInImg(IPtl, IPbr, quadT)
    S = full(quadT);
    subS = sparse(S(IPtl(1):IPbr(1),IPtl(2):IPbr(2)));
    [blkr,blkc,subSv] = find(subS);
    subSv = sort(subSv, 'descend');
    lvl = nodeSizeToLvl(subSv(1));
end

function isContained = isInsideBlk(pos, blkPos, blkSize)
    isContained = false;
    if pos(1) - blkPos(1) < blkSize && pos(2) - blkPos(2) < blkSize
        isContained = true;
    end
end

function blkSize = lvlToNodeSize(lvl)
   blkSize = 2^(lvl-1);
end

function lvl = nodeSizeToLvl(nodeSize)
   lvl = log2(nodeSize) + 1;
end