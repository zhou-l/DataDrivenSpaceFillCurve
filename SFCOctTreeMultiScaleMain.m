
% computes space-filling curve for octree data
%% NOTE: VERY important! Matlab uses a different axis direction for 3D coordinates and 2D: y goes into depth in 3D, while y goes up in 2D.

function [clLT, clVisitOrder, fullLT] = SFCOctTreeMultiScaleMain(finestLevelFileName, attr)

close all;
mode = 1; % 1- volume data; 2 - particle
plotOctree = true;% false;
V=[];
testCaseNum = -1;
if nargin >= 1
    %filename = './data/heartCrop.nhdr';
    [filepath, baseFileName, ext] = fileparts(finestLevelFileName);
    if strcmp(ext,'.nrrd') || strcmp(ext, '.nhdr')
        % load scalar volume data
        headerInfo = nhdr_nrrd_read(finestLevelFileName, true);
        % headerInfo = nhdr_nrrd_read('./data/heartCrop.nhdr', true);
        V = headerInfo.data;
        dimX = size(V,2);
        dimY = size(V,1);
        dimZ = size(V,3);
        mode = 1;
    elseif strcmp(ext,'.csv')
        % load particle data
        X = readtable(finestLevelFileName);
        PTS = [X.posX, X.posY, X.posZ];
        % original code:
        PTminmax = [min(PTS,[],1) max(PTS,[],1)];
        PTr = PTminmax(4:6) - PTminmax(1:3);
        %         dimX = 128; dimY = 128; dimZ = 128;
        %% normalize the point locations
        PTS = (PTS - PTminmax(1:3)) ./ PTr;
        numPtsPerBin = 1;
        %         % transform to [0,dim]
        %         plot3(PTS(:,1),PTS(:,2),PTS(:,3), 'o');
        %         Indx = int16(floor(PTS));
        %         tId = transpose(1:size(PTS,1));
        %         for i = 1:length(tId)
        %             V(Indx(i,1),Indx(i,2),Indx(i,3)) = X.pressure(i);%tId(i);%norm(PTS);
        %         end
        %         volshow(V);
        mode = 2;
    end
else
    testCaseNum = 3;
    % use test cases
    numPtsPerBin = 1;
    mode = 1;
    switch testCaseNum
        case 1        % Case 1: simple
            [ V,PTS, dimX, dimY, dimZ] = buildOctTreeTestImage(1);
        case 2        % Case 2: harder
            [ V,PTS, dimX, dimY, dimZ] = buildOctTreeTestImage(2);
        case 3        % Case 3: randomly located points
            [ V,PTS, dimX, dimY, dimZ] = buildOctTreeTestImage(3);
            mode = 2;
            numPtsPerBin = 2;%10;
        case 4  % Case 4: spheres
            % Case 4: randomly placed spheres
            dimX = 128; dimY = 128; dimZ = 128;
            nSphere = 5;
            V = testVolCreate(dimX, dimY, dimZ, nSphere);
%             volshow(V);
            volumeViewer(V);
    end
    
    baseFileName = 'testCase';
end

if nargin < 2
    if mode == 2
        zattr = 'pressure';
    else
        zattr = '';
    end
else
    if mode == 2
        zattr = attr;
    else
        zattr = '';
    end
end

LTfilename = sprintf('LT%s%s.csv', baseFileName, zattr);
VOfilename = sprintf('VO%s%s.csv', baseFileName, zattr);

% generate the octree with #numPtsPerBin pt per bin

corsLvl = 1000;
fineLvl = -1;

if mode == 1
    %% Build volume-based octree
    valDiff = 1;%30;
    % use sphere test case
    maxSize = max(size(V));
    nextPow2 = 2^ceil(log2(maxSize));
    padDimFirstHalf = zeros(3,1);
    padDimSecondHalf = zeros(3,1);
    for i = 1:3
        if mod(size(V,i),2) == 0
            padDimFirstHalf(i) = (nextPow2 - size(V,i))/2;
            padDimSecondHalf(i) = (nextPow2 - size(V,i))/2;
        else
            padDimFirstHalf(i) = nextPow2/2 - floor(size(V,i)/2);
            padDimSecondHalf(i) = nextPow2 - size(V,i) - padDimFirstHalf(i);
        end
    end
    % pad in 3D
    VP = zeros(nextPow2,nextPow2,nextPow2);
    VP(1+padDimFirstHalf(1):padDimFirstHalf(1)+dimY,...
        1+padDimFirstHalf(2):padDimFirstHalf(2)+dimX,...
        1+padDimFirstHalf(3):padDimFirstHalf(3)+dimZ) = V;
    V = VP;
    dimX = nextPow2; dimY = nextPow2; dimZ = nextPow2;
    volumeViewer(VP);

    OT = VolOctree(V, [dimX dimY dimZ], 'binValDiff', valDiff);
    dim = nextPow2;
    % check leaf nodes, i.e., nodes containing data points
    
    minBlockSize = realmax;
    maxBlockSize = -realmax;
    
    binChildren = arrayfun(@(i)find(OT.BinParents==i),1:OT.BinCount,'Un',0)';
    binIsLeaf = cellfun(@isempty, binChildren);
    leafNodeCnt = sum((binIsLeaf));
    %     corsLvl = OT.LeafDepthMin;
    %     fineLvl = OT.DepthMax;
    
    for i = 1:OT.BinCount
        hasChildren = find(OT.BinParents == i,1);
        if isempty(hasChildren)
            bounds = OT.BinBoundaries(i,:);
            if mode == 1
                boundSize = diff(bounds([1:3;4:6])) + [1 1 1];
            else
                boundSize = diff(bounds([1:3;4:6]));
            end
            minBlockSize = min(minBlockSize, boundSize(1));
            maxBlockSize = max(maxBlockSize, boundSize(1));
            
            corsLvl = min(corsLvl, OT.BinDepths(i));
            fineLvl = max(fineLvl, OT.BinDepths(i));
        end
    end
    
    corsBlockSize = dim/ pow2(corsLvl);
    fineBlockSize= dim / pow2(fineLvl);
    blockSizeFactor = maxBlockSize / corsBlockSize;
    disp(fineBlockSize);
    disp(corsBlockSize);
    nlevels = int16(log2(corsBlockSize) - log2(fineBlockSize)+1);
    
    %%
    %% Generate data for each level by aggregation
    allLvls = log2(dim)+1;
    Vtmp = cell(log2(dim)+1,1);
    Vtmp{end} = V;
    for i = allLvls-1:-1:1
        Vtmp{i} = AggregateIma(Vtmp{i+1},2);
    end
    Vlvls = cell(nlevels,1);
    % Caution: V may not be the finest level due to the setting of number of point
    cnt = nlevels;
    for l = 0:fineLvl
        if l >= corsLvl && l <= fineLvl
            Vlvls{cnt} = Vtmp{l+1};%AggregateIma(Vlvls{i-1}, 2);
            corsLvl = min(corsLvl, OT.BinDepths(i));
            fineLvl = max(fineLvl, OT.BinDepths(i));
            cnt = cnt - 1;
        end
    end
else
    %% Build point-based octree
    OT = OcTree(PTS,'binCapacity',numPtsPerBin, 'maxDepth', 8);
    binChildren = arrayfun(@(i)find(OT.BinParents==i),1:OT.BinCount,'Un',0)';
    binIsLeaf = cellfun(@isempty, binChildren);
    leafNodeCnt = sum((binIsLeaf));
    for i = 1:OT.BinCount
        if binIsLeaf(i)
            corsLvl = min(corsLvl, OT.BinDepths(i));
            fineLvl = max(fineLvl, OT.BinDepths(i));
        end
    end
    dim = pow2(fineLvl);
    V = zeros(dim,dim,dim);
    CntV = zeros(dim,dim,dim);
    cnt = 0;
    if testCaseNum == -1
        for i = 1:OT.BinCount
            if binIsLeaf(i)
    %             doplot3(pts(OT.PointBins==i,:),'.','Color',cols(i,:))
                val = X.pressure(OT.PointBins==i,:);
                if isempty(val) 
                    continue;
                end
                if length(val) > 1
                    disp(length(val));
                    cnt = cnt + 1;
                end
                binMinMax = OT.BinBoundaries(i,:);
                binMinMax = binMinMax;
                V(binMinMax(1)+1:binMinMax(4),binMinMax(2)+1:binMinMax(5),binMinMax(3)+1:binMinMax(6)) = mean(val);
            end
        end
    else
        for i = 1:OT.BinCount
            if binIsLeaf(i)
                %  

                val = double(i);%PTS(OT.PointBins==i,:);
                if isempty(val)
                    continue;
                end
                if length(val) > 1
                    disp(length(val));
                    cnt = cnt + 1;
                end
                binMinMax = OT.BinBoundaries(i,:);
                binMinMax = binMinMax .* dim;
                V(binMinMax(1)+1:binMinMax(4),binMinMax(2)+1:binMinMax(5),binMinMax(3)+1:binMinMax(6)) = mean(val);
            end
        end
    end
    disp(cnt);
    %% PTS should be normalized
%     for i = 1:length(PTS)
%         idx = int16(PTS(i,:) .* (dim-1)) + 1;
%         CntV(idx(1),idx(2),idx(3)) = CntV(idx(1),idx(2),idx(3)) + 1;
% %         V(idx(1),idx(2),idx(3)) = V(idx(1),idx(2),idx(3)) + (X.pressure(i) - V(idx(1),idx(2),idx(3)))/CntV(idx(1),idx(2),idx(3)); % compute the running mean
%         if(strcmp(zattr,'pressure') == 1)
%             V(idx(1),idx(2),idx(3)) = V(idx(1),idx(2),idx(3)) + (X.pressure(i) - V(idx(1),idx(2),idx(3)))/CntV(idx(1),idx(2),idx(3)); % compute the running mean
%         elseif(strcmp(zattr, 'density')==1) % subtract 1000
%             V(idx(1),idx(2),idx(3)) = V(idx(1),idx(2),idx(3)) + (X.density(i) - 1000 - V(idx(1),idx(2),idx(3)))/CntV(idx(1),idx(2),idx(3)); % compute the running mean
%         elseif(strcmp(zattr, 'velocityX')==1)
%             V(idx(1),idx(2),idx(3)) = V(idx(1),idx(2),idx(3)) + (X.velocityX(i) -...
%                 V(idx(1),idx(2),idx(3)))/CntV(idx(1),idx(2),idx(3)); % compute the running mean
%         elseif(strcmp(zattr, 'velocityY')==1)
%             V(idx(1),idx(2),idx(3)) = V(idx(1),idx(2),idx(3)) + (X.velocityY(i) -...
%                 V(idx(1),idx(2),idx(3)))/CntV(idx(1),idx(2),idx(3)); % compute the running mean
%         elseif(strcmp(zattr, 'velocityZ')==1)
%             V(idx(1),idx(2),idx(3)) = V(idx(1),idx(2),idx(3)) + (X.velocityZ(i) -...
%                 V(idx(1),idx(2),idx(3)))/CntV(idx(1),idx(2),idx(3)); % compute the running mean
%         elseif(strcmp(zattr, 'speed')==1)
%             V(idx(1),idx(2),idx(3)) = V(idx(1),idx(2),idx(3)) + (norm([X.velocityX(i),X.velocityY(i),X.velocityZ(i)]) -...
%                 V(idx(1),idx(2),idx(3)))/CntV(idx(1),idx(2),idx(3)); % compute the running mean
%         end
% 
%     end
%     maxV = max(max(max(V)));
%     disp(maxV);
%     meanCnt = mean(mean(CntV));
%     disp(meanCnt);
%     maxCnt = max(max(max(CntV)));
%     disp(maxCnt);
%     volumeViewer(V);
    nlevels = fineLvl - corsLvl + 1;
    Vlvls = cell(nlevels,1);
    % Caution: V may not be the finest level due to the setting of number of point
    Vlvls{1} = V;
    Fct = 5;
    for l = 2:nlevels
        Vlvls{l} = AggregateIma(Vlvls{l-1},2,Fct);
        disp(max(max(max(Vlvls{l}))));
    end
    
end
%% plot the octree

if plotOctree
    figure;
    boxH = OT.plot;
    cols = lines(OT.BinCount);
    doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
    for i = 1:OT.BinCount
        if OT.BinDepths(i) == 3
            set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
            %        doplot3(PTS(OT.PointBins==i,:),'.','Color',cols(i,:))
        end
    end
    axis image, view(3)
    hold on;
end



%% Calculate the SFC
if mode == 1
    [clLT, clVisitOrder] = SFCOctTree(Vlvls, OT, blockSizeFactor, minBlockSize(1)/blockSizeFactor);
else
    [clLT, clVisitOrder] = SFCOctTreePoint(Vlvls, OT, 1, [0 0 0]);
end
%        [clLT, clVisitOrder] = SFCOctTree(Vlvls, OT, blockSizeFactor, minBlockSize(1));
%% visualization
% plot3(clVisitOrder(:,1),clVisitOrder(:,2),clVisitOrder(:,3),'r-');
idOrder = 1:length(clLT);
lineColorCoded(clVisitOrder(:,1)', clVisitOrder(:,2)', clVisitOrder(:,3)', idOrder);
hold off;

figure;
% plot3(clVisitOrder(:,1),clVisitOrder(:,3),clVisitOrder(:,2),'b-');
lineColorCoded(clVisitOrder(:,1)', clVisitOrder(:,2)', clVisitOrder(:,3)', idOrder);

figure, hold on;
subplot(2,1,1), plot(1:length(clLT),clLT);title('Multiscale SFC');
fprintf('total octree (leaf) nodes = %d, total found nodes = %d\n', leafNodeCnt, length(clLT));

% reconstruct the full sfc
fullLT = fullSFC3D(clLT);
subplot(2,1,2), plot(1:length(fullLT),fullLT);title('Reconstructed SFC');
hold off;

% draw the traversal order


%     orderBlocks = repmat(uint64(0),);
%     dim = nextPow2;
%     totalNodes = 0;
%     while(dim>=1)
%         [vals,r,c] = qtgetblk(V,S,dim);
%
%         if ~isempty(vals)
%             totalNodes = totalNodes + length(r);
%             for i = 1:length(r)
%                 v = [r(i),c(i)];
%                 orderId1 = find(clVisitOrder(:,1) == v(1));
%                 orderId2 = find(clVisitOrder(:,2) == v(2));
%                 orderId = intersect(orderId1, orderId2);
%                 if ~isempty(orderId)
%                     orderBlocks(r(i):r(i)+dim-1,c(i):c(i)+dim-1,:) = orderId;
% %                 else
% %                     errStr = sprintf('missed block [%d %d]', v(1),v(2));
% %                     disp(errStr);
%                 end
%             end
%         end
%         dim = dim / 2;
%
%     end
%     figure;
%     fprintf('total quadtree nodes = %d, total found nodes = %d\n', totalNodes, length(clLT));

% writeout results

csvwrite(LTfilename, clLT);
csvwrite(VOfilename, clVisitOrder);
end
