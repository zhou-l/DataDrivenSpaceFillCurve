%coarse level image file name: clFilename
%next level image file name: nlFilename
function [clLT, clVisitOrder, fullLT] = vis20Fig8Gen
    close all;
    
    % use test cases
    dimX = 32;%64;
    dimY = 32;%64;
    dimZ = 1;
    nSphere = 5;
%           V = testVolCreate(dimX, dimY, dimZ, nSphere);
    [V, dimX, dimY] = buildQuadTreeTestImage();
    baseFileName = 'testCase';
    figure;fV = flip(V,1); imagesc(fV);

    LTfilename = sprintf('LT%s.csv', baseFileName);
    VOfilename = sprintf('VO%s.csv', baseFileName);
    global zalpha;
    zalpha = 0.1;

% use sphere test case
    maxSize = max(size(V));
    nextPow2 = 2^ceil(log2(maxSize));
    padDimFirstHalf = zeros(2,1);
    padDimSecondHalf = zeros(2,1);
    for i = 1:2
        if mod(size(V,i),2) == 0
            padDimFirstHalf(i) = (nextPow2 - size(V,i))/2;
            padDimSecondHalf(i) = (nextPow2 - size(V,i))/2;
        else
            padDimFirstHalf(i) = nextPow2/2 - floor(size(V,i)/2);
            padDimSecondHalf(i) = nextPow2 - size(V,i) - padDimFirstHalf(i);
        end
    end
    V = padarray(V, [padDimFirstHalf(1),padDimFirstHalf(2)], 'pre');
    V = padarray(V, [padDimSecondHalf(1),padDimSecondHalf(2)],'post');

    thres = 0.08 * (max(max(V)) - min(min(V)));
%     thres = 0; % no threshold
    % generate quadtree
    S = qtdecomp(V, thres);
%      S = qtdecomp(V);
    % detect the finest and coarsest levels
    fineBlockDim = nextPow2;
    corsBlockDim = 1;
    dim = nextPow2;
%     for dim = nextPow2:nextPow2/2:1
    blocks = repmat(uint8(0),size(S));
    while(dim>=1)
        numblocks = length(find(S==dim));
        if(numblocks > 0)
            % find smallest block size
            if corsBlockDim < dim
                corsBlockDim = dim;
            end
            if fineBlockDim > dim
                fineBlockDim = dim;
            end
             values = repmat(uint8(1),[dim dim numblocks]);
             values(2:dim,2:dim,:) = 0;
             blocks = qtsetblk(blocks,S,dim,values);
        end
        dim = dim / 2;
    end
    disp(fineBlockDim);
    disp(corsBlockDim);
    figure, imshow(blocks, []);
    nlevels = log2(corsBlockDim) - log2(fineBlockDim)+1;
    % generate data for each level by aggregation
    Vlvls = cell(nlevels,1);
    Vlvls{1} = V;
    for i = 2:nlevels
%         Vlvls{i} = imresize(Vlvls{i-1}, 0.5, 'nearest');
        Vlvls{i} = AggregateIma(Vlvls{i-1}, 2);
%          Vlvls{i} = impyramid(Vlvls{i-1}, 'reduce');
    end
    % Calculate the SFC
   
    [clLT, clVisitOrder] = SFCQuadTree(Vlvls, S, [dimY, 1]);
%      clLT = allNodesLinearFunc(clVisitOrder, cI);% visiting all nodes without aggregation

    figure, hold on;
    subplot(2,1,1), plot(1:length(clLT),clLT);title('Multiscale SFC');

    
    % reconstruct the full sfc
    fullLT = fullSFC(clLT);
    subplot(2,1,2), plot(1:length(fullLT),fullLT);title('Reconstructed SFC');
    hold off;
    
    % draw the traversal order
    orderBlocks = repmat(uint64(0),size(S));
    dim = nextPow2;
    totalNodes = 0;
    
    ncVisitOrder = clVisitOrder;
    while(dim>=1)
        [vals,r,c] = qtgetblk(V,S,dim);
        
        if ~isempty(vals)
            totalNodes = totalNodes + length(r);
            for i = 1:length(r)
                v = [r(i),c(i)];
                orderId1 = find(clVisitOrder(:,1) == v(1));
                orderId2 = find(clVisitOrder(:,2) == v(2));
                orderId = intersect(orderId1, orderId2);
                if ~isempty(orderId)
                    blkCtr = [r(i)+dim/2-0.5,c(i)+dim/2-0.5];
                    orderBlocks(r(i):r(i)+dim-1,c(i):c(i)+dim-1,:) = orderId;
                    ncVisitOrder(orderId,:) = blkCtr;
%                 else
%                     errStr = sprintf('missed block [%d %d]', v(1),v(2));
%                     disp(errStr);
                end
            end
        end
        dim = dim / 2;
        
    end
    figure; hold on;
    fprintf('total quadtree nodes = %d, total found nodes = %d\n', totalNodes, length(clLT));
%     blender = vision.AlphaBlender;  
    imagesc(V);
    % compute block center
    line(ncVisitOrder(:,2), ncVisitOrder(:,1), 'Color', 'blue', 'LineWidth', 2);
%    line(clVisitOrder(:,2), clVisitOrder(:,1), 'Color', 'white', 'LineWidth', 2);
    title('Traverse order');
    hold off;
    % writeout results

    csvwrite(LTfilename, clLT);
    csvwrite(VOfilename, clVisitOrder);
end