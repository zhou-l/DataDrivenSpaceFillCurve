
function [contextLT, contextVO, myLT, myVO, HilLT, HilVO, lineLT, lineVO] = SFCMethodsTestMain(V, isPadded, isPlot)
if nargin < 2
    zisPad = true;
else
    zisPad = isPadded;
end

if nargin < 3
    zisPlot = true;
else
    zisPlot = isPlot;
end
%     testHamilton;
dimX = size(V,2);
dimY = size(V,1);
dimZ = size(V,3);
maxLags = 100;
refPt = [1 1]; %[dimY/2, dimX/2];
close all;

if zisPlot
    figure, imagesc(V);
end
global useLocWeight;
useLocWeight = false;
global zalpha;
zalpha = 0.1; %0.9;
%% Original Context-based
[WpXR, WpXL, WpYU, WpYD] = buildSmallCircsDualGraph(V);
% figure, imagesc(WpXR);
% figure, imagesc(WpXL);
% figure, imagesc(WpYD);
% figure, imagesc(WpYU);

[T,mstSet] = findMinSpanTree(WpXR, WpXL, WpYU, WpYD);
% for reference


figure;
hold on;

%   [LT,visitOrder] = linearizeHamCycle(T, V, [dimY,1]);
[LT, visitOrder] = linearizeHamCycleMerge(T, V, [dimY, 1]);
contextLT = LT;
contextVO = visitOrder;
avgAutoCorrContext = compAvgAutoCorr(LT, maxLags);
% refOrder = 1:length(LT);
% refOrder = transpose(refOrder);
refOrder = zeros(length(LT),1);
dist = vecnorm(visitOrder(:,:) - refPt, 2, 2);
avgAutoCorrDistContext = compAvgAutoCorr(dist - refOrder, maxLags);

dist2 = vecnorm(visitOrder(:,:) - [dimX/2,dimY/2], 2, 2);


subplot(1,4, 2);hold on;
plot(1:length(LT),LT(:,1));  %plot(1:length(LT),LT(:,2));
hold off;
title('Context SFC Hamilton');
csfcHamOrder = visitOrder;
%% Ours
useLocWeight = true;
zalpha = 0.1;
[Gmag,Gdir] = imgradient(V);
[GpXR, GpXL, GpYU, GpYD] = buildSmallCircsDualGraph(Gmag);

% [GpXR, GpXL, GpYU, GpYD] = buildSmallCircsDualGraph(Gmag);
[T,mstSet] = findMinSpanTree2(WpXR, WpXL, WpYU, WpYD, GpXR, GpXL, GpYU, GpYD);
%     [T,mstSet] = findMinSpanTree(WpXR, WpXL, WpYU, WpYD);
%     [LT,visitOrder] = linearizeHamCycle(T, V, [dimY,1]);
[LT, visitOrder] = linearizeHamCycleMerge(T, V, [dimY, 1]);
myLT = LT;
myVO = visitOrder;

subplot(1,4, 1); hold on;
plot(1:length(LT),LT(:,1));
%     plot(1:length(LT),LT(:,2));
hold off;
title('Locality Hamilton');

%     [LT,visitOrder] = linearizeWithSpaceFillCurveDST(T,V);
%     subplot(3,2,4), plot(1:length(LT),LT);title('DFS Context SFC');
%     dfsOrder = visitOrder;
avgAutoCorrOurs = compAvgAutoCorr(LT, maxLags);
dist = vecnorm(visitOrder(:,:) - refPt, 2, 2);
avgAutoCorrDistOurs = compAvgAutoCorr(dist - refOrder, maxLags);
%% Scanline
[LT,visitOrder] = linearizeScanline(V);
lineLT = LT;
lineVO = visitOrder;
subplot(1,4,4);plot(1:length(LT),LT);title('Scanline');

avgAutoCorrScanline = compAvgAutoCorr(LT, maxLags);
dist = vecnorm(visitOrder(:,:) - refPt, 2, 2);
avgAutoCorrDistScanline = compAvgAutoCorr(dist - refOrder, maxLags);
%% Hilbert
HilLT = [];
HilVO = [];
avgAutoCorrHilbert = [];
avgAutoCorrDistHilbert = [];
if zisPad
    [rD, rI] = hilbertCurve2D(V);
    subplot(1,4,3);plot(1:length(rD),rD);title('Hilbert');
    HilLT = rD;
    avgAutoCorrHilbert = compAvgAutoCorr(HilLT, maxLags);
    % get 2D coords from rI
    HilVO = zeros(length(rI),2);
    for i = 1:length(rI)
        HilVO(i,1) = floor((rI(i)-1)/dimX) + 1;
        HilVO(i,2) = mod(rI(i) - 1, dimX)+1;
    end
    dist = vecnorm(HilVO(:,:) - refPt, 2, 2);
    avgAutoCorrDistHilbert = compAvgAutoCorr(dist - refOrder, maxLags);
end
hold off;

% disp(visitOrder)
travOrder = zeros(dimY,dimX);
travOrder2 = zeros(dimY,dimX);
travOrderHil = zeros(dimY, dimX);
travOrderLine = zeros(dimY, dimX);
for i = 1:length(myVO)
    travOrder(myVO(i,1),myVO(i,2)) = i;
    travOrder2(contextVO(i,1),contextVO(i,2)) = i;
%    travOrderHil(HilVO(i,1),HilVO(i,2)) = i;
    travOrderLine(lineVO(i,1),lineVO(i,2)) = i;
    %         if showAnim
    %             imagesc(travOrder);
    %             pause(0.04);
    %         end
end

figure, subplot(1,4,1); 
imagesc(travOrder); line(myVO(:,2), myVO(:,1)); title('Ours');
subplot(1,4,2); 
imagesc(travOrder2); line(contextVO(:,2),contextVO(:,1));title('Context-based');
subplot(1,4,3); 
imagesc(travOrderHil);
line(HilVO(:,2),HilVO(:,1));title('Hilbert');
subplot(1,4,4);
imagesc(travOrderLine);
line(lineVO(:,2),lineVO(:,1));title('Scanline');
% plot avg corrs
figure; hold on;
plot(1:maxLags, avgAutoCorrScanline, 'k-', 'LineWidth', 5);
plot(1:maxLags, avgAutoCorrContext, 'g-', 'LineWidth', 5);
plot(1:maxLags, avgAutoCorrOurs, 'b-', 'LineWidth', 5);
if~isempty(avgAutoCorrHilbert)
    plot(1:maxLags, avgAutoCorrHilbert, 'r-', 'LineWidth', 5);
end
legend('Scanline', 'Context', 'Ours', 'Hilbert');
title('AutoCorr values');
hold off;

% plot avg corrs
figure; hold on;
plot(1:maxLags, avgAutoCorrDistScanline, 'k-');
plot(1:maxLags, avgAutoCorrDistContext, 'g-');
plot(1:maxLags, avgAutoCorrDistOurs, 'b-');
if~isempty(avgAutoCorrDistHilbert)
    plot(1:maxLags, avgAutoCorrDistHilbert, 'r-');
end
legend('Scanline', 'Context', 'Ours', 'Hilbert');
title('AutoCorr Eucledian Dist');
hold off;

end

function [LT, visitOrder] = linearizeWithSpaceFillCurveMSTorder(mstSet, I)
seqLen= length(mstSet);
LT = zeros(seqLen,1); % the linearized data
visitOrder = zeros(seqLen,2); %coordinate visiting order
for i = 1:seqLen
    currNode = mstSet(i,:);
    visitOrder(i,:) = currNode;
    xx = 2*(currNode(2)-1)+1; % index of the original grid
    yy = 2*(currNode(1)-1)+1;
    LT(i,:) = 0.25 * (I(yy,xx)+I(yy+1,xx)+I(yy,xx+1)+I(yy+1,xx+1));
end
end


%Depth first search
function dfs(zmstT, node)
global order;
global V;
global LT;
global visitOrder;
stack = [];
stack(end+1,:) = node;

numChildren = 0;

while ~isempty(stack)
    cnode = stack(1,:);
    xx = 2*(cnode(2)-1)+1; % index of the original grid
    yy = 2*(cnode(1)-1)+1;
    LT(order,:) = 0.25 * (V(yy,xx)+V(yy+1,xx)+V(yy,xx+1)+V(yy+1,xx+1));
    visitOrder(order,:) = cnode;
    order = order + 1;
    % remove first
    stack(1,:) = [];
    children = zmstT(cnode(1),cnode(2),:);
    for i = 1:4
        if children(i) > 0
            numChildren = numChildren + 1;
            %                switch(children(i)) % this is for older priority list
            %                version
            switch(i) % for the current channel==direction version
                case 1%up
                    nextNode = cnode+[-1,0];
                    
                case 2%down
                    nextNode = cnode+[1,0];
                    
                case 3%right
                    nextNode = cnode+[0,1];
                    
                case 4%left
                    nextNode= cnode+[0,-1];
            end
            
            %            dfs(zmstT, nextNode);
            stack = [nextNode;stack];
        end
    end
end
end

% DFS main function
function [LT, visitOrder] = linearizeWithSpaceFillCurveDST(mstT, I)
seqLen= size(mstT,1)*size(mstT,2);
global LT;
global visitOrder;
visitOrder = zeros(seqLen,2); %coordinate visiting order
LT = zeros(seqLen,1); % the linearized data
currNode = [size(I,1)/2,1];
global order;
order = 1;
dfs(mstT, currNode);
end


%Breadth First Search ordering for space filling curve generation
function [LT,visitOrder] = linearizeWithSpaceFillCurveBST(mstT, I)
seqLen= size(mstT,1)*size(mstT,2);
LT = zeros(seqLen,1); % the linearized data
% do depth first traversal of the tree
currNode = [size(I,1)/2,1];
visitOrder = zeros(seqLen,2); %coordinate visiting order

order = 1;
nodeStack = [];
visitOrder(1,:) = currNode;
nodeStack(end+1,:) = currNode;

while ~isempty(nodeStack)
    currNode = nodeStack(1,:);
    visitOrder(order,:) = currNode;
    % compute average value for the linearized function
    xx = 2*(currNode(2)-1)+1; % index of the original grid
    yy = 2*(currNode(1)-1)+1;
    LT(order,:) = 0.25 * (I(yy,xx)+I(yy+1,xx)+I(yy,xx+1)+I(yy+1,xx+1));
    nodeStack(1,:) = [];%remove the first element
    children = mstT(currNode(1),currNode(2),:);
    for i = 1:4 %check 4 children
        nn = children(1,1,i);
        switch(nn)
            case 1%up
                nodeStack(end+1,:) = currNode+[-1,0];
            case 2%down
                nodeStack(end+1,:) = currNode+[1,0];
            case 3%right
                nodeStack(end+1,:) = currNode+[0,1];
            case 4%left
                nodeStack(end+1,:) = currNode+[0,-1];
            otherwise
                break;
        end
    end
    order = order+1;
end
disp(order);
end

% for reference: scan line order
function [LT, visitOrder] = linearizeScanline(I)
blockDimY = size(I,1);
blockDimX = size(I,2);
seqLen= blockDimY * blockDimX;
visitOrder = zeros(seqLen,2); %coordinate visiting order
LT = zeros(seqLen,1); % the linearized data
order = 1;
for y = blockDimY:-1:1
    for x = 1:blockDimX
        visitOrder(order,:) = [y,x];
        LT(order,:) = I(y,x);
        order = order + 1;
    end
end
end

function [LT, visitOrder] = linearizeScanlineBlk(I)
blockDimY = size(I,1)/2;
blockDimX = size(I,2)/2;
seqLen= blockDimY * blockDimX;
visitOrder = zeros(seqLen,2); %coordinate visiting order
LT = zeros(seqLen,1); % the linearized data
currNode = [size(I,1)/2,1];
order = 1;
for y = blockDimY:-1:1
    for x = 1:blockDimX
        xx = 2*(x-1)+1; % index of the original grid
        yy = 2*(y-1)+1;
        
        visitOrder(order,:) = [y,x];
        LT(order,:) = 0.25 * (I(yy,xx)+I(yy+1,xx)+I(yy,xx+1)+I(yy+1,xx+1));
        order = order + 1;
    end
end
end