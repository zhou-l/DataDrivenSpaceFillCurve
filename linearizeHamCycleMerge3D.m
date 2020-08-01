%Traversing the image I using a Hamilton cycle generated from a precomputed minimum spanning tree
% use tree as input
function [zLT,zVisitOrder] = linearizeHamCycleMerge3D(mstT, I, entryPix)
    currPix = entryPix;
    currCirc = [floor((currPix(1)+1)/2),floor((currPix(2)+1)/2),floor((currPix(3)+1)/2)];
  
    global nneighbors;
    global neighborDirs;
    LT = zeros(size(I,1)*size(I,2)*size(I,3),1);
    hamCyc = zeros(size(I,1)*size(I,2)*size(I,3),3);


    nneighbors = 6;
    neighborDirs = [
        -1,0,0; % up
        1,0,0; % down
        0,1,0; %right
        0,-1,0;%left
        0,0,1;%back
        0,0,-1%front
        ];
%     dirSeq = genDirSeq(mstSet);
    toDir = [-1,0,0];
    stack = currCirc; 
    dirStack = toDir;
    cnt = 0;
    while ~isempty(stack)
        cnt = cnt + 1;
        currCirc = stack(1,1:3);   
        toDir = dirStack(1,1:3);
        % remove first
        stack(1,:) = [];
        dirStack(1,:) = [];
        
        children = mstT(currCirc(1),currCirc(2),currCirc(3),:);
        children = reshape(children, [nneighbors,1]);
        
        hamCyc = mergeToCyc3D(hamCyc, currCirc, toDir);
 
        for i = 1:nneighbors
            if children(i) > 0
                   % for the current channel==direction version
                   tryDir = neighborDirs(i,:);
                   nextCirc = currCirc + tryDir;        
                   toDir = tryDir;
                  % if next pixel is not in circuit, we visit
                  % current pixel and skip all
                   stack = [nextCirc;stack];
                   dirStack = [toDir;dirStack];
            end
        end
    end     

 % traverse using the hamCycle
    zLT = zeros(length(hamCyc),1);
    for i = 1:length(hamCyc)
        zLT(i,:) = I(hamCyc(i,1),hamCyc(i,2),hamCyc(i,3),:);
    end
    zVisitOrder = hamCyc;
    
end

% get element ham cube given the indexing of cube vertices
function  circGraphAdjMat = getElementHamCube(hamMode)
    circGraphAdjMat = zeros(8,8);
    switch hamMode
        case 1
               %   p6   p7
               %  /|   /|
               % p1   p2|
               % | p5-|p8
               % p4---p3
            circGraphAdjMat =...
                [0 0 0 1 0 1 0 0; %p1<->p4, p1<->p6
                 0 0 1 0 0 0 1 0; %p2<->p3, p2<->p7
                 0 1 0 1 0 0 0 0; %p3<->p2, p3<->p4
                 1 0 1 0 0 0 0 0; %p4<->p1, p4<->p3
                 0 0 0 0 0 1 0 1; %p5<->p6, p5<->p8
                 1 0 0 0 1 0 0 0; %p6<->p1, p6<->p5
                 0 1 0 0 0 0 0 1; %p7<->p2, p7<->p8
                 0 0 0 0 1 0 1 0];
        case 2
               %   p6--p7
               %   |    |
               % p1---p2|
               % | /p5| p8
               % p4  p3/
              circGraphAdjMat =...
                [0 1 0 1 0 0 0 0; %p1<->p4, p1<->p2
                 1 0 1 0 0 0 0 0; %p2<->p3, p2<->p1
                 0 1 0 0 0 0 0 1; %p3<->p2, p3<->p8
                 1 0 0 0 1 0 0 0; %p4<->p1, p4<->p5
                 0 0 0 1 0 1 0 0; %p5<->p6, p5<->p4
                 0 0 0 0 1 0 1 0; %p6<->p7, p6<->p5
                 0 0 0 0 0 1 0 1; %p7<->p6, p7<->p8
                 0 0 1 0 0 0 1 0];
        case 3
               %   p6--- p7
               %  /      |
               % p1---p2 |
               %  /p5-|- p8
               % p4--p3/
               circGraphAdjMat =...
                [0 1 0 0 0 1 0 0; %p1<->p2, p1<->p6
                 1 0 1 0 0 0 0 0; %p2<->p3, p2<->p1
                 0 1 0 1 0 0 0 0; %p3<->p2, p3<->p4
                 0 0 1 0 1 0 0 0; %p4<->p3, p4<->p5
                 0 0 0 1 0 0 0 1; %p5<->p8, p5<->p4
                 1 0 0 0 0 0 1 0; %p6<->p7, p6<->p1
                 0 0 0 0 0 1 0 1; %p7<->p6, p7<->p8
                 0 0 0 0 1 0 1 0]; %p8<->p5, p8<->p7
         case 4
               %   p6--- p7
               %    |   / 
               % p1---p2 
               %  | p5-- p8
               % p4--p3/
               circGraphAdjMat =...
                [0 1 0 1 0 0 0 0; %p1<->p2, p1<->p4
                 1 0 0 0 0 0 1 0; %p2<->p3, p2<->p7
                 0 0 0 1 0 0 0 1; %p3<->p4, p3<->p8
                 1 0 1 0 0 0 0 0; %p4<->p3, p4<->p1
                 0 0 0 0 0 1 0 1; %p5<->p8, p5<->p6
                 0 0 0 0 1 0 1 0; %p6<->p7, p6<->p5
                 0 1 0 0 0 1 0 0; %p7<->p6, p7<->p2
                 0 0 1 0 1 0 0 0]; %p8<->p5, p8<->p3
        case 5 
               %   p6--- p7
               %   /    / 
               % p1    p2 
               %  | p5-|- p8
               % p4/   p3/
               circGraphAdjMat =...
                [0 0 0 1 0 1 0 0; %p1<->p6, p1<->p4
                 0 0 1 0 0 0 1 0; %p2<->p3, p2<->p7
                 0 1 0 0 0 0 0 1; %p3<->p2, p3<->p8
                 1 0 0 0 1 0 0 0; %p4<->p5, p4<->p1
                 0 0 0 1 0 0 0 1; %p5<->p8, p5<->p4
                 1 0 0 0 0 0 1 0; %p6<->p7, p6<->p1
                 0 1 0 0 0 1 0 0; %p7<->p6, p7<->p2
                 0 0 1 0 1 0 0 0]; %p8<->p5, p8<->p3
         case 6 
               %   p6    p7
               %   /|    /| 
               % p1-|---p2|
               %   p5    p8
               % p4/---p3/
               circGraphAdjMat =...
                [0 1 0 0 0 1 0 0; %p1<->p6, p1<->p2
                 1 0 0 0 0 0 1 0; %p2<->p1, p2<->p7
                 0 0 0 1 0 0 0 1; %p3<->p4, p3<->p8
                 0 0 1 0 1 0 0 0; %p4<->p5, p4<->p3
                 0 0 0 1 0 1 0 0; %p5<->p6, p5<->p4
                 1 0 0 0 1 0 0 0; %p6<->p5, p6<->p1
                 0 1 0 0 0 0 0 1; %p7<->p8, p7<->p2
                 0 0 1 0 0 0 1 0]; %p8<->p7, p8<->p3
    end
end
% helper functions
function hamCycNew = mergeToCyc3D(hamCycOld, currCirc, mergeDir)
   % use the association rule in paper "3D Hardware Canaries"
   hamCycNew = hamCycOld;
   % number of existing pixels in the hampath
   lastN = find(~any(hamCycOld,3),1)-1;
 
   % Six cases of Ham cubes
   % vertex encoding
   %   p6---p7
   %  /|   /|
   % p1---p2|
   % | p5 -|p8
   % p4---p3
   % get pixels in the circuit to be merged
   
   pCC = zeros(8,3); 
   pCC(1,:) = [2*(currCirc(1)-1)+1, 2*(currCirc(2)-1)+1, 2*(currCirc(3)-1)+1]; % top left front
   pCC(2,:) = pCC(1,:) + [0,1,0]; % top right front
   pCC(3,:) = pCC(2,:) + [1,0,0]; % bottom right front
   pCC(4,:) = pCC(3,:) + [0,-1,0];% bottom left front
   
   pCC(5,:) = pCC(4,:) + [0,0,1]; % bottom left back
   pCC(6,:) = pCC(5,:) + [-1,0,0];% top left back
   pCC(7,:) = pCC(6,:) + [0,1,0]; % top right back
   pCC(8,:) = pCC(7,:) + [1,0,0]; % bottom right back
   
   pCCMin = pCC(1,:);
   pCCMax = pCC(8,:);
   
 % no pixels in the cycle yet!  
  if lastN == 0 
      edgeBreakId = 0;
      tDir = mergeDir;
      % first node in hamCyc
        if tDir(1) ~= 0
            % merge in y direction
            if tDir(1) == 1 % going down % Case 1
                    hamCycNew(edgeBreakId+1,:) = pCC(1,:); % top left front
                    hamCycNew(edgeBreakId+4,:) = pCC(2,:); % top right front
                    hamCycNew(edgeBreakId+3,:) = pCC(3,:); % bottom right front
                    hamCycNew(edgeBreakId+2,:) = pCC(4,:); % bottom left front
                    hamCycNew(edgeBreakId+7,:) = pCC(5,:); % bottom left back
                    hamCycNew(edgeBreakId+8,:) = pCC(6,:); % top left back
                    hamCycNew(edgeBreakId+5,:) = pCC(7,:); % top right back
                    hamCycNew(edgeBreakId+6,:) = pCC(8,:); % bottom right back
            elseif tDir(1) == -1 % going up % Case 2
                    hamCycNew(edgeBreakId+2,:) = pCC(1,:); % top left front
                    hamCycNew(edgeBreakId+3,:) = pCC(2,:); % top right front
                    hamCycNew(edgeBreakId+4,:) = pCC(3,:); % bottom right front
                    hamCycNew(edgeBreakId+1,:) = pCC(4,:); % bottom left front
                    hamCycNew(edgeBreakId+8,:) = pCC(5,:); % bottom left back
                    hamCycNew(edgeBreakId+7,:) = pCC(6,:); % top left back
                    hamCycNew(edgeBreakId+6,:) = pCC(7,:); % top right back
                    hamCycNew(edgeBreakId+5,:) = pCC(8,:); % bottom right back
            end
        elseif tDir(2)~= 0
            % merge in x direction
            if tDir(2) == 1 % going right
                    hamCycNew(edgeBreakId+4,:) = pCC(1,:); % top left front
                    hamCycNew(edgeBreakId+3,:) = pCC(2,:); % top right front
                    hamCycNew(edgeBreakId+2,:) = pCC(3,:); % bottom right front
                    hamCycNew(edgeBreakId+1,:) = pCC(4,:); % bottom left front
                    hamCycNew(edgeBreakId+8,:) = pCC(5,:); % bottom left back
                    hamCycNew(edgeBreakId+5,:) = pCC(6,:); % top left back
                    hamCycNew(edgeBreakId+6,:) = pCC(7,:); % top right back
                    hamCycNew(edgeBreakId+7,:) = pCC(8,:); % bottom right back
            elseif tDir(2) == -1 % going left
                    hamCycNew(edgeBreakId+8,:) = pCC(1,:); % top left front
                    hamCycNew(edgeBreakId+7,:) = pCC(2,:); % top right front
                    hamCycNew(edgeBreakId+2,:) = pCC(3,:); % bottom right front
                    hamCycNew(edgeBreakId+1,:) = pCC(4,:); % bottom left front
                    hamCycNew(edgeBreakId+4,:) = pCC(5,:); % bottom left back
                    hamCycNew(edgeBreakId+5,:) = pCC(6,:); % top left back
                    hamCycNew(edgeBreakId+6,:) = pCC(7,:); % top right back
                    hamCycNew(edgeBreakId+3,:) = pCC(8,:); % bottom right back
            end
        else
            % merge in z direction
           if tDir(3) == -1 % going front
                    hamCycNew(edgeBreakId+8,:) = pCC(1,:); % top left front
                    hamCycNew(edgeBreakId+5,:) = pCC(2,:); % top right front
                    hamCycNew(edgeBreakId+4,:) = pCC(3,:); % bottom right front
                    hamCycNew(edgeBreakId+1,:) = pCC(4,:); % bottom left front
                    hamCycNew(edgeBreakId+2,:) = pCC(5,:); % bottom left back
                    hamCycNew(edgeBreakId+7,:) = pCC(6,:); % top left back
                    hamCycNew(edgeBreakId+6,:) = pCC(7,:); % top right back
                    hamCycNew(edgeBreakId+3,:) = pCC(8,:); % bottom right back
            elseif tDir(3) == 1 % going back
                    hamCycNew(edgeBreakId+6,:) = pCC(1,:); % top left front
                    hamCycNew(edgeBreakId+5,:) = pCC(2,:); % top right front
                    hamCycNew(edgeBreakId+2,:) = pCC(3,:); % bottom right front
                    hamCycNew(edgeBreakId+1,:) = pCC(4,:); % bottom left front
                    hamCycNew(edgeBreakId+8,:) = pCC(5,:); % bottom left back
                    hamCycNew(edgeBreakId+7,:) = pCC(6,:); % top left back
                    hamCycNew(edgeBreakId+4,:) = pCC(7,:); % top right back
                    hamCycNew(edgeBreakId+3,:) = pCC(8,:); % bottom right back
            end
        end
        return;
  end
  
   
   % find adjacent pixels in hamCyc and currCirc
   adjPixInHamCyc = zeros(4,4); % first 3 channels are coordinate, the last element is index in ham cycle
   adjPixInCurrCirc = zeros(4,4);
   cnt = 1;
   allPixFound = false;
   for i = 1:lastN
       circInHamCyc = [floor((hamCycOld(i,1)+1)/2),floor((hamCycOld(i,2)+1)/2),floor((hamCycOld(i,3)+1)/2)];
       tDir = currCirc - circInHamCyc;
       if tDir(1) == mergeDir(1) && tDir(2) == mergeDir(2) && tDir(3) == mergeDir(3)
           for j = 1:8
               if norm(hamCycOld(i,:) - pCC(j,:)) == 1
                   % found all nodes in the adjacent circuit of the hamCyc
                   adjPixInHamCyc(cnt,1:3) = hamCycOld(i,:);
                   adjPixInHamCyc(cnt,4) = i;
                   adjPixInCurrCirc(cnt,1:3) = pCC(j,:);
                   adjPixInCurrCirc(cnt,4) = j;
                   cnt = cnt + 1;
                   if cnt > 4
                       allPixFound = true;
                       break;
                   end
               end
           end
       end
       
       if allPixFound
           break;
       end
   end

   % check edges in adjPixInHamCyc
   eH = zeros(4,8); 
   cc = 1;
   for i = 1:size(adjPixInHamCyc,1)
       for j = i+1:size(adjPixInHamCyc,1)
            if abs(adjPixInHamCyc(i,4) - adjPixInHamCyc(j,4)) == 1
                eH(cc,1:4) = adjPixInHamCyc(i,1:4);
                eH(cc,5:8) = adjPixInHamCyc(j,1:4);
                cc = cc + 1;
            end
        end
   end
   lastE = find(~any(eH,2),1)-1;
   if lastE < 4 % remove invalid items
    eH(lastE+1:end,:) = [];
   end
   % generate a random configuration for the circuit to be merged
   allHamModes = randperm(6);
   for m = 1:length(allHamModes)
       hamMode = allHamModes(m); % get a random mode of the 3D cube
       circGraphAdjMat = getElementHamCube(hamMode);
       % find adjacent circuit
       edgesToBreak = findAdjParallelEdges(eH, adjPixInCurrCirc, pCC, circGraphAdjMat);

       if ~isempty(edgesToBreak)
           hamCycNew = parallelEdgeCombine(hamCycOld, edgesToBreak, pCC, circGraphAdjMat);
           return;
       end
   end 
   % there's really no way out
   % TODO: If there's branching at vertices, we cannot resolve now
    hamCycNew = nonParallelCombine(hamCycOld, adjPixInHamCyc, adjPixInCurrCirc, pCCMin, pCCMax);
end

% find adjacent parallel edges in the currCirc and hamCyc
function edgesToBreak = findAdjParallelEdges(eH, adjPixInCurrCirc, pixInCurrCirc, circGraphAdjMat)
    hasParallel = false;
    edgesToBreak = []; % records the actual edges and ids of their corresponding vertices
    for i = 1:size(adjPixInCurrCirc,1)
        t = adjPixInCurrCirc(i,4);
        circE = find(circGraphAdjMat(t,:));
        % get none zero edge id
        for j = 1:length(circE)
           s = circE(j);
           eC = pixInCurrCirc(t,:) - pixInCurrCirc(s,:);
           for kk = 1:size(eH,1) % check for each edge in eH
               % distance requirement
               
             %%%TODO: Make sure that the vertices in currCirc and hamCyc are
             %%not twisting!!!
               eCctr = 0.5 * (pixInCurrCirc(t,:) + pixInCurrCirc(s,:));
               eHkctr = 0.5 * (eH(kk,1:3) + eH(kk,5:7));
               
               if norm(eCctr - eHkctr) == 1 % adjacent
                   eHk = eH(kk,1:3)-eH(kk,5:7);
                   if norm(cross(eC, eHk)) == 0 % are they parallel 
                       if norm(pixInCurrCirc(t,:) - eH(kk,1:3)) == 1
                        edgesToBreak(end+1,:) = [pixInCurrCirc(t,:),t, pixInCurrCirc(s,:), s]; % First row is the edge in curr circuit!
                       else
                        edgesToBreak(end+1,:) = [pixInCurrCirc(s,:),s, pixInCurrCirc(t,:), t];  
                       end
                       edgesToBreak(end+1,:) = eH(kk,:); % Second row is the edge in hamCyc!
                       hasParallel = true;
                       return; % found one pair 
                   end
               end
           end
        end
    end
    
    if ~hasParallel
        edgesToBreak = [];
        return;
    end
end

% add curr circ if there're parallel edges!
function hamCycComb = parallelEdgeCombine(hamCyc, edgesToBreak, pixInCurrCirc, circGraphAdjMat)
    hamCycComb = hamCyc;
    startVCircId = edgesToBreak(1,4);
    endVCircId = edgesToBreak(1,8);
    currCycLen = find(~any(hamCyc,3),1)-1;
    
    edgeBrokeAdjMat = circGraphAdjMat;
    edgeBrokeAdjMat(startVCircId,endVCircId) = 0;
    edgeBrokeAdjMat(endVCircId,startVCircId) = 0;
    
    startVCycId = edgesToBreak(2,4);
    endVCycId = edgesToBreak(2,8);
    
    hamCycComb(1:startVCycId,:) = hamCyc(1:startVCycId,:);
    currId = startVCycId;
    %merge with the curr circ
    % find the hamPath
    hamPath = hamiltonian(edgeBrokeAdjMat, startVCircId, endVCircId);
    for i = 1:length(hamPath)
        vk = hamPath(i);
        currId = currId + 1;
        hamCycComb(currId,:) = pixInCurrCirc(vk,:);
    end
    hamCycComb(currId+1:currCycLen+8,:) = hamCyc(endVCycId:currCycLen,:);
end

function hamCycComb = nonParallelCombine(hamCyc, adjPixInHamCyc, adjPixInCurrCirc, pCCMin, pCCMax)
    hamCycComb = hamCyc;
    startId = 1;
    wid = 1;
    currCircSearchDir = cross(adjPixInCurrCirc(2,1:3) - adjPixInCurrCirc(1,1:3), adjPixInCurrCirc(3,1:3) - adjPixInCurrCirc(1,1:3)); 
    tryVcurrCirc = adjPixInCurrCirc(1,1:3) + currCircSearchDir;
    distToCircMin = tryVcurrCirc - pCCMin;
    if min(distToCircMin) < 0 || max(distToCircMin) > 1
        currCircSearchDir = -currCircSearchDir;
    end
    
    for breakPt = 1:2
%         startVCyc = adjPixInHamCyc(2*breakPt-1,1:3);
        startVCycId = adjPixInHamCyc(2*breakPt-1,4);
        
        % copy part before the merged region
        hamCycComb(wid:wid+startVCycId-startId,:) = hamCyc(startId:startVCycId,:);
        
        inPixInCurrCirc = adjPixInCurrCirc(2*breakPt-1,1:3);
        outPixInCurrCirc = adjPixInCurrCirc(2*breakPt,1:3);
        wid = wid+startVCycId-startId+1;
        
        hamCycComb(wid,:) = inPixInCurrCirc;
        hamCycComb(wid+1,:) = inPixInCurrCirc + currCircSearchDir; %pCC(); % next adjacent pix in currCirc
        hamCycComb(wid+2,:) = outPixInCurrCirc + currCircSearchDir; % next adjacent pix in currCirc
        hamCycComb(wid+3,:) = outPixInCurrCirc;  % next adjacent pix in currCirc
        
        % copy part after the merged region
%         endVCyc = adjPixInHamCyc(2*breakPt,1:3);
        endVCycId = adjPixInHamCyc(2*breakPt,4);
        
        copyLen = endVCycId - (startVCycId + 1) + 1;
        wid = wid+4;
        hamCycComb(wid:wid+copyLen-1,:) = hamCyc(startVCycId+1:endVCycId,:);
        
        startId = endVCycId+1;
        wid = wid+copyLen;
    end
    currCycLen = find(~any(hamCyc,3),1)-1;
    % copy the residual part in hamCyc if any
    if startId < currCycLen
        hamCycComb(wid:wid+currCycLen-startId,:) = hamCyc(startId:currCycLen,:);
    end
end
