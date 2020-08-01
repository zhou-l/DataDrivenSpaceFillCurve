
%% iterative version that finds a low cost path among several paths using nearest neighbor heuristics
function [Status, hamPath, LT, exitPix] = hamPathGridRepNN(I, Source, Dsts)
    global neighbors;
    % Ending Condition check
    neighbors = [-1,0; %T
        0,1;%R
        1,0;%B
        0,-1;%L
        ];
    nodeDeg = length(neighbors); % max degrees of a vertex
    totalNodes = size(I,1)*size(I,2); 
       
    numPaths = length(Dsts);
    totalHamPaths = cell(numPaths,2);
    
        if Source(1)>size(I,1)
            Source(1) = Source(1) - size(I,1);
        end
        if Source(2)>size(I,2)
            Source(2) = Source(2) - size(I,2);
        end
    dimI = size(I);
    alpha = 0.9;%1; %blend term for value and locality terms


    for t = 1:size(Dsts,1)    
%     for t = 1:2:size(Dsts,1)
%     for t = 1:1
        %Todo: use totalWeights for better path
        Destination = Dsts(t,:);
        hamPath = zeros(totalNodes,2); %change dimensionality to 3 for 3D case!
        totalWeights = Inf .* ones(totalNodes,1);  
        Status = 0;
        stack = [];
        stack(end+1,:) = Source;
        LT = zeros(totalNodes,1);
        visited = zeros(size(I,1),size(I,2),nodeDeg);
    
        while ~isempty(stack)
            noWhereToGo = true;
            nodesFound = find(~any(hamPath,2), 1)-1;
            node = stack(1,:);
            stack(1,:) = [];
            if nodesFound >= 1
                prevNode = hamPath(nodesFound,:);
                toCur = node - prevNode;
                dir = -1;
                for jj = 1:nodeDeg
                    if toCur(1) == neighbors(jj,1) && toCur(2) == neighbors(jj,2)
                        dir = jj;
                        break;
                    end
                end
                if dir >= 1
                    visited(prevNode(1), prevNode(2), dir) = 1;
                else
                    disp('Something is wrong!');
                end
            end
            
            hamPath(nodesFound+1,:) = node;   
            % get weight from this node to previous node
            if nodesFound == 0 
                totalWeights(nodesFound+1,:) = 0;
            else
                totalWeights(nodesFound+1,:) = alpha * pathWeightFunc(I, hamPath(nodesFound,:), node) +...
           (1 - alpha) * locTerm(nodesFound+1, node, dimI);
            end
            
            LT(nodesFound+1,:) = I(node(1),node(2));
    %         totalWeights = totalWeights + currWeight;
             % sort next node by the cost
             weights = Inf * ones(nodeDeg,1);
             for i = 1:nodeDeg
                 nextNode = node+neighbors(i,:);
                 if nextNode(1)>=1 && nextNode(1)<=size(I,1) && nextNode(2)>=1 && nextNode(2)<=size(I,2)
                    weights(i) = alpha * pathWeightFunc(I, node, nextNode)+...
                       (1 - alpha) * locTerm(nodesFound+1, node, dimI); 
                 end
             end
             [weights, sInd] = sort(weights,'descend');

             % Greedy strategy to check next node
             for i= 1:1:length(sInd) %1:length(sInd) % check four neighbors only
    %           for i = 1:1:totalNodes
                if weights(i) == Inf
                    continue; % if weight is infinity, skip the loop
                end
                nextNode = node+ neighbors(sInd(i),:);

                if nextNode(1)==Destination(1) && nextNode(2)==Destination(2)
                    visited(node(1), node(2), sInd(i)) = 1; % set visit direction
                    lastN = find(~any(hamPath,2),1)-1;%nnz(hamPath);

                    if ( (lastN == totalNodes-1 && (Source(1)~=Destination(1)||Source(2)~=Destination(2)) ) ||...
                          (lastN == totalNodes && Source(1)==Destination(1) && Source(2) == Destination(2)) ) % check finish first!
    %                     if ( Graph(hamPath(lastN), Destination) ~= 0)
                            hamPath(lastN+1, :) = Destination;
                            LT(lastN+1,:) = I(Destination(1),Destination(2));
                            totalWeights(lastN+1,:) = alpha * pathWeightFunc(I, hamPath(lastN,:), Destination) +...
                                (1 - alpha) * locTerm(lastN+1, Destination, dimI);
                            
                            Status = 1;
                            totalHamPaths{t,1} = hamPath;
                            totalHamPaths{t,2} = sum(totalWeights,1);
                            break; % should always have a path?
    %                     else
    %                         Status = 0;
    %                         return;
    %                     end
                    else
                        continue;
                    end
                end
                
                if Status == 1
                    break;
                end
                
                if isPathSafeGrid(hamPath, nextNode)
                     noWhereToGo = false;
                     stack = [nextNode;stack];            
                else
                     visited(node(1),node(2), sInd(i)) = 1;
                end

             end
             
            if Status == 1
                break; % done for this path
            end
             
            if noWhereToGo 
    %             k = nnz(hamPath);%remove last node from path
                k = find(~any(hamPath,2),1)-1;
                ln = hamPath(k,:);
                visited(ln(1),ln(2),:) = 0;%reset flags
                hamPath(k, :) = [0,0]; % this path is deadend
                LT(k,:) = 0;
                Status = 0;
                % back track to node that's not fully visited
                for l = k - 1:-1:1
                    nk = hamPath(l,:);
                    deg = getDeg(nk, size(I,1), size(I,2));
                    % corner and boundaries have different
                    sv = sum(visited(nk(1),nk(2),:));
                    if sv >= deg % four connectivity 
                        visited(nk(1), nk(2),:) = 0;%reset flags
                        hamPath(l,:) = [0,0];
                        LT(l,:) = 0;
                        % remove weight
                        totalWeights(l,1) = 0;
                    else
                        break;
                    end
                end
%                 disp(hamPath);
            end
        end
        if Status == 0
            totalHamPaths{t,1} = [];
            totalHamPaths{t,2} = Inf;
        end
    end
    
    
    % find min cost path
    minCost = Inf;
    minId = -1;
    for i = 1:length(totalHamPaths)
        if totalHamPaths{i,2} < minCost
            minCost = totalHamPaths{i,2};
            minId = i;
        end
    end
    % get result
    if minCost ~= Inf
        hamPath = totalHamPaths{minId, 1};
        Status = 1;
        exitPix = Dsts(minId,:);
    else
        hamPath = [];
        Status = 0;
        exitPix = Dsts(minId,:);
    end
end

function deg = getDeg(node, gridDimY, gridDimX)
    if node(1) == 1 || node(1) == gridDimY
        if node(2) == 1 || node(2) == gridDimX
            deg = 2;
        else
            deg = 3;
        end
        return;
    end
    if node(2) == 1 || node(2) == gridDimX
        if node(1) == 1 || node(1) == gridDimY
            deg = 2;
        else
            deg = 3;
        end
        return;
    end
    deg = 4;
end

function w = pathWeightFunc(I, node, nextNode)
   w = abs(-I(node(1),node(2)) + I(nextNode(1),nextNode(2)));
   
end

function w = locTerm(visitId, pos, dim)
   
%     w = abs(visitId - ((dim(1) - pos(1)+1)*dim(2)+pos(2)))/(dim(1)*dim(2))*256;
    
    % get center of the volume
    normT = 1/sqrt(dim(1)*dim(2));
    Ctr = [dim(1)/2,dim(2)/2];
    CtrId = Ctr(1)*dim(2)+Ctr(2);
    % get distance to center of the volume
    distId2Ctr = visitId - CtrId;
    distPos2Ctr = norm(pos - Ctr); %(pos(1)-Ctr(1))*dim(2)+pos(2)-Ctr(2);
    w = abs(distId2Ctr - distPos2Ctr) * normT ;
     w = w;
end

function Flag = isPathSafeGrid(hamPath, testNode)
    lastN = find(~any(hamPath,2),1)-1;
    diff = testNode - hamPath(lastN,:);
    if norm(diff) ~= 1
        Flag = 0;
        return;
    end

    for ii=1:lastN%nnz(hamPath)
        if hamPath(ii,1) == testNode(1) && hamPath(ii,2) == testNode(2)
            Flag = 0;
            return;
        end
    end
    Flag = 1;
end