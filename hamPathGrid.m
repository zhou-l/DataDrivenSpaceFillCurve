%% iterative version
function [Status, hamPath, LT] = hamPathGrid(I, Source, Destination)
global neighbors;
    % Ending Condition check
    neighbors = [-1,0; %T
        0,1;%R
        1,0;%B
        0,-1;%L
%         -1,1;%TR
%         1,1;%BR
%         1,-1;%BL
%         -1,-1%TL
        ];
      nodeDeg = length(neighbors); % need to loose the restrictions!
     totalNodes = size(I,1)*size(I,2);
    hamPath = zeros(totalNodes,2); %change dimensionality to 3 for 3D case!
    LT = zeros(totalNodes,1);
    visited = zeros(size(I,1),size(I,2),nodeDeg);
   
    stack = [];
    stack(end+1,:) = Source;
    %Todo: use totalWeights for better path
    totalWeights = 0;  
    while ~isempty(stack)
        noWhereToGo = true;
        nodesFound = find(~any(hamPath,2), 1)-1;
%         if ( (nodesFound == totalNodes-1 && Source~=Destination) || (nodesFound == totalNodes && Source==Destination) )
%             if ( Graph(hamPath(nodesFound), Destination) ~= 0)
%                 hamPath(nodesFound+1) = Destination;
%                 Status = 1;
%                 return;
%             else
%                 Status = 0;
%                 return;
%             end
%         end
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
        LT(nodesFound+1,:) = I(node(1),node(2));
%         totalWeights = totalWeights + currWeight;
         % sort next node by the cost
         weights = Inf * ones(nodeDeg,1);
         for i = 1:nodeDeg
             nextNode = node+neighbors(i,:);
             if nextNode(1)>=1 && nextNode(1)<=size(I,1) && nextNode(2)>=1 && nextNode(2)<=size(I,2)
                weights(i) = abs(-I(node(1),node(2)) + I(nextNode(1),nextNode(2)));
             end
         end
         [weights, sInd] = sort(weights,'descend');
         
         % Greedy strategy to check next node
         for i= 1:length(sInd) % check four neighbors only
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
                        Status = 1;
                        return; % should always have a path?
%                     else
%                         Status = 0;
%                         return;
%                     end
                end

                continue;
            end
            
            
            if isPathSafeGrid(hamPath, nextNode)
                 noWhereToGo = false;
                 stack = [nextNode;stack];            
            else
                 visited(node(1),node(2), sInd(i)) = 1;
            end

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
                else
                    break;
                end
            end
        end
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