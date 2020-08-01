%%
% Let us create the following graph
%       (1)--(2)--(3)-------(4)
%        |   / \   |		 |
%        |  /   \  |		 |
%        | /     \ |		 |
%       (5)-------(6)		 |
%		 |		 			 |
%		 |					 |
%		 |					 |
%		(7)-----------------(8)
%   
% g=[0 1 0 0 1 0 0 0;
%    1 0 1 0 1 1 0 0;
%    0 1 0 1 0 1 0 0;
%    0 0 1 0 0 0 0 1;
%    1 1 0 0 0 1 1 0;
%	 0 1 1 0 1 0 0 0;
%	 0 0 0 0 1 0 0 1;
%	 0 0 0 1 0 0 1 0]
% s=5; % Source
% d=1; % Destination
%
% P = hamiltonianPath(g,s,d);
%
% P will be an array mentioning the path/cycle, if path/cycle found; or a
% string: 'No Path Found', if path/cycle not found
%
% #Note: This code can be used for finding Hamiltonian cycle also. For
% that, make sure Source and Destination are same.
%%
%{
    Main Function
%}
function hamPath = hamiltonian(Graph, Source, Destination)
% Input Checking
if ~isreal(Graph)
    error('Graph must be in real form');
elseif ~isnumeric(Graph)
    error('Matrix must be numeric');
elseif ~ismatrix(Graph)
    error('Check Matrix Dimensions');
else
    [r, c] = size (Graph);
    if r~=c
        error('Matrix must be square matrix');
    end
end
if ~(isreal(Source)||isreal(Destination)||(Source>0 && Source<=r) || (Destination>0 && Destination<=r))
    error('improper Source/Destination');
end
clear c;
% Function call
hamPath = findHam(Graph, Source, Destination, r);
end
%%
%{
    This functions sets some initial parameters, and calls the actual
    function.
%}
function hamPath = findHam(Graph, Source, Destination, totalNodes)
hamPath = zeros(size(Graph(1,:)));
hamPath(1) = Source;
global recCnt;
recCnt = 0;
global MAX_CNT;
MAX_CNT = size(Graph,1)*size(Graph,2);%power(max(degrees),size(Graph,1));
MAX_CNT = MAX_CNT * 15;
%    [Status, hamPath] = hamRec(Graph, hamPath, Source, Destination, totalNodes, 1);
%  hamPath
% iterative version
[Status, hamPath] = hamRecNR(Graph, hamPath, Source, Destination, totalNodes, 1);
if Status == 0
    if Source ~= Destination
        hamPath =[];
        disp('No Path Found');
    else
        hamPath = [];
        disp('No Cycle Found');
    end
    return;
end
end
%%
%{
    This function recursively call itself, hence finding the solution
%}
function [Status, hamPath] = hamRec(Graph, hamPath, Source, Destination, totalNodes, nodesFound)
% Ending Condition check
global recCnt;
global MAX_CNT;
recCnt = recCnt + 1;
if ( (nodesFound == totalNodes-1 && Source~=Destination) || (nodesFound == totalNodes && Source==Destination) )
    if ( Graph(hamPath(nodesFound), Destination) ~= 0)
        hamPath(nodesFound+1) = Destination;
        Status = 1;
        return;
    else
        Status = 0;
        return;
    end
end
for i=1:totalNodes
    if i==Destination
        continue;
    end
    if recCnt > MAX_CNT
        Status = 0;
        return;
    end
    if isSafe(Graph, hamPath, nodesFound, i)
        hamPath(nodesFound+1) = i;
        
        [Status, hamPath] = hamRec(Graph, hamPath, Source, Destination, totalNodes, nodesFound+1);
        if Status
            return;
        end
        
        hamPath(nodesFound+1) = 0;
    end
end
Status = 0;
end

%% iterative version
function [Status, hamPath] = hamRecNR(Graph, hamPath, Source, Destination, totalNodes, nodesFound)
    % Ending Condition check
   global MAX_CNT;
    hamPath = zeros(size(Graph(1,:)));
    degrees = sum(Graph,2);
    % set an empirical upper bound of iterations
%     MAX_CNT = size(Graph,1)*size(Graph,2);%power(max(degrees),size(Graph,1));
%     MAX_CNT = MAX_CNT * 10;
    visited = zeros(size(Graph,1),size(Graph,2));
     stack = [];
     stack(end+1,:) = Source;
     cnt = 0;
    while ~isempty(stack)
        cnt = cnt + 1;
        if cnt > MAX_CNT 
           Status = 0;
           return;
        end
        noWhereToGo = true;
        nodesFound = nnz(hamPath);
        if ( (nodesFound == totalNodes-1 && Source~=Destination) || (nodesFound == totalNodes && Source==Destination) )
            if ( Graph(hamPath(nodesFound), Destination) ~= 0)
                hamPath(nodesFound+1) = Destination;
                Status = 1;
                return;
            else
                Status = 0;
                return;
            end
        end
        node = stack(1,:);
        stack(1,:) = [];
        if nodesFound >= 1
            visited(hamPath(nodesFound), node) = 1;
        end
        hamPath(nodesFound+1) = node;   

        
         for i= totalNodes:-1:1
%           for i = 1:1:totalNodes
            if isChild(Graph, node, i)
                if i==Destination
                    visited(node, i) = 1;
                    lastN = nnz(hamPath);
                    if ( (lastN == totalNodes-1 && Source~=Destination) || (lastN == totalNodes && Source==Destination) ) % check finish first!
                        if ( Graph(hamPath(lastN), Destination) ~= 0)
                            hamPath(lastN+1) = Destination;
                            Status = 1;
                            return;
                        else
                            Status = 0;
                            return;
                        end
                    end

                    continue;
                end
                if isPathSafe(hamPath, i)
                     noWhereToGo = false;
                     stack = [i;stack];            
                else
                      visited(node, i) = 1;
                end
            end
        end
        if noWhereToGo 
            k = nnz(hamPath);%remove last node from path
            visited(hamPath(k),:) = 0;%reset flags
            hamPath(k) = 0; % this path is deadend

            Status = 0;
            % back track to node that's not fully visited
            for l = k - 1:-1:1
                nk = hamPath(l);
                sv = sum(visited(nk,:));
                if sv >= degrees(nk) 
                    visited(nk,:) = 0;%reset flags
                    hamPath(l) = 0;
                else
                    break;
                end
            end
        end
    end

end

function Flag = isChild(Graph, curNode, testNode)
    if Graph(curNode,testNode) == 0
        Flag = 0;
    else
        Flag = 1;
    end
end

function Flag = isPathSafe(hamPath, testNode)
for ii=1:nnz(hamPath)
    if hamPath(ii) == testNode
        Flag = 0;
        return;
    end
end
Flag = 1;
end
%%
%{
    This function is used to check whether the current node can be added
    or not for making the path/cycle.
%}
function Flag = isSafe(Graph, hamPath, nodesFound, i)
if Graph(hamPath(nodesFound),i) == 0
    Flag = 0;
    return;
end
for ii=1:nodesFound
    if hamPath(ii) == i
        Flag = 0;
        return;
    end
end
Flag = 1;
end