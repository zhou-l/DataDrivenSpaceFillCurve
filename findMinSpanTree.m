function [T,mstSet] = findMinSpanTree(WpXR, WpXL, WpYU, WpYD, firstNodePos)
    [dimY,dimXx,dimZ] = size(WpXR);  
    [dimYy,dimX,dimZ] = size(WpYU);
    T = zeros(dimY,dimX,4); %quadtree records the MST graph node ordering--each node has 4 children: parent->next
    posTaken = false(dimY,dimX,dimZ);
    mstSet = zeros(dimX*dimY*dimZ,2);
    if nargin < 5
        nextid = [dimY,1]; %first node
    else
        nextid = firstNodePos;
    end
    mstSet(1,:) = nextid;% first node
    posTaken(nextid(1),nextid(2),1) = true;
    for mstLast = 1:dimX * dimY-1 %remaining nodes
        
       [nextid, parentid, minDir] = minKey(mstSet, mstLast, posTaken,...
           WpXR, WpXL, WpYU, WpYD, dimX, dimY);
       mstSet(mstLast+1,:) = nextid;
       posTaken(nextid(1),nextid(2),1) = true;
       
       % set the direction for parent id, set the minDir channel to 1
       T(parentid(1),parentid(2),minDir) = 1; 
       
%        % set the direction for parent id
%        % set after existing directions
%        for d=1:4
%            if T(parentid(1),parentid(2),d)== 0
%                T(parentid(1),parentid(2),d) = minDir;
%                break;
%            end
%        end
       
%        T(mstLast,:) = [parentid,nextid];
%        
%        fprintf('[%i,%i]->[%i,%i]\n', T(mstLast,1), T(mstLast,2), T(mstLast,3), T(mstLast,4));
    end
end