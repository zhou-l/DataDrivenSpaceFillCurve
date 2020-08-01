function [T,mstSet] = findMinSpanTree3D2(WpXR, WpXL, WpYU, WpYD, WpZB, WpZF, GpXR, GpXL, GpYU, GpYD, GpZB, GpZF, firstNodePos)
    [dimY,dimXx,dimZ] = size(WpXR);  
    [dimYy,dimX,dimZ] = size(WpYU);
    T = zeros(dimY,dimX,dimZ,6); %quadtree records the MST graph node ordering--each node has 4 children: parent->next
    posTaken = false(dimY,dimX,dimZ);
    mstSet = zeros(dimX*dimY*dimZ,3);
    if nargin < 13
        nextid = [dimY,1,1]; %first node
    else
        nextid = firstNodePos;
    end
    mstSet(1,:) = nextid;% first node
    posTaken(nextid(1),nextid(2),nextid(3)) = true;
    for mstLast = 1:dimX*dimY*dimZ-1 %remaining nodes
        
       [nextid, parentid, minDir] = minKey3D2(mstSet, mstLast, posTaken,...
           WpXR, WpXL, WpYU, WpYD, WpZB, WpZF, GpXR, GpXL, GpYU, GpYD, GpZB, GpZF, dimX, dimY, dimZ);
       mstSet(mstLast+1,:) = nextid;
       posTaken(nextid(1),nextid(2),nextid(3)) = true;
       
       % set the direction for parent id, set the minDir channel to 1
       T(parentid(1),parentid(2),parentid(3),minDir) = 1; 
       
    end
end