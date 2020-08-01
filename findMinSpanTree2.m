function [T,mstSet] = findMinSpanTree2(WpXR, WpXL, WpYU, WpYD, GpXR, GpXL, GpYU, GpYD, firstNodePos)

    [dimY,dimXx,dimZ] = size(WpXR);  
    [dimYy,dimX,dimZ] = size(WpYU);
    T = zeros(dimY,dimX,4); %quadtree records the MST graph node ordering--each node has 4 children: parent->next
    posTaken = false(dimY,dimX,dimZ);
    mstSet = zeros(dimX*dimY*dimZ,2);
    if nargin < 9
        nextid = [dimY,1]; %first node
    else
        nextid = firstNodePos;
    end
    mstSet(1,:) = nextid;% first node
    posTaken(nextid(1),nextid(2),1) = true;
    for mstLast = 1:dimX * dimY-1 %remaining nodes
        
       [nextid, parentid, minDir] = minKey2(mstSet, mstLast, posTaken,...
           WpXR, WpXL, WpYU, WpYD, GpXR, GpXL, GpYU, GpYD, dimX, dimY);
       mstSet(mstLast+1,:) = nextid;
       posTaken(nextid(1),nextid(2),1) = true;
       
       % set the direction for parent id, set the minDir channel to 1
       T(parentid(1),parentid(2),minDir) = 1; 
       
    end
end