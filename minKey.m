
function [minidx,parentidx, minDir] = minKey(zmstSet, zmstLast, zTaken,...
    WXR, WXL, WYU, WYD, dimX, dimY)
    global useLocWeight;
% find the smallest weight of current location idx
    minVal = Inf;
    minidx = [0,0];
    minDir = 0; %minvalue direction from current node
    
        
    % with location weight
%     useLocWeight = true;
    if useLocWeight
        alpha = 0.1; % blend term for location and other weights
    else
        alpha = 0;
    end
    
    
    for i = zmstLast:-1:1
       iid = zmstSet(i,:);
       idX = iid(2); idY = iid(1);
           % check 4 neighbors
           if idY-1 < 1% up 
           else
                if zTaken(idY-1,idX) == false &&...
                   (1-alpha)* WYU(idY-1,idX) + alpha*locTerm(zmstLast+1, [idY-1,idX], [dimY, dimX]) < minVal
               % WYU(idY-1,idX) < minVal
                   minVal = (1-alpha)* WYU(idY-1,idX) + alpha*locTerm(zmstLast+1, [idY-1,idX], [dimY, dimX]);%WYU(idY-1,idX);
                   minidx = [idY-1,idX];
                   parentidx = [idY,idX];
                   minDir = 1;
               end
              end    
           
           
           if idX+1 > dimX % right
           else
               if zTaken(idY,idX+1) == false &&...
                    (1-alpha)* WXR(idY,idX) + alpha * locTerm(zmstLast+1, [idY,idX], [dimY, dimX]) < minVal
%                    WXR(idY,idX) < minVal
                   minVal = (1-alpha)* WXR(idY,idX) + alpha * locTerm(zmstLast+1, [idY,idX], [dimY, dimX]);%WXR(idY,idX);
                   minidx = [idY,idX+1];
                   parentidx = [idY,idX];
                   minDir = 3;
               end
           end
           
           
           if idY+1 > dimY  % down
           else
                if zTaken(idY+1,idX) == false &&...
                    (1-alpha)* WYD(idY,idX) + alpha*locTerm(zmstLast+1, [idY,idX], [dimY, dimX]) < minVal
%                 WYD(idY,idX) < minVal
                   minVal =  (1-alpha)* WYD(idY,idX) + alpha*locTerm(zmstLast+1, [idY,idX], [dimY, dimX]);%WYD(idY,idX);
                   minidx = [idY+1,idX];
                   parentidx = [idY,idX];
                   minDir = 2;
               end
           end
           
           if idX-1<1 %left
           else
                if zTaken(idY,idX-1) == false &&...
                   (1-alpha)*WXL(idY,idX-1) + alpha*locTerm(zmstLast+1, [idY,idX-1], [dimY, dimX]) < minVal
%                 && WXL(idY,idX-1) < minVal
                   minVal =  (1-alpha)*WXL(idY,idX-1) + alpha*locTerm(zmstLast+1, [idY,idX-1], [dimY, dimX]);%WXL(idY,idX-1);
                   minidx = [idY,idX-1];
                   parentidx = [idY,idX];
                   minDir = 4;
               end
           end
    end
end


% the position weight term
function w = locTerm(visitId, pos, dim)
   
%     w = abs(visitId - ((dim(1) - pos(1)+1)*dim(2)+pos(2)))/(dim(1)*dim(2))*256;
    
    % get center of the volume
    normT = 1/(dim(1)*dim(2));
    Ctr = [dim(1)/2,dim(2)/2];
    CtrId = Ctr(1)*dim(2)+Ctr(2);
    % get distance to center of the volume
    distId2Ctr = visitId - CtrId;
    distPos2Ctr = (pos(1)-Ctr(1))*dim(2)+pos(2)-Ctr(2);
    w = abs(distId2Ctr - distPos2Ctr) * normT ;
     w = w * 255;
end