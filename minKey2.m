
function [minidx,parentidx, minDir] = minKey2(zmstSet, zmstLast, zTaken,...
    WXR, WXL, WYU, WYD, GXR, GXL, GYU, GYD, dimX, dimY)
    global useLocWeight;
    global zalpha;
    global gmstSet;
    gmstSet = zmstSet;
% find the smallest weight of current location idx
    minVal = Inf;
    minidx = [0,0];
    minDir = 0; %minvalue direction from current node
    
        
    % with location weight
%     useLocWeight = true;
    if useLocWeight
        alpha = zalpha; % blend term for location and other weights
    else
        alpha = 0;
    end
%     beta = 0.9;
    beta = 1;
    global factor;
    factor = 1;%0.95;
    for i = zmstLast:-1:1
       iid = zmstSet(i,:);
       idX = iid(2); idY = iid(1);
           % check 4 neighbors
           if idY-1 < 1% up 
           else
                if zTaken(idY-1,idX) == false &&...
                   (1-alpha)*factor* (beta * WYU(idY-1,idX) + (1-beta)*GYU(idY-1,idX)) + alpha*locTerm(zmstLast+1, [idY-1,idX], [dimY, dimX]) < minVal
               % WYU(idY-1,idX) < minVal
                   minVal = (1-alpha)* (beta * WYU(idY-1,idX) + (1-beta)*GYU(idY-1,idX)) + alpha*locTerm(zmstLast+1, [idY-1,idX], [dimY, dimX]);%WYU(idY-1,idX);
                   minidx = [idY-1,idX];
                   parentidx = [idY,idX];
                   minDir = 1;
               end
              end    
           
           
           if idX+1 > dimX % right
           else
               if zTaken(idY,idX+1) == false &&...
                    (1-alpha)*factor* (beta * WXR(idY,idX)+(1-beta)*GXR(idY,idX)) + alpha * locTerm(zmstLast+1, [idY,idX], [dimY, dimX]) < minVal
%                    WXR(idY,idX) < minVal
                   minVal = (1-alpha)* (beta * WXR(idY,idX)+(1-beta)*GXR(idY,idX)) + alpha * locTerm(zmstLast+1, [idY,idX], [dimY, dimX]);%WXR(idY,idX);
                   minidx = [idY,idX+1];
                   parentidx = [idY,idX];
                   minDir = 3;
               end
           end
           
           
           if idY+1 > dimY  % down
           else
                if zTaken(idY+1,idX) == false &&...
                    (1-alpha)*factor* (beta*WYD(idY,idX)+(1-beta)*GYD(idY,idX)) + alpha*locTerm(zmstLast+1, [idY,idX], [dimY, dimX]) < minVal
%                 WYD(idY,idX) < minVal
                   minVal = (1-alpha)* (beta*WYD(idY,idX)+(1-beta)*GYD(idY,idX)) + alpha*locTerm(zmstLast+1, [idY,idX], [dimY, dimX]);%WYD(idY,idX);
                   minidx = [idY+1,idX];
                   parentidx = [idY,idX];
                   minDir = 2;
               end
           end
           
           if idX-1<1 %left
           else
                if zTaken(idY,idX-1) == false &&...
                   (1-alpha)*factor*(beta*WXL(idY,idX-1)+(1-beta)*GXL(idY,idX-1)) + alpha*locTerm(zmstLast+1, [idY,idX-1], [dimY, dimX]) < minVal
%                 && WXL(idY,idX-1) < minVal
                   minVal =  (1-alpha)*(beta*WXL(idY,idX-1)+(1-beta)*GXL(idY,idX-1)) + alpha*locTerm(zmstLast+1, [idY,idX-1], [dimY, dimX]);
                   minidx = [idY,idX-1];
                   parentidx = [idY,idX];
                   minDir = 4;
               end
           end
    end
end


% the position weight term
function w = locTerm(visitId, pos, dim)
   global factor;
   global gmstSet;
%    w = abs(visitId - ((dim(1) - pos(1)+1)*dim(2)+pos(2)))/(dim(1)*dim(2))*256;
%     neighborhood = 4 * 4;
%     maxDist = [0 0];
%     minId = max(visitId - neighborhood,1);
%     idO = [1 1];
%     for i = minId:visitId
%         iid = gmstSet(i,:);
%         if i == minId
%             idO = iid;
%         end
%        maxDist(1) = max(iid(1) - idO(1), maxDist(1));
%        maxDist(2) = max(iid(2) - idO(2), maxDist(2));
%     end
%     w = max(norm(maxDist), norm(pos - idO))/ neighborhood;
%     return;
    
    % get center of the volume
    posF = pos ./ dim;
    bkSize = 1/4; % 1/8 of the size
    portY = floor(posF(1) / bkSize);
    portX = floor(posF(2) / bkSize);
    CtrThis = [dim(1)*(portY+0.5)*bkSize, dim(2)*(portX+0.5)*bkSize];
    
    CtrThisId = CtrThis(1)*dim(2) + CtrThis(2);
    
    % get the closest block center ID
    numBkY = 1.0/bkSize;
    numBkX = 1.0/bkSize;
    
%     minDist = realmax;
%     for yy = 0:numBkY-1
%         for xx = 0:numBkX-1
%             bkCtr = [dim(1)*(yy+0.5)*bkSize, dim(2)*(xx+0.5)*bkSize];
%             bkCtrId = bkCtr(1)*dim(1)+bkCtr(2);
%             if (abs(visitId - bkCtrId) < minDist)
%                 minDist = abs(visitId - bkCtrId);
%                 minDistBlkId = yy*numBkX + xx;
%             end
%         end
%     end
%     
%     normT = 1/(dim(1)*dim(2));
%     distId2PosId = abs(visitId - CtrThisId);
%     distId2PosId = minDist;
    
    distPos2Ctr = norm(pos - CtrThis); %(pos(1)-Ctr(1))*dim(2)+pos(2)-Ctr(2);
%     normC = 1/norm(dim(1)*dim(2));
    normC = 1/sqrt(bkSize*dim(1)*bkSize*dim(2));
    normT = normC;
%     w = abs(distId2PosId * normT - distPos2Ctr*normC)/factor;
%      w = abs(distId2PosId * normT -  distPos2Ctr*normC)/factor;
    w = abs(distPos2Ctr*normC)/factor;
    return;
    
    %% old codes
%     normT = 1/(dim(1)*dim(2));
%     Ctr = [dim(1)/2,dim(2)/2];
%     CtrId = Ctr(1)*dim(2)+Ctr(2);
%     % get distance to center of the volume
%     distId2Ctr = abs(visitId - CtrId);
%         posId = (pos(1)*dim(2)+pos(2));
%         distId2PosId = abs(visitId -  CtrId);
%     distPos2Ctr = norm(pos - Ctr); %(pos(1)-Ctr(1))*dim(2)+pos(2)-Ctr(2);
%     normC = 1/norm(Ctr);
% %     w = abs(distId2PosId * normT - distPos2Ctr*normC)/factor;
%     w = abs(distId2PosId * normT -  distPos2Ctr*normC)/factor;
% %     w = abs(distId2Ctr* normT - distPos2Ctr*normC)/factor;
% %     w = abs(distPos2Ctr*normC);
% %     w = abs(distId2Ctr* normT);
% %      w = w * 255;
end