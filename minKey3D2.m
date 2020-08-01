
function [minidx,parentidx, minDir] = minKey3D2(zmstSet, zmstLast, zTaken,...
    WXR, WXL, WYU, WYD, WZB, WZF, GXR, GXL, GYU, GYD, GZB, GZF, dimX, dimY, dimZ)
    global useLocWeight;
% find the smallest weight of current location idx
    minVal = Inf;
    minidx = [0,0,0];
    minDir = 0; %minvalue direction from current node
    
        
    % with location weight
%     useLocWeight = true;
    if useLocWeight
        alpha = 0.1;%0.1; % blend term for location and other weights
    else
        alpha = 0;
    end
%     alpha = 0.3;
    alpha = 0.1;%05;
    alpha = 0;
    beta = 1; %0.9;
    for i = zmstLast:-1:1
       iid = zmstSet(i,:);
       idX = iid(2); idY = iid(1);idZ = iid(3);
           % check 6 neighbors
           if idY-1 < 1% up 
           else
                if zTaken(idY-1,idX,idZ) == false &&...
                   (1-alpha)* (beta * WYU(idY-1,idX,idZ) + (1-beta)*GYU(idY-1,idX,idZ)) + alpha*locTerm3D(zmstLast+1, [idY-1,idX,idZ], [dimY, dimX,dimZ]) < minVal
               % WYU(idY-1,idX) < minVal
                   minVal = (1-alpha)* (beta * WYU(idY-1,idX,idZ) + (1-beta)*GYU(idY-1,idX,idZ)) + alpha*locTerm3D(zmstLast+1, [idY-1,idX,idZ], [dimY, dimX,dimZ]);%WYU(idY-1,idX);
                   minidx = [idY-1,idX,idZ];
                   parentidx = [idY,idX,idZ];
                   minDir = 1;
               end
              end    
           
           
           if idX+1 > dimX % right
           else
               if zTaken(idY,idX+1,idZ) == false &&...
                    (1-alpha)* (beta * WXR(idY,idX,idZ)+(1-beta)*GXR(idY,idX,idZ)) + alpha * locTerm3D(zmstLast+1, [idY,idX,idZ], [dimY, dimX,dimZ]) < minVal
%                    WXR(idY,idX) < minVal
                   minVal = (1-alpha)* (beta * WXR(idY,idX,idZ)+(1-beta)*GXR(idY,idX,idZ)) + alpha * locTerm3D(zmstLast+1, [idY,idX,idZ], [dimY, dimX,dimZ]);%WXR(idY,idX);
                   minidx = [idY,idX+1,idZ];
                   parentidx = [idY,idX,idZ];
                   minDir = 3;
               end
           end
           
           
           if idY+1 > dimY  % down
           else
                if zTaken(idY+1,idX,idZ) == false &&...
                    (1-alpha)* (beta*WYD(idY,idX,idZ)+(1-beta)*GYD(idY,idX,idZ)) + alpha*locTerm3D(zmstLast+1, [idY,idX,idZ], [dimY, dimX,dimZ]) < minVal
%                 WYD(idY,idX) < minVal
                   minVal = (1-alpha)* (beta*WYD(idY,idX,idZ)+(1-beta)*GYD(idY,idX,idZ)) + alpha*locTerm3D(zmstLast+1, [idY,idX,idZ], [dimY, dimX,dimZ]);%WYD(idY,idX);
                   minidx = [idY+1,idX,idZ];
                   parentidx = [idY,idX,idZ];
                   minDir = 2;
               end
           end
           
           if idX-1<1 %left
           else
                if zTaken(idY,idX-1,idZ) == false &&...
                   (1-alpha)*(beta*WXL(idY,idX-1,idZ)+(1-beta)*GXL(idY,idX-1,idZ)) + alpha*locTerm3D(zmstLast+1, [idY,idX-1,idZ], [dimY, dimX,dimZ]) < minVal
%                 && WXL(idY,idX-1) < minVal
                   minVal =  (1-alpha)*(beta*WXL(idY,idX-1,idZ)+(1-beta)*GXL(idY,idX-1,idZ)) + alpha*locTerm3D(zmstLast+1, [idY,idX-1,idZ], [dimY, dimX,dimZ]);
                   minidx = [idY,idX-1,idZ];
                   parentidx = [idY,idX,idZ];
                   minDir = 4;
               end
           end
           
           if idZ+1>dimZ %back
           else
                if zTaken(idY,idX,idZ+1) == false &&...
                   (1-alpha)*(beta*WZB(idY,idX,idZ)+(1-beta)*GZB(idY,idX,idZ)) + alpha*locTerm3D(zmstLast+1, [idY,idX,idZ], [dimY, dimX,dimZ]) < minVal
%                 && WXL(idY,idX-1) < minVal
                   minVal =  (1-alpha)*(beta*WZB(idY,idX,idZ)+(1-beta)*GZB(idY,idX,idZ)) + alpha*locTerm3D(zmstLast+1, [idY,idX,idZ], [dimY, dimX,dimZ]);
                   minidx = [idY,idX,idZ+1];
                   parentidx = [idY,idX,idZ];
                   minDir = 5;
               end
           end
           
           if idZ-1<1 %left
           else
                if zTaken(idY,idX,idZ-1) == false &&...
                   (1-alpha)*(beta*WZF(idY,idX,idZ-1)+(1-beta)*GZF(idY,idX,idZ-1)) + alpha*locTerm3D(zmstLast+1, [idY,idX,idZ-1], [dimY, dimX,dimZ]) < minVal
%                 && WXL(idY,idX-1) < minVal
                   minVal =  (1-alpha)*(beta*WZF(idY,idX,idZ-1)+(1-beta)*GZF(idY,idX,idZ-1)) + alpha*locTerm3D(zmstLast+1, [idY,idX,idZ-1], [dimY, dimX,dimZ]);
                   minidx = [idY,idX,idZ-1];
                   parentidx = [idY,idX,idZ];
                   minDir = 6;
               end
           end
    end
end


% the position weight term
function w = locTerm3D(visitId, pos, dim)
   
%     w = abs(visitId - ((dim(1) - pos(1)+1)*dim(2)+pos(2)))/(dim(1)*dim(2))*256;
    
%     % get center of the volume
%     normT = 1/(dim(1)*dim(2)*dim(3));
%     Ctr = [dim(1)/2,dim(2)/2,dim(3)/2];
%     CtrId = (Ctr(3)*dim(3) + Ctr(1))*dim(2)+Ctr(2);
%     % get distance to center of the volume
%     distId2Ctr = visitId - CtrId;
%     distPos2Ctr = ((pos(3)-Ctr(3))*dim(1)+(pos(1)-Ctr(1)))*dim(2)+pos(2)-Ctr(2);
%     w = abs(distId2Ctr - distPos2Ctr) * normT ;
%      w = w * 255;
     
      % get center of the volume
    posF = pos ./ dim;
    bkSize = 1/2; % 1/8 of the size
    portY = floor(posF(1) / bkSize);
    portX = floor(posF(2) / bkSize);
    portZ = floor(posF(3) / bkSize);
    CtrThis = [dim(1)*(portY+0.5)*bkSize, dim(2)*(portX+0.5)*bkSize, dim(3)*(portZ+0.5)*bkSize];
    CtrThisId = (CtrThis(3) * dim(1) +CtrThis(1))*dim(2) + CtrThis(2);
    normT = 1/(bkSize*bkSize*bkSize*dim(1)*dim(2)*dim(3));
    distId2PosId = abs(visitId - CtrThisId);
    distPos2Ctr = norm(pos - CtrThis); %(pos(1)-Ctr(1))*dim(2)+pos(2)-Ctr(2);
    normC = 1/norm(bkSize*bkSize*bkSize*dim(1)*dim(2)*dim(3));
%     w = abs(distId2PosId * normT - distPos2Ctr*normC)/factor;
%     w = abs(distId2PosId * normT -  distPos2Ctr*normC)*255;
    w = abs(distPos2Ctr*normC)*255;
    return;
    
end