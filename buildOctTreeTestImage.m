% build the octree test image
% V - the aggregated volume
% PTS - the point data
function [V,PTS, dimX,dimY,dimZ] = buildOctTreeTestImage(caseNum)



    PTS = [];
    blksize = 8;
    blk = zeros(blksize,blksize,blksize);
  
    for i = 1:4
        blk(:,:,i) =... 
                [ 1 3 4 3 0 0 0 0;
                  0 1 1 2 0 0 0 0;
                  0 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 0 0;
                  0 0 0 0 3 1 4 4;
                  0 0 0 0 1 2 4 4;
                  0 0 0 0 2 2 3 3;
                  0 0 0 0 2 2 3 3
                ];
    end
    for i = 5:8
        blk(:,:,i)   =...
                [ 0 0 4 3 3 1 4 4;
                  0 0 1 2 1 2 4 4;
                  0 0 0 0 2 2 3 3;
                  0 0 0 0 1 2 3 4;
                  0 0 0 0 3 1 4 4;
                  0 0 0 0 1 2 4 4;
                  0 0 0 0 2 2 3 3;
                  0 0 0 0 2 2 3 3
                ];
    end
    
    switch caseNum
        case 1 % simple case
        dimX = 8;%128;
        dimY = 8;%128;
        dimZ = 8;
        case 2
        dimX = 64;
        dimY = 64; 
        dimZ = 64;
        case 3
        dimX = 128;
        dimY = 128;
        dimZ = 128;
    end
    V = zeros(dimY,dimX,dimZ); % the aggregated volume
    
    %%Case 1:
    % simple case
      V(1:8,1:8,1:8) = blk;
%     %%Case 2:
    switch caseNum
        case 2     
%         % harder case
%         dimX = 64;
%         dimY = 64; 
%         dimZ = 64;
%         V = zeros(dimY,dimX,dimZ); % the aggregated volume
        V(9:16,1:8,9:16) = blk;
        V(1:8,9:16,1:8) = 3;
        V(9:16,9:16,9:16) = blk;
        randBlkSize = 16;
          V(10:10+randBlkSize-1,7:7+randBlkSize-1,7:7+randBlkSize-1) = rand(randBlkSize,randBlkSize,randBlkSize);
    %        V(10:10+randBlkSize-1,7:7+blk-1,) = rand(randBlkSize,randBlkSize,randBlkSize);
        case 3 % the case in the OcTree library: generate randomly positioned points
           PTS = rand(100,3);%abs(rand(100,3).*  0.5 .* dimX) ;
           Indx = int16(floor(PTS)+[1 1 1]);
           tId = transpose(1:size(PTS,1));
           for i = 1:length(tId)
             V(Indx(i,1),Indx(i,2),Indx(i,3)) = tId(i);%norm(PTS);
           end
    end
    

    
    if caseNum == 1|| caseNum == 2
      for z = 1:dimZ
          for x = 1:dimX
              for y = 1:dimY
                  if V(x,y,z) > 0
                      PTS(end+1,:) = [x,y,z]; 
                  end
              end
          end
      end
    end
      figure;
      subplot(2,1,1); plot3(PTS(:,1),PTS(:,2),PTS(:,3), 'o');
      subplot(2,1,2); volshow(V);
%     V=rand(32,32);
end