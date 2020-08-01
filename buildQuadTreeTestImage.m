function [V,dimX,dimY] = buildQuadTreeTestImage()
    dimX = 32;%64;%128;
    dimY = 32;%64;%128;
    dimX = 8; dimY = 8;
     V = zeros(dimY,dimX);
  
   
  blk   = [ 0 0 4 3 0 0 0 0;
          0 0 1 2 0 0 0 0;
          0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0;
          0 0 0 0 3 1 4 4;
          0 0 0 0 1 2 4 4;
          0 0 0 0 2 2 3 3;
          0 0 0 0 2 2 3 3
        ];
      blk   = [ 0 0 4 3 3 1 4 4;
          0 0 1 2 1 2 4 4;
          0 0 0 0 2 2 3 3;
          0 0 0 0 1 2 3 4;
          0 0 0 0 3 1 4 4;
          0 0 0 0 1 2 4 4;
          0 0 0 0 2 2 3 3;
          0 0 0 0 2 2 3 3
        ];
      V(1:8,1:8,:) = blk;
%     V(9:16,1:8,:) = blk;
%     V(1:8,9:16,:) = 3;
%     V(9:16,9:16,:) = blk;
    randBlkSize = 32;
%       V(10:10+randBlkSize-1,7:7+randBlkSize-1,:) = rand(randBlkSize,randBlkSize);
%           V(10:10+blk-1,7:7+blk-1,:) = rand(randBlkSize,randBlkSize);
%     V=rand(32,32);
end