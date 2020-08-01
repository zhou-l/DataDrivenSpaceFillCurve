function genQuadTreeV
dimX = 64;
dimY = 64;
dimZ = 1;
nSphere = 5;
V = testVolCreate(dimX, dimY, dimZ, nSphere);
S = qtdecomp(V);
blocks = repmat(uint8(0),size(S));
for dim = [512 256 128 64 32 16 8 4 2 1];    
  numblocks = length(find(S==dim));    
  if (numblocks > 0)        
    values = repmat(uint8(1),[dim dim numblocks]);
    values(2:dim,2:dim,:) = 0;
    blocks = qtsetblk(blocks,S,dim,values);
  end
end


blocks(end,1:end) = 1;
blocks(1:end,end) = 1;

imshow(V);
figure
imshow(blocks,[]);