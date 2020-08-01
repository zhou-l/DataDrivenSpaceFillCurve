function [LT, visitOrder] = linearizeScanline3D(I)
    blockDimY = size(I,1);
    blockDimX = size(I,2);
    blockDimZ = size(I,3);
    seqLen= blockDimY * blockDimX*blockDimZ;
    visitOrder = zeros(seqLen,3); %coordinate visiting order 
    LT = zeros(seqLen,1); % the linearized data 
    order = 1;
    for z = 1:blockDimZ
        for y = blockDimY:-1:1
%         for y = 1:blockDimY
            for x = 1:blockDimX
               visitOrder(order,:) = [y,x,z];
               LT(order,:) = I(y,x,z);
               order = order + 1;
            end
        end
    end
end