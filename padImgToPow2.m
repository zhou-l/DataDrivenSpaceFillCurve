function pV = padImgToPow2(V)
    maxSize = max(size(V));
    nextPow2 = 2^ceil(log2(maxSize));
    padDimFirstHalf = zeros(2,1);
    padDimSecondHalf = zeros(2,1);
    for i = 1:2
        if mod(size(V,i),2) == 0
            padDimFirstHalf(i) = (nextPow2 - size(V,i))/2;
            padDimSecondHalf(i) = (nextPow2 - size(V,i))/2;
        else
            padDimFirstHalf(i) = nextPow2/2 - floor(size(V,i)/2);
            padDimSecondHalf(i) = nextPow2 - size(V,i) - padDimFirstHalf(i);
        end
    end
    pV = padarray(V, [padDimFirstHalf(1),padDimFirstHalf(2)], 'pre');
    pV = padarray(pV, [padDimSecondHalf(1),padDimSecondHalf(2)],'post');
end