function Vp = padToNextPow2(V)
% use sphere test case
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
    Vp = padarray(V, [padDimFirstHalf(1),padDimFirstHalf(2)], 'pre');
    Vp = padarray(V, [padDimSecondHalf(1),padDimSecondHalf(2)],'post');
end