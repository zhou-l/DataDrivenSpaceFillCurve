function [reducedData, hilbertInd] = hilbertCurve2D(data)
%hilbertCurve For now only 2D to 1D
%   data = the 2D matrix. Needs to be square and each dimension needs to be
%   2^n length.
% Most of algorithm from: https://blogs.mathworks.com/steve/2012/01/25/generating-hilbert-curves/
rowNumel = sqrt(numel(data));
order = log(numel(data))/log(4);
% convert original index to x,y format
a = 1 + 1i;
b = 1 - 1i;
% Generate point sequence
z = 0;
for k = 1:order
    w = 1i*conj(z);
    z = [w-a; z-b; z+a; b-w]/2;
end
newCol = real(z);
newRow = imag(z);
newCol = rowNumel*newCol/2 + rowNumel/2 + 0.5;
newRow = rowNumel*newRow/2 + rowNumel/2 + 0.5;
hilbertInd = (newCol-1)*rowNumel+newRow;
reducedData = data(hilbertInd);
% reducedData = imresize(reducedData, 0.25);
% reducedData(reducedData<0) = 0;