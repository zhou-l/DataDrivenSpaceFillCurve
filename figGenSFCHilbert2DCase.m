function figGenSFCHilbert2DCase(filename, isPad)
global V;
V = imread(filename);
dV = size(V);
if length(dV)==3
    V = double(rgb2gray(V));
else
    V = double(V(:,:,1));
end
  V = V ./ 255;
figure, hold on;
subplot(1,2,1); imagesc(V);
  V = imdiffusefilt(V);
V = V .* 255;

if isPad
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
    V = padarray(V, [padDimFirstHalf(1),padDimFirstHalf(2)], 'pre');
    V = padarray(V, [padDimSecondHalf(1),padDimSecondHalf(2)],'post');
end
subplot(1,2,2); imagesc(V);
[folder, baseFileName, ext] = fileparts(filename);
% %% Hilbert SFC:

    dimX = size(V,2);
    dimY = size(V,1);
    dimZ = size(V,3);
 [rD, rI] = hilbertCurve2D(V);
    subplot(2,2,4);plot(1:length(rD),rD);title('Hilbert');
    HLT = rD;
    % get 2D coords from rI
    HVO = zeros(length(rI),2);
    for i = 1:length(rI)
        HVO(i,1) = floor(rI(i)/dimX) + 1;
        HVO(i,2) = mod(rI(i)-1, dimX)+1;
    end
    hLTfilename = sprintf('HilbertLT%s.csv', baseFileName);
    hVOfilename = sprintf('HilbertVO%s.csv', baseFileName);
    csvwrite(hLTfilename, HLT);
    csvwrite(hVOfilename, HVO);
return;