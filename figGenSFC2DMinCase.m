function [ctLT, ctVO, myLT, myVO, HLT, HVO, lineLT, lineVO] = figGenSFC2DMinCase
global V;
 zisPlot = true;
V = [0 0 0 0 0 0 0 0
     0 1 1 1 0 0 0 0
     0 1 1 1 1 1 0 0
     0 0 0 0 1 1 0 0];
 isPad = false;
baseFileName = 'minimum2Dcase';
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
if zisPlot
    subplot(1,2,2); imagesc(V);
end

%% All SFCs
[ctLT, ctVO, myLT, myVO, HLT, HVO, lineLT, lineVO] = SFCMethodsTestMain(V, isPad);
LTfilename = sprintf('LT%s.csv', baseFileName);
VOfilename = sprintf('VO%s.csv', baseFileName);
% 
% % write out SFC
csvwrite(LTfilename, myLT);
csvwrite(VOfilename, myVO);

if isPad
    hLTfilename = sprintf('HilbertLT%s.csv', baseFileName);
    hVOfilename = sprintf('HilbertVO%s.csv', baseFileName);
    csvwrite(hLTfilename, HLT);
    csvwrite(hVOfilename, HVO); 
end
end