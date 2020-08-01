function evalAutoCorrOctTree
testFiles ={
'data/testVols/fuel32.nrrd',...
'data/testVols/heart32.nrrd',...
'data/testVols/neghip32.nrrd',...
'data/testVols/nucleon32.nrrd',...
'data/testVols/tangle32.nrrd'};

avgOursQuadTreeACL = [];
avgOursRegGridACL = [];
avgContextRegGridACL = [];
avgHilbertACL = [];
avgRegGridLineACL = [];
avgDistOursRegGridACL = [];
avgDistContextRegGridACL = [];
avgDistHilbertACL = [];
avgDistScanlineACL = [];
avgIdOursRegGridACL = [];
avgIdCtRegGridACL = [];
avgIdHilbertACL = [];
avgIdScanlineACL = [];
L = length(testFiles);
for i = 1:L
    [a1, a2, a3, a4, d1, d2, d3] = compareRegGridandOctTree(testFiles{i});
    if i == 1
        avgOursQuadTreeACL = a1;
        avgOursRegGridACL = a2;
        avgHilbertACL = a3;
         avgRegGridLineACL = a4;
        avgDistOursRegGridACL = d1;
        avgDistHilbertACL = d2;
        avgDistScanlineACL = d3;

    else
        avgOursQuadTreeACL = avgOursQuadTreeACL + a1;
        avgOursRegGridACL = avgOursRegGridACL + a2;
        avgHilbertACL = avgHilbertACL + a3;
        avgRegGridLineACL = avgRegGridLineACL + a4;
        avgDistOursRegGridACL = avgDistOursRegGridACL + d1;
        avgDistHilbertACL = avgDistHilbertACL + d2;
        avgDistScanlineACL = avgDistScanlineACL + d3;
    end
    close all;
end
avgOursQuadTreeACL = avgOursQuadTreeACL ./ L;
avgOursRegGridACL = avgOursRegGridACL ./ L;
avgHilbertACL = avgHilbertACL ./ L;
avgRegGridLineACL = avgRegGridLineACL ./ L;
avgDistOursRegGridACL = avgDistOursRegGridACL ./ L;
avgDistHilbertACL = avgDistHilbertACL ./ L;
avgDistScanlineACL = avgDistScanlineACL ./ L;


lagCnt = length(avgOursQuadTreeACL);
%% plot autocorr of values
figure; hold on;
plot(1:lagCnt, avgOursQuadTreeACL, 'k-', 'Linewidth', 5);
plot(1:lagCnt, avgOursRegGridACL, 'color',[0 0.5 0], 'linestyle', '-','Linewidth', 5);
if ~isempty(avgHilbertACL)
    plot(1:lagCnt, avgHilbertACL, 'b-','Linewidth', 5);
end

plot(1:lagCnt, avgRegGridLineACL, 'r-','Linewidth', 5);
legend('Ours-Quadtree', 'Ours-RegGrid', 'Hilbert', 'Scanline');
title('Averaged AutoCorr OctTree v.s. RegGrid Value');
hold off;
%% plot autocorr of Eucledian distance
figure; hold on;
% plot(1:lagCnt, avgDistOursRegGridACL, 'ko-');
plot(1:lagCnt, avgDistOursRegGridACL, 'color',[0 0.5 0], 'linestyle', '-','Linewidth', 5);
if ~isempty(avgHilbertACL)
    plot(1:lagCnt, avgDistHilbertACL, 'b-','Linewidth', 5);
end
plot(1:lagCnt, avgDistScanlineACL, 'r-','Linewidth', 5);
legend('Ours-RegGrid', 'Hilbert', 'Scanline');
title('Averaged AutoCorr RegGrid Distance');
hold off;
