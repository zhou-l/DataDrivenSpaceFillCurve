function evalAutoCorrQuadTree
testFiles ={
'aneurism64.png',...
'beattle64.png',...
'bonsai64.png',...
'brainBruce64.png',...
'disk5.png',...
'engine_64.png',...
'foot64.png',...
'fuel.png',...
'isabelSmSm.png',...
'neghip.png',...
'nucleon64.png'};

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
    close all;
    [a1, a2, a3, a4, a5, d1, d2, d3, d4, o1, o2, o3, o4] = compareRegGridandQuadTree(testFiles{i});
    if i == 1
        avgOursQuadTreeACL = a1;
        avgOursRegGridACL = a2;
        avgContextRegGridACL = a3;
        avgHilbertACL = a4;
        avgRegGridLineACL = a5;
        avgDistOursRegGridACL = d1;
        avgDistContextRegGridACL = d2;
        avgDistHilbertACL = d3;
        avgDistScanlineACL = d4;
         avgIdOursRegGridACL = o1;
         avgIdCtRegGridACL = o2;
         avgIdHilbertACL = o3;
         avgIdScanlineACL = o4;
    else
        avgOursQuadTreeACL = avgOursQuadTreeACL + a1;
        avgOursRegGridACL = avgOursRegGridACL + a2;
        avgContextRegGridACL = avgContextRegGridACL + a3;
        avgHilbertACL = avgHilbertACL + a4;
        avgRegGridLineACL = avgRegGridLineACL + a5;
        avgDistOursRegGridACL = avgDistOursRegGridACL + d1;
        avgDistContextRegGridACL = avgDistContextRegGridACL + d2;
        avgDistHilbertACL = avgDistHilbertACL + d3;
        avgDistScanlineACL = avgDistScanlineACL + d4;
            avgIdOursRegGridACL = avgIdOursRegGridACL + o1;
         avgIdCtRegGridACL = avgIdCtRegGridACL + o2;
         avgIdHilbertACL = avgIdHilbertACL + o3;
         avgIdScanlineACL = avgIdScanlineACL + o4;
    end
end
avgOursQuadTreeACL = avgOursQuadTreeACL ./ L;
avgOursRegGridACL = avgOursRegGridACL ./ L;
avgContextRegGridACL = avgContextRegGridACL ./ L;
avgHilbertACL = avgHilbertACL ./ L;
avgRegGridLineACL = avgRegGridLineACL ./ L;
avgDistOursRegGridACL = avgDistOursRegGridACL ./ L;
avgDistHilbertACL = avgDistHilbertACL ./ L;
avgDistScanlineACL = avgDistScanlineACL ./ L;
avgDistContextRegGridACL= avgDistContextRegGridACL ./ L;

%% calculates aggregated values
sumOursQTACL = sum(avgOursQuadTreeACL);
sumOursRGACL = sum(avgOursRegGridACL);
sumContextRGACL = sum(avgContextRegGridACL);
sumHilbertACL = sum(avgHilbertACL);
sumRGLineACL = sum(avgRegGridLineACL);
sumDistOursRGACL = sum(avgDistOursRegGridACL); 
sumDistHilbertACL = sum(avgDistHilbertACL);
sumDistScanlineACL = sum(avgDistScanlineACL);
sumDistContextRGACL = sum(avgDistContextRegGridACL); 

% A = [sumOursQTACL, sumOursRGACL, sumContextRGACL, sumHilbertACL, sumRGLineACL, sumDistOursRGACL, sumDistHilbertACL, sumDistScanlineACL, sumDistContextRGACL]; 
% T = array2table(A,'VariableNames',{'sumOursQTACL','sumOursRGACL','sumContextRGACL', 'sumHilbertACL', 'sumRGLineACL', 'sumDistOursRGACL', 'sumDistOursRGACL',...
%     'sumDistHilbertACL', 'sumDistScanlineACL', 'sumDistContextRGACL'});
% writetable(T, 'aggAutoCorrQuadTree.csv');

%%
avgIdOursRegGridACL = avgIdOursRegGridACL ./ L;
avgIdCtRegGridACL = avgIdCtRegGridACL ./ L;
avgIdHilbertACL = avgIdHilbertACL ./ L;
avgIdScanlineACL = avgIdScanlineACL ./ L;
lagCnt = length(avgOursQuadTreeACL);
%% plot autocorr of values
figure; hold on;
plot(1:lagCnt, avgOursQuadTreeACL, 'k-', 'LineWidth', 4);
plot(1:lagCnt, avgOursRegGridACL, 'color',[0 0.5 0], 'linestyle', '-', 'LineWidth', 4);
plot(1:lagCnt, avgContextRegGridACL, 'm-', 'LineWidth', 4);
if ~isempty(avgHilbertACL)
    plot(1:lagCnt, avgHilbertACL, 'b-', 'LineWidth', 4);
end

plot(1:lagCnt, avgRegGridLineACL, 'r-', 'LineWidth', 4);
legend('Ours-Quadtree', 'Ours-RegGrid', 'Context-RegGrid', 'Hilbert', 'Scanline');
title('Averaged AutoCorr QuadTree v.s. RegGrid Value');
hold off;
%% plot autocorr of Eucledian distance
figure; hold on;
% plot(1:lagCnt, avgDistOursRegGridACL, 'ko-');
plot(1:lagCnt, avgDistOursRegGridACL, 'color',[0 0.5 0], 'linestyle', '-', 'LineWidth', 4);
plot(1:lagCnt, avgDistContextRegGridACL, 'm-', 'LineWidth', 4);
if ~isempty(avgHilbertACL)
    plot(1:lagCnt, avgDistHilbertACL, 'b-', 'LineWidth', 4);
end
plot(1:lagCnt, avgDistScanlineACL, 'r-', 'LineWidth', 4);
legend('Ours-RegGrid', 'Context-RegGrid', 'Hilbert', 'Scanline');
title('Averaged AutoCorr RegGrid Distance');
hold off;
% %% plot autocorr of id ordering
% figure; hold on;
% % plot(1:lagCnt, avgDistOursRegGridACL, 'ko-');
% plot(1:lagCnt, avgIdOursRegGridACL, 'go-');
% plot(1:lagCnt, avgIdCtRegGridACL, 'mo-');
% if ~isempty(avgHilbertACL)
%     plot(1:lagCnt, avgIdHilbertACL, 'bo-');
% end
% plot(1:lagCnt, avgIdScanlineACL, 'ro-');
% legend('Ours-RegGrid', 'Context-RegGrid', 'Hilbert', 'Scanline');
% title('Averaged AutoCorr RegGrid Id Ordering');
% hold off;