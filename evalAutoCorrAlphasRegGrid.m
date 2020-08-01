function [cummulateACvals, cummulateACdist] = evalAutoCorrAlphasRegGrid
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


L = length(testFiles);
avgACvals = [];
avgACdist = [];
numAlpha = 1;
len = 1;
for i = 1:L
    close all;
    [listACvals, listACdist]= figRegGrid2DalphaEffect(testFiles{i});
  
    if i == 1
        [numAlpha, len] = size(listACvals);
        avgACvals = listACvals;
        avgACdist = listACdist;
    else
        avgACvals = avgACvals + listACvals;
        avgACdist = avgACdist + listACdist;
    end
end
figure; 
title('Autocorr of data val');
hold on;
for i = 1:numAlpha
    plot(1:len, avgACvals(i,:));
end
legend('\alpha = 0', '\alpha = 0.2', '\alpha = 0.4', '\alpha = 0.6', '\alpha = 0.8', '\alpha = 1.0');
hold off;

figure; 
title('Autocorr of dist');
hold on;
for i = 1:numAlpha
    plot(1:len, avgACdist(i,:));
end
legend('\alpha = 0', '\alpha = 0.2', '\alpha = 0.4', '\alpha = 0.6', '\alpha = 0.8', '\alpha = 1.0');
hold off;

cummulateACvals = sum(avgACvals, 2);
cummulateACdist = sum(avgACdist, 2);

avgCumACvals = cummulateACvals ./ L;
avgCumACdist = cummulateACdist ./ L;

