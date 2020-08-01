close all;
V = testVolCreate(4, 4, 4, 1);
entryPix = [3 1 1];
exitFace = 5;
[LT, visitOrder, exitPix] = linearizeHamPath3D(V, entryPix, exitFace);
id = 1:length(LT);
lineColorCoded(visitOrder(:,2)', visitOrder(:,1)', visitOrder(:,3)', id);