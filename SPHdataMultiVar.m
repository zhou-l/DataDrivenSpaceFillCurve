function SPHdataMultiVar(SPHdataFile, SPHvoFile, SPHltFile)
close all;
VO = dlmread(SPHvoFile);
LT = dlmread(SPHltFile);
dim = max(VO,[],1);
dimX = max(dim);
dimY = dimX;
dimZ = dimX;
% load particle data
    X = readtable(SPHdataFile);
    PTS = [X.posX, X.posY, X.posZ];
                % original code:
    PTminmax = [min(PTS,[],1) max(PTS,[],1)];   
    PTr = PTminmax(4:6) - PTminmax(1:3);
    PTS = (PTS - PTminmax(1:3)) .* [dimX-1 dimY-1 dimZ-1] ./ PTr + [1 1 1];
    VP = zeros(dimX, dimY, dimZ);
    VD = zeros(dimX, dimY, dimZ);
    VV = zeros(dimX, dimY, dimZ);
    % transform to [0,dim]
    plot3(PTS(:,1),PTS(:,3),PTS(:,2), 'ro');
   Indx = int16(floor(PTS));
   tId = transpose(1:size(PTS,1));
   for i = 1:length(tId)
     VP(Indx(i,1),Indx(i,2),Indx(i,3)) = X.pressure(i);
     VD(Indx(i,1),Indx(i,2),Indx(i,3)) = X.density(i);
     VV(Indx(i,1),Indx(i,2),Indx(i,3)) = sqrt(X.velocityX(i).^2 + X.velocityY(i).^2 + X.velocityZ(i).^2);
   end
   % plot
   figure; hold on;
   mvarLT = zeros(length(VO),3);
   for i = 1:length(VO)
       pt = VO(i,:);
       vp = VP(pt(1),pt(2),pt(3));
       if VP(pt(1),pt(2),pt(3)) > 0 
       vp = VP(pt(1),pt(2),pt(3)); % / 1000;
       end
       
       vd = VD(pt(1),pt(2),pt(3));
%        if VD(pt(1),pt(2),pt(3)) > 0 
%            vd = VD(pt(1),pt(2),pt(3)) - 1000;
%        end
       vv = VV(pt(1),pt(2),pt(3));
       mvarLT(i,:) = [vp,vd,vv];
   end
   
   
   
   dlmwrite('LT_db4k_pressure.csv',mvarLT(:,1));
   dlmwrite('LT_db4k_density.csv', mvarLT(:,2));
   dlmwrite('LT_db4k_velocity.csv', mvarLT(:,3));
   
      fullLTp = fullSFC([mvarLT(:,1),LT(:,2)]);
     subplot(3,1,1); plot(1:length(fullLTp),(fullLTp+1));
    fullLTd = fullSFC([mvarLT(:,2),LT(:,2)]);
   subplot(3,1,2); plot(1:length(fullLTd),(fullLTd));ylim([990 1020]);
    fullLTv = fullSFC([mvarLT(:,3),LT(:,2)]);
   subplot(3,1,3); plot(1:length(fullLTv),fullLTv);
   hold off;
   
   
   dlmwrite('LT_db4k_pressure.csv',VP);
   dlmwrite('LT_db4k_density.csv', VD);
   dlmwrite('LT_db4k_velocity.csv', VV);
   
   % compute averaged autocorrelation
   avgACL = zeros(10,1);
   for i = 1:10
       [c, lags] = xcorr(fullLTp, fullLTp, i, 'normalized');
        avgAC = mean(c);
        avgACL(i,:) = avgAC;
%        stem(lags, c);
   end
   figure, stem(1:10, avgACL);
end