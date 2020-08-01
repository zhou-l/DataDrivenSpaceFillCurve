function [avgOursOctTreeACL, avgOursRegGridACL, avgHilbertACL, avgRegGridLineACL,...
     avgDistOursRegGridACL, avgDistHilbertACL, avgDistScanlineACL] = compareRegGridandOctTree(filename)
    [clLT, clVisitOrder, fullLT] = SFCOctTreeMultiScaleMain(filename);
    isPad = true;
    [myLT, myVO, HLT, HVO, lLT, lVO, dimX, dimY, dimZ] = figGenSFC3DCase(filename);
    % compare quadtree-reconstructed autocorr and regular grid autocorr
    lagCnt = 500;
    %% quad tree - ours
    avgOursOctTreeACL = compAvgAutoCorr(fullLT, lagCnt);
    %% regular grid - ours
    centraledIds = -length(fullLT)/2:length(fullLT)/2-1;
    centraledIds = abs(centraledIds);
    centraledIds= transpose(centraledIds);
    idOrder = centraledIds;
    avgOursRegGridACL = compAvgAutoCorr(myLT, lagCnt);
    dist = vecnorm(myVO(:,:) - [1 1 1], 2, 2);
    avgDistOursRegGridACL = compAvgAutoCorr(dist, lagCnt);
%     idOrder = vecnorm(myVO(:,:) - [dimY/2, dimX/2, dimZ/2],2,2);
%     idOrder = idOrder - centraledIds;
%     avgIdOursRegGridACL = compAvgAutoCorr(idOrder,lagCnt);
%     %% regular grid - context-based
%     avgCtRegGridACL = compAvgAutoCorr(ctLT, lagCnt);
%     dist = vecnorm(ctVO(:,:) - [1 1 1], 2, 2);
%     avgDistCtRegGridACL = compAvgAutoCorr(dist, lagCnt);
% %     idOrder = vecnorm(ctVO(:,:) - [dimX/2, dimY/2],2,2); %(ctVO(:,1)-1) * dimX + ctVO(:,2);
% %     idOrder = idOrder - centraledIds;
%     avgIdCtRegGridACL = compAvgAutoCorr(idOrder,lagCnt);
    %% scanline
    avgRegGridLineACL = compAvgAutoCorr(lLT, lagCnt);
    dist = vecnorm(lVO(:,:) - [1 1 1], 2, 2);
    avgDistScanlineACL = compAvgAutoCorr(dist, lagCnt);
%     idOrder = vecnorm(lVO(:,:) - [dimX/2, dimY/2],2,2);%(lVO(:,1)-1) * dimX + lVO(:,2);
% %     idOrder = idOrder - centraledIds;
%     avgIdScanlineACL= compAvgAutoCorr(idOrder,lagCnt);
    %% hilbert
    if~isempty(HLT)
        avgHilbertACL = compAvgAutoCorr(HLT, lagCnt);
        dist = vecnorm(HVO(:,:) - [1 1 1], 2, 2);
        avgDistHilbertACL = compAvgAutoCorr(dist, lagCnt);
%         idOrder =  vecnorm(HVO(:,:) - [dimX/2, dimY/2],2,2); %(HVO(:,1)-1) * dimX + HVO(:,2);
% %         idOrder = idOrder - centraledIds;
%         avgIdHilbertACL  = compAvgAutoCorr(idOrder,lagCnt);
    end
%     plot
%     plot avg corrs
%     figure; hold on;
%     plot(1:lagCnt, avgIdOursRegGridACL, 'ko-');
% %     plot(1:lagCnt, avgIdCtRegGridACL, 'go-');
%     if ~isempty(avgIdHilbertACL)
%         plot(1:lagCnt, avgIdHilbertACL, 'bo-');
%     end
%     
%     plot(1:lagCnt, avgIdScanlineACL, 'ro-');
%     legend('RegGrid-Ours',, 'Hilbert', 'Scanline');
%     title('AutoCorr QuadTree v.s. RegGrid Dist');
%     hold off;
% end