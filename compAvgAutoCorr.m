function avgACL = compAvgAutoCorr(LT, lagCnt)   
% compute averaged autocorrelation
   avgACL = zeros(lagCnt,1);
   for i = 1:lagCnt
%        [c, lags] = xcorr(LT, LT, i);
       [c, lags] = xcorr(LT, LT, i, 'normalized');
        avgAC = mean(c);
        avgACL(i,:) = avgAC;
%        stem(lags, c);
   end
end