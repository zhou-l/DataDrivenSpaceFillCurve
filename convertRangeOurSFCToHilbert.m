function rangeHilbert = convertRangeOurSFCToHilbert(rangeOurs, VOours, VOhilbert, minRange)
    rangeHilbert = [];
    if nargin < 4
        minRange = 20;
    end
for jj = 1:size(rangeOurs,1)
    Rr = rangeOurs(jj,:);% [20248.4, 21449.4];
    PosOurs = VOours(int16(Rr(1)):int16(Rr(2)), :);
    % search the pos in VOhilbert

    Hidx = zeros(length(PosOurs),1);
    for i = 1:length(PosOurs)
        
        idX = find(VOhilbert(:,1) == PosOurs(i,1) & VOhilbert(:,2) == PosOurs(i,2) & VOhilbert(:,3) == PosOurs(i,3));
        Hidx(i,:) = idX;
        
    end
    Hidx = sort(Hidx);
    for i = 1:length(Hidx)
        idX = Hidx(i);
        if i == 1
            Rmin = idX;
            Rmax = idX;
            
            rangeHilbert(end+1,:) = [Rmin, Rmax];
        else
            if idX - Rmax > 5 % 20 % start a new range
                    Rmin = idX;
                    Rmax = idX;
%                 if rangeHilbert(end,2) - rangeHilbert(end,1) > minRange
%                     rangeHilbert(end+1,:) = [Rmin, Rmax];
%                 else
%                     rangeHilbert(end,:) = [Rmin, Rmax];
%                 end
                 rangeHilbert(end+1,:) = [Rmin, Rmax];
            else
                Rmax = idX;
                rangeHilbert(end,2) = Rmax; % update rmax
            end
        end
    end
end

xfRangeHilbert = rangeHilbert;
    toRmInd = rangeHilbert(:,2) - rangeHilbert(:,1) < minRange;
    xfRangeHilbert(toRmInd,:) = [];
rangeHilbert = xfRangeHilbert;

end