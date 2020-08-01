function I = plotDisks(dimI, rawCtrs, R)
dimY = dimI(1);
dimX = dimI(2);
Ctrs = round(rawCtrs);
I = zeros(dimY, dimX);

for i = 1:length(Ctrs)
    for y = Ctrs(i,1) - R:Ctrs(i,1)+R
        for x = Ctrs(i,2)-R:Ctrs(i,2)+R
            if x < 1 || y < 1 || x > dimX || y > dimY
                continue;
            end
            dist = norm([y,x] - Ctrs(i,:));
            if dist <= R 
                I(y,x) = 255;
%                if i == 2
%                     I(y,x) = -200/R * dist+ 255;
%                else
%                I(y,x) = -255/R * dist+ 255;
%                end
            end
       end
    end
end
