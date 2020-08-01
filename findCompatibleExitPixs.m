
function exitPixs = findCompatibleExitPixs(imSize, entryPix, exitEdge)
    colormat = zeros(imSize(1),imSize(2));
    for i = 1:imSize(1)
        for j = 1:imSize(2)
            s = mod(i + j,2);
            colormat(i,j) = s;
        end
    end
    exitPixs = [];
    % even cases
%     if imSize(1) == imSize(2) && 
    if mod(imSize(1) * imSize(2),2) == 0 % 2*n * 2*n case
        if entryPix(1)>size(colormat,1)
            entryPix(1) = entryPix(1) - size(colormat,1);
        end
        if entryPix(2)>size(colormat,2)
            entryPix(2) = entryPix(2) - size(colormat,2);
        end
        entryColor = colormat(entryPix(1),entryPix(2));
        cnt = 1;
        switch(exitEdge)
        case 1 % bottom
            exitPixs = zeros(imSize(2)/2, 2);% one entry point can have two 
            exitPixs(:,1) = imSize(1); 
            
            for i = 1:size(colormat,2)
                if colormat(exitPixs(cnt,1),i)~= entryColor
                    exitPixs(cnt,2) = i;
                    cnt = cnt + 1;
                    if cnt >= size(exitPixs,1)+1
                        return;
                    end
                end
            end
        case 2 % right
            exitPixs = zeros(imSize(1)/2, 2);% one entry point can have two 
            exitPixs(:,2) = imSize(2);
            for i = 1:size(colormat,1)
                if colormat(i,exitPixs(cnt,2))~= entryColor
                    exitPixs(cnt,1) = i;
                    cnt = cnt + 1;
                    if cnt >= size(exitPixs,1)+1
                        return;
                    end
                end
            end
        case 3 % top
            exitPixs = zeros(imSize(2)/2, 2);% one entry point can have two 
           exitPixs(:,1) = 1;
           for i = 1:size(colormat,2)
                if colormat(exitPixs(cnt,1),i)~= entryColor
                    exitPixs(cnt,2) = i;
                    cnt = cnt + 1;
                    if cnt >= size(exitPixs,1)+1
                        return;
                    end
                end
            end
        case 4 % left
            exitPixs = zeros(imSize(1)/2, 2);% one entry point can have two 
            exitPixs(:,2) = 1;
            for i = 1:size(colormat,1)
             if colormat(i,exitPixs(cnt,2))~= entryColor
                    exitPixs(cnt,1) = i;
                    cnt = cnt + 1;
                    if cnt >= size(exitPixs,1)+1
                        return;
                    end
                end
            end   
        end   
    else
        return;
    end
end