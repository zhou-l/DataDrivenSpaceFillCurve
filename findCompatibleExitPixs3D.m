
function [exitPixs,isEntryPixInExitFace, entryPixInBlock] = findCompatibleExitPixs3D(imSize, entryPix, exitFace)
    colormat = zeros(imSize(1),imSize(2),imSize(3));
    for i = 1:imSize(1)
        for j = 1:imSize(2)
            for k = 1:imSize(3)
            s = mod(i + j+k,2);
            colormat(i,j,k) = s;
            end
        end
    end
    isEntryPixInExitFace= false;
    exitPixs = [];
    entryPixInBlock = entryPix;
    % even cases
%     if imSize(1) == imSize(2) && 
    if mod(imSize(1) * imSize(2) * imSize(3),2) == 0 % 2*n * 2*n case
        if entryPix(1)>size(colormat,1)
            entryPix(1) = entryPix(1) - size(colormat,1);
        end
        if entryPix(2)>size(colormat,2)
            entryPix(2) = entryPix(2) - size(colormat,2);
        end
        if entryPix(3)>size(colormat,3)
            entryPix(3) = entryPix(3) - size(colormat,3);
        end
        entryPixInBlock = entryPix;
        if entryPix(1)>size(colormat,1) || entryPix(2) > size(colormat,2)||entryPix(3)>size(colormat,3)
            disp('Something went wrong!');
        end
        entryColor = colormat(entryPix(1),entryPix(2),entryPix(3));
        switch(exitFace)
        case 1 % front y = 1 
            j = 1;
           for i = 1:size(colormat,1)
                for k = 1:size(colormat,3)
                    if i == entryPix(1) && j == entryPix(2) && k == entryPix(3)
                        isEntryPixInExitFace = true;
                    end
                    if colormat(i,j,k)~= entryColor
%                         if j == entryPix(1) && i == entryPix(2) && k == entryPix(3)
%                             continue;
%                         else
                            exitPixs(end+1,:) = [i,j,k];
%                         end
                    end
                end
            end  
        case 2 % left x = 1
            i = 1;
            for j = 1:size(colormat,2)
                for k = 1:size(colormat,3)
                    if i == entryPix(1) && j == entryPix(2) && k == entryPix(3)
                        isEntryPixInExitFace = true;
                    end
                    if colormat(i,j,k)~= entryColor
%                         if j == entryPix(1) && i == entryPix(2) && k == entryPix(3)
%                             continue;
%                         else
                            exitPixs(end+1,:) = [i,j,k];
%                         end
                    end
                end
            end
        case 3 % top z = imSize(3)
            k = imSize(3);
            for i = 1:size(colormat,1)
                for j = 1:size(colormat,2)
                    if i == entryPix(1) && j == entryPix(2) && k == entryPix(3)
                        isEntryPixInExitFace = true;
                    end

                    if colormat(i,j,k)~= entryColor
%                         if j == entryPix(1) && i == entryPix(2) && k == entryPix(3)
%                             continue;
%                         else
                            exitPixs(end+1,:) = [i,j,k];
%                         end
                    end
                end
            end
             
         
        case 4 % right x = imSize(2)
            i = imSize(1);
            for j = 1:size(colormat,2)
                for k = 1:size(colormat,3)
                    if i == entryPix(1) && j == entryPix(2) && k == entryPix(3)
                        isEntryPixInExitFace = true;
                    end
                    if colormat(i,j,k)~= entryColor
%                         if j == entryPix(1) && i == entryPix(2) && k == entryPix(3)
%                             continue;
%                         else
                            exitPixs(end+1,:) = [i,j,k];
%                         end
                    end
                end
            end
        case 5 % bottom z = 1
           k = 1;
            for j = 1:size(colormat,2)
                for i = 1:size(colormat,1) % scan the front face
                    if i == entryPix(1) && j == entryPix(2) && k == entryPix(3)
                        isEntryPixInExitFace = true;
                    end
                     if colormat(i,j,k)~= entryColor   
%                         if j == entryPix(1) && i == entryPix(2) && k == entryPix(3)
%                             continue;
%                         else
                            exitPixs(end+1,:) = [i,j,k];
%                         end
                     end
                end
            end
        case 6 % back y = imSize(1)
              j = imSize(2);
             for i = 1:size(colormat,1)
                for k = 1:size(colormat,3)
                    if i == entryPix(1) && j == entryPix(2) && k == entryPix(3)
                        isEntryPixInExitFace = true;
                    end
                     if colormat(i,j,k)~= entryColor
%                         if j == entryPix(1) && i == entryPix(2) && k == entryPix(3)
%                             continue;
%                         else
                            exitPixs(end+1,:) = [i,j,k];
%                         end
                    end
                end
             end
        end   
    else
        return;
    end
end