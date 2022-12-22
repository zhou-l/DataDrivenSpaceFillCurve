
function AdjMat = createAdjMatrixVol(V, sizeV)
w = sizeV(1);
h = sizeV(3);
d = sizeV(2);
%     h = size(V,1); w = size(V,2); d = size(V,3);
    numVox = w*h*d;
    AdjMat = zeros(numVox, numVox);
   
    % TODO: ensure the correct indexing using Matlab's convention!!!
    for z = 1:h
        for y = 1:d
            for x = 1:w
                id = ((z-1)*d+(y-1))*w+x-1 + 1;
                
                % check neighborhood
                for nz = -1:1
                    for ny = -1:1
                        for nx = -1:1
                            distToCtr = norm([nx,ny,nz],2);
                            if distToCtr > 1 || distToCtr == 0% skip if L2-norm is greater than 1 or equal to 0
                                continue;
                            end
                            zz = z + nz;
                            yy = y + ny;
                            xx = x + nx;
                            if zz < 1 || yy < 1 || xx < 1 ||...
                                    zz > h || yy > d || xx > w
                                continue; % out of bound
                            end
                            
                            idn = ((zz-1)*d+(yy-1))*w+xx-1 + 1;
                            
                            % set connection
                            AdjMat(id,idn) = 1;
                        end
                    end
                end
            end
        end
    end
end




