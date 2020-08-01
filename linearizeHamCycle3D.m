%Traversing the image I using a Hamilton cycle generated from a precomputed minimum spanning tree
% use tree as input
function [zLT,zVisitOrder] = linearizeHamCycle3D(mstT, I, entryPix)
    currPix = entryPix;
    currCirc = [floor((currPix(1)+1)/2),floor((currPix(2)+1)/2),floor((currPix(3)+1)/2)];
    lastCirc = [0 0 0];
    circVisitCnt = zeros(size(mstT,1),size(mstT,2),size(mstT,3));
    global pixOcc;%pixel occupancy map
    global LT;
    global visitOrder;
    global pixId;
    global nneighbors;
    global neighborDirs;
    LT = zeros(size(I,1)*size(I,2)*size(I,3),2);
    visitOrder = zeros(size(I,1)*size(I,2)*size(I,3),3);
    pixOcc = false(size(I,1),size(I,2),size(I,3),1);
    pixId = 1;
    nneighbors = 6;
    neighborDirs = [
        -1,0,0; % up
        1,0,0; % down
        0,1,0; %right
        0,-1,0;%left
        0,0,1;%back
        0,0,-1%front
        ];
%     dirSeq = genDirSeq(mstSet);
    fromDir = [0,0,0];
    toDir = fromDir;
    stack = [currCirc]; 

    while ~isempty(stack)
        currCirc = stack(1,1:3);   
        numvisited = circVisitCnt(currCirc(1),currCirc(2),currCirc(3));
        % remove first
        stack(1,:) = [];
        fromDir = toDir;
        isCurrCircDone = false;
        children = mstT(currCirc(1),currCirc(2),currCirc(3),:);
        children = reshape(children, [nneighbors,1]);
        childCnt = nnz(children);

        if numvisited < childCnt % not all paths are visited
            for i = 1:nneighbors
                if children(i) > 0
                       % for the current channel==direction version
                       tryDir = neighborDirs(i,:);
                       tryPix = currPix + tryDir;
                       if pixOcc(tryPix(1),tryPix(2),tryPix(3))
                           continue;
                       end

                       if checkDirValid3D(tryDir, currPix, currCirc,...
                               lastCirc, children) %check which direction can go
                         nextCirc = currCirc + tryDir;
                          % if next pixel is not in circuit, we visit
                          % current pixel and skip all
                         trynextPix = currPix+tryDir;
                         inCirc = checkPixInCirc3D(trynextPix, currCirc);
                         if ~inCirc
                            currPix=visitPix3D(currPix,currCirc,[0,0,0],I);
                            lastCirc = currCirc;
                            isCurrCircDone = true;      
                            currPix = trynextPix;
                         end
                          toDir = tryDir;
                          if pixId == 1
                              fromDir = toDir;
                          end
                          stack = [nextCirc; currCirc;stack];
                          circVisitCnt(currCirc(1),currCirc(2),currCirc(3)) =...
                              circVisitCnt(currCirc(1),currCirc(2),currCirc(3)) + 1;
                          break;
                       end
                end
            end
        else % has been visited
        % trace back 
%             fromDir = toDir; % use last todir
            if childCnt == 0 % leaf node: flip direction
                toDir = lastCirc - currCirc; %new todir
            else
                %none leaf: get dir to next node in stack 
                if~isempty(stack)
                    toDir = stack(1,:) - currCirc;
                else
                    % we are done! 
                    toDir = fromDir;
                end
            end
            circVisitCnt(currCirc(1),currCirc(2),currCirc(3)) =...
            circVisitCnt(currCirc(1),currCirc(2),currCirc(3)) + 1;
        end
      
        if isCurrCircDone
%          currCirc = currCirc + toDir; % go to next circuit
         lastCirc = currCirc;
        end                   
         
        while ~isCurrCircDone % finish circuit by circuit

            if fromDir == toDir
                % fill side
                currPix= visitPix3D(currPix,currCirc,[0,0,0], I);
%                 currPix = currPix + fromDir;
                currPix=visitPix3D(currPix,currCirc,fromDir, I);
                nextPix = currPix + toDir;
            else
                % fill side and turn side
                currPix=visitPix3D(currPix, currCirc,[0,0,0], I);
                if fromDir == -toDir
                    currPix=visitPix3D(currPix,currCirc,fromDir,I);
                    %if directions are opposite->this is the turn from
                    %forward to backward
                    %Add one more perpendicular edge!
                    if fromDir(1) ~= 0 % y dir
                        for t = 1:2
                            if t == 1
                                dir = [0,1, 0];
                            else
                                dir = [0,-1, 0];
                            end
                        
                            isInCirc = checkPixInCirc3D(currPix + dir, currCirc);
                            if isInCirc
%                                 currPix = currPix + dir;
                                currPix=visitPix3D(currPix,currCirc,dir, I);
                                break;
                            end
                        end
                    elseif fromDir(2) ~= 0 % x direction
                        for t = 1:2
                            if t == 1
                                dir = [1,0,0];
                            else
                                dir = [-1,0,0];
                            end
                        
                            isInCirc = checkPixInCirc3D(currPix + dir, currCirc);
                            if isInCirc
%                                 currPix = currPix + dir;
                               currPix=visitPix3D(currPix,currCirc,dir,I);
                                break;
                            end
                        end
                    else % z  direction
                         for t = 1:2
                            if t == 1
                                dir = [0,0,1];
                            else
                                dir = [0,0,-1];
                            end
                        
                            isInCirc = checkPixInCirc3D(currPix + dir, currCirc);
                            if isInCirc
%                                 currPix = currPix + dir;
                               currPix=visitPix3D(currPix,currCirc,dir,I);
                                break;
                            end
                        end
                    end
                    
                    currPix=visitPix3D(currPix,currCirc,toDir,I);
                    nextPix = currPix + toDir;
                else    
                    
                    % regular turn point
                   currPix= visitPix3D(currPix,currCirc,fromDir,I);       
                   currPix= visitPix3D(currPix,currCirc,toDir,I);
                    nextPix = currPix + toDir;
                end
            end

            % if next pixel is not in circuit
            if pixId >= length(LT)
                break;
            end
            inCirc = checkPixInCirc3D(nextPix, currCirc);
            if ~inCirc
                lastCirc = currCirc;
                isCurrCircDone = true;
                currPix = nextPix;
                
            end

        end

    end
    zLT = LT;
    zVisitOrder = visitOrder;
end
% helper functions

function isInside = checkPixInCirc3D(pix, circ)
   node = 2* circ - 1;
   pixInCirc = zeros(8,3);
   cnt = 1;
   for z = 1:2
       for y = 1:2
           for x = 1:2
               pixInCirc(cnt,:) = node + [y,x,z];
               cnt = cnt + 1;
           end
       end
   end

    isInside = false;
    for nn = 1:size(pixInCirc,1)
       if pix == pixInCirc(nn,:)
           isInside = true;
           break;
       end
    end   
end

function lastPix = visitPix3D(currPix, currCirc, dir, I)
global pixId
global LT;
global visitOrder;
global pixOcc;
modPix = currPix + dir;
if ~pixOcc(modPix(1),modPix(2),modPix(3)) 
    if checkPixInCirc3D(modPix, currCirc)
        LT(pixId,1) = I(modPix(1),modPix(2),modPix(3));
        % distance to center
        ctr = [size(pixOcc,1)/2, size(pixOcc,2)/2, size(pixOcc,3)/2];
        LT(pixId,2) = norm(ctr - modPix)/norm(ctr);
        visitOrder(pixId,:) = modPix;
        pixId = pixId + 1;
        pixOcc(modPix(1),modPix(2),modPix(3)) = true;
        lastPix = modPix;
    else
        lastPix = currPix;
    end
else
    lastPix = currPix;
end
end

% test if 2D projected line segments between 3D line segs aa->bb and cc->dd
function isProjIntersect = projIntersect(aa, bb, cc, dd)
   if max(aa(3),bb(3)) < min(cc(3),dd(3))
        isProjIntersect = false;
        return;
    end
    
   if max(aa(2),bb(2)) < min(cc(2),dd(2))
        isProjIntersect = false;
        return;
    end
    
    if max(aa(1),bb(1)) < min(cc(1),dd(1))
        isProjIntersect = false;
        return;
    end
    
   if max(cc(3),dd(3))<min(aa(3),bb(3))
        isProjIntersect = false;
        return;
   end
     
     if max(cc(2),dd(2))<min(aa(2),bb(2))
        isProjIntersect = false;
        return;
     end

    if max(cc(1),dd(1))<min(aa(1),bb(1))
        isProjIntersect = false;
        return;
    end
    
    val = cross(aa-bb, cc-dd);

    
    isProjIntersect = true;
end

function isValid = checkDirValid3D(dir, currPix, currCirc, lastCirc, children)
     isValid = true;
    
     global neighborDirs;
     global nneighbors;
     ctrPos = [2*(currCirc(1)-1)+1.5, 2*(currCirc(2)-1)+1.5, 2*(currCirc(3)-1)+1.5];
     children = reshape(children,[nneighbors,1]);
     if nnz(children) == 1
         isValid = true;
         return;
     else
         for i = 1:nneighbors
             if children(i) > 0
                  tdir = neighborDirs(i,:);
                  if dir(1)~=tdir(1) || dir(2)~=tdir(2) || dir(3)~=tdir(3)
                       if projIntersect(currPix, currPix+dir, ctrPos, ctrPos+tdir)
                         isValid = false;
                         return;
                       end
                  end
             end
         end
         %check link to the last circ
         tdir = lastCirc - currCirc;
         if dir(1)~=tdir(1) || dir(2)~=tdir(2) || dir(3)~=tdir(3)
               if projIntersect(currPix, currPix+dir, ctrPos, ctrPos+tdir)
                 isValid = false;
                 return;
               end
          end
     end
end
