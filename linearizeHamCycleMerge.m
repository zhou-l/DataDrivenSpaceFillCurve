%Traversing the image I using a Hamilton cycle generated from a precomputed minimum spanning tree
% use tree as input
function [zLT,zVisitOrder] = linearizeHamCycleMerge(mstT, I, entryPix)
    currPix = entryPix;
    currCirc = [floor((currPix(1)+1)/2),floor((currPix(2)+1)/2)];

%     dirSeq = genDirSeq(mstSet);
    toDir = [-1,0];
    stack = [currCirc]; 
    dirStack = [toDir];
    
    hamCyc = zeros(size(I,1)*size(I,2),2);
%     hamCyc(1,:) = [2* (currCirc(1) - 1)+1, 2*(currCirc(2) - 1)+1];
%     hamCyc(2,:) = hamCyc(1,:) + [0,1];
%     hamCyc(3,:) = hamCyc(2,:) + [-1,0];
%     hamCyc(4,:) = hamCyc(3,:) + [0,-1];
%     
    cnt = 0;
    while ~isempty(stack)
        cnt = cnt + 1;
        currCirc = stack(1,1:2);   
        toDir = dirStack(1,1:2);
%         numvisited = circVisitCnt(currCirc(1),currCirc(2));
        % remove first
        stack(1,:) = [];
        dirStack(1,:) = [];
        children = mstT(currCirc(1),currCirc(2),:);
        children = reshape(children, [4,1]);
%         childCnt = nnz(children);
        
        hamCyc = mergeToCyc(hamCyc, currCirc, toDir);
        
        for i = 1:4
            if children(i) > 0
                   switch(i) % for the current channel==direction version
                   case 1%up
                      tryDir = [-1,0];
                   case 2%down
                      tryDir = [1,0];
                   case 3%right
                       tryDir = [0,1];
                   case 4%left
                       tryDir = [0,-1];  
                   end
                   nextCirc = currCirc + tryDir;        
                   toDir = tryDir;
                      % if next pixel is not in circuit, we visit
                      % current pixel and skip all
                  stack = [nextCirc;stack];
                  dirStack = [toDir;dirStack];
            end
        end
       
    end
    cnt
    zVisitOrder = hamCyc;
    
    % traverse using the hamCycle
    zLT = zeros(length(hamCyc),1);
    for i = 1:length(hamCyc)
        zLT(i,:) = I(hamCyc(i,1),hamCyc(i,2),:);
    end
%     zVisitOrder = visitOrder;
end
% helper functions
function hamCycNew = mergeToCyc(hamCycOld, currCirc, mergeDir)
   hamCycNew = hamCycOld;
   % number of existing pixels in the hampath
   lastN = find(~any(hamCycOld,2),1)-1;
 
   edgeBreakId = 0;

   % get pixels in the circuit to be merged
   pCC = zeros(4,2); 
   pCC(1,:) = [2*(currCirc(1)-1)+1, 2*(currCirc(2)-1)+1]; % top left
   pCC(2,:) = pCC(1,:) + [1,0]; % bottom left
   pCC(3,:) = pCC(2,:) + [0,1]; % bottom right
   pCC(4,:) = pCC(3,:) + [-1,0]; %top right
   
  % no pixels in the cycle yet!  
  if lastN == 0 
     
      tDir = mergeDir;
      % first node in hamCyc
        if tDir(1) ~= 0
            % merge in y direction
            if tDir(1) == 1 % going down
                    hamCycNew(edgeBreakId+1,:) = pCC(1,:); % top left
                    hamCycNew(edgeBreakId+2,:) = pCC(2,:); % bottom left
                    hamCycNew(edgeBreakId+3,:) = pCC(3,:); % bottom right
                    hamCycNew(edgeBreakId+4,:) = pCC(4,:); % top right

            elseif tDir(1) == -1 % going up
                    hamCycNew(edgeBreakId+1,:) = pCC(2,:); % bottom left 
                    hamCycNew(edgeBreakId+2,:) = pCC(1,:); % top left
                    hamCycNew(edgeBreakId+4,:) = pCC(3,:); % bottom right
                    hamCycNew(edgeBreakId+3,:) = pCC(4,:); % top right
            end
        else 
            % merge in x direction
            if tDir(2) == 1 % going right
                    hamCycNew(edgeBreakId+1,:) = pCC(2,:);  % bottom left
                    hamCycNew(edgeBreakId+4,:) = pCC(1,:); % top left
                    hamCycNew(edgeBreakId+2,:) = pCC(3,:); % bottom right
                    hamCycNew(edgeBreakId+3,:) = pCC(4,:); % top right
            elseif tDir(2) == -1 % going left
                     hamCycNew(edgeBreakId+1,:) = pCC(4,:); % top right
                    hamCycNew(edgeBreakId+2,:) = pCC(1,:); % top left
                    hamCycNew(edgeBreakId+3,:) = pCC(2,:); % bottom left
                    hamCycNew(edgeBreakId+4,:) = pCC(3,:); % bottom right
            end
        end
        return;
  end
  
     % find circ in hamCyc that adjacents to currCirc
   edgeFound = false;
   firstNodeInCirc = 0;
   for i = 1:lastN
       if edgeFound
           break;
       end
       for k = 1:4
           if norm(hamCycOld(i,:) - pCC(k,:)) == 1 
               tDir = (pCC(k,:) - hamCycOld(i,:));
                if tDir(1) == mergeDir(1) && tDir(2) == mergeDir(2)
                   edgeBreakId = i;
                   firstNodeInCirc = k;
                   edgeFound = true;
                   break;
                end
           end
       end
   end
   
   % copy existing cycle
   if edgeBreakId ~= 0
    hamCycNew(1:edgeBreakId,:) = hamCycOld(1:edgeBreakId,:);
    hamCycNew(edgeBreakId+5:lastN+4,:) = hamCycOld(edgeBreakId+1:lastN,:);
   end
   
      % merge at edgeBreakId

  
    hamCycNew(edgeBreakId+1,:) = pCC(firstNodeInCirc,:); % first node
     if tDir(1) ~= 0
        % merge in y direction
        if tDir(1) == 1 % going down
            if firstNodeInCirc == 1 % first node is top left
                hamCycNew(edgeBreakId+2,:) = pCC(2,:); % bottom left
                hamCycNew(edgeBreakId+3,:) = pCC(3,:); % bottom right
                hamCycNew(edgeBreakId+4,:) = pCC(4,:); % top right
            else % or top right
                hamCycNew(edgeBreakId+3,:) = pCC(2,:); % bottom left
                hamCycNew(edgeBreakId+2,:) = pCC(3,:); % bottom right
                hamCycNew(edgeBreakId+4,:) = pCC(1,:); % top left
            end
        elseif tDir(1) == -1 % going up
            if firstNodeInCirc == 2 % bottom left
                hamCycNew(edgeBreakId+2,:) = pCC(1,:); % top left
                hamCycNew(edgeBreakId+4,:) = pCC(3,:); % bottom right
                hamCycNew(edgeBreakId+3,:) = pCC(4,:); % top right
            else % bottom right
                hamCycNew(edgeBreakId+3,:) = pCC(1,:); % top left
                hamCycNew(edgeBreakId+4,:) = pCC(2,:); % bottom left
                hamCycNew(edgeBreakId+2,:) = pCC(4,:); % top right
            end
        end
   else 
        % merge in x direction
        if tDir(2) == 1 % going right
            if firstNodeInCirc == 2 % bottom left
                hamCycNew(edgeBreakId+4,:) = pCC(1,:); % top left
                hamCycNew(edgeBreakId+2,:) = pCC(3,:); % bottom right
                hamCycNew(edgeBreakId+3,:) = pCC(4,:); % top right
            else % top left
                hamCycNew(edgeBreakId+4,:) = pCC(2,:); % bottom left
                hamCycNew(edgeBreakId+3,:) = pCC(3,:); % bottom right
                hamCycNew(edgeBreakId+2,:) = pCC(4,:); % top right
            end
        elseif tDir(2) == -1 % going left
            if firstNodeInCirc == 4 % top right
                hamCycNew(edgeBreakId+2,:) = pCC(1,:); % top left
                hamCycNew(edgeBreakId+3,:) = pCC(2,:); % bottom left
                hamCycNew(edgeBreakId+4,:) = pCC(3,:); % bottom right
            else % bottom right
                hamCycNew(edgeBreakId+3,:) = pCC(1,:); % top left
                hamCycNew(edgeBreakId+2,:) = pCC(2,:); % bottom left
                hamCycNew(edgeBreakId+4,:) = pCC(4,:); % top right
            end
        end
    end
end
