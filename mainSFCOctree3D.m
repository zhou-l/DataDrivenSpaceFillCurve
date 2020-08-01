% main function of Data-Driven Space-Filling Curves on octree 
% input: 
% filename---the input volume or point-cloud file; 
%            a test case will be used if no filename provided.
function mainSFCOctree3D(filename)
    if nargin < 1
        [clLT, clVisitOrder, fullLT] = SFCOctTreeMultiScaleMain();
    else
        [clLT, clVisitOrder, fullLT] = SFCOctTreeMultiScaleMain(filename);
    end
end