function figGenSFCOctree3D(filename)
    if nargin < 1
        [clLT, clVisitOrder, fullLT] = SFCOctTreeMultiScaleMain();
    else
        [clLT, clVisitOrder, fullLT] = SFCOctTreeMultiScaleMain(filename);
    end
end