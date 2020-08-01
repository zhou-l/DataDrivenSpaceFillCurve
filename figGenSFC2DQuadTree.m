function figGenSFC2DQuadTree(filename)
    if nargin < 1
        SFCQuadTreeMultiScaleMain();
    else
        SFCQuadTreeMultiScaleMain(filename);
    end
end