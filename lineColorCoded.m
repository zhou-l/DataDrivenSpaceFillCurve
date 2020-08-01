function lineColorCoded(x, y, z, col)
    % x, y, z, col are row vectors
     surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    