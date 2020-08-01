function [I,rc] = line_plane_intersection(u, N, n, M, verbose)
% line_plane_intersection : function to compute the intersection point
% between the (N,u) line and the (M,n) plane of the 3D space.
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2020.
%
%% Syntax
%
% [I,rc] = line_plane_intersection(u, N, n, M);
% [I,rc] = line_plane_intersection(u, N, n, M, verbose);
%
%% Description
%
% [I,rc] = line_plane_intersection(u, N, n, M) computes the coordinates of I,
% the intersection point between the line (u,N) and the plane (n,M).
% In the most generic case, I is a point in the 3D space, but when
% the line is stricly parallel to the plane, I is the empty set, and when
% the line is included in the plane, I is a function handle corresponding
% to the system of parametric equations of the line.
%
% [I,rc] = line_plane_intersection(u, N, n, M, verbose) displays a message in
% console when verbose is set either to logical true or real numeric 1, and
% doesn't when it is set to logical false or real numeric 0.
%
%% Principle
%
% Based on solving Descartes plane equation :
%
% ax + by + cz + d = 0, where n = [a, b, c] is a vector normal to the plane,
%
% combined with the parametric equations system of a 3D line :
%
% x(t) = x0 + at 
% y(t) = y0 + bt
% z(t) = z0 + ct
%
% where N0 = [x0, y0, z0] is a point belonging to the line, and u = [a, b, c], a vector directing this line.
%
%% Input arguments
%
% - u : real numeric row or column vector. numel(u) = 3. One director vector of the parametric line.
%
% - N : real numeric row or column vector. numel(N) = 3. One point belonging to the line.
%
% - n : real numeric row or column vector. numel(n) = 3. One normal vector to the plane.
%
% - M : real numeric row or column vector. numel(M) = 3. One point belonging to the plane.
%
% - verbose : either logical *true/false or real numeric *1/0. 
%
%% Output arguments
%
% - I = [xI yI zI], real numeric row or column vector, the intersection point.
%
% - rc : return code, numeric integer in the set {1,2,3}.
%        0 = void / [] intersection
%        1 = point intersection
%        2 = line intersection
%
%        rc return code is necessary to distinguish between cases where
%        (N,u) line and the (M,n) plane intersection is a single point
%        and where it is the line itself.
%
%% Example #1 : one unique intersection point
%
% n = [1 1 1];
% M = n;
% u = [1 0 0];
% N = u; % (N,u) = (OX) axis
% [I,rc] = line_plane_intersection(u, N, n, M) % one unique intersction point expected : I = [3 0 0], rc = 1
%
%% Example #2 : line and plane are strictly // ; no intersection
%
% n = [0 0 1];
% M = [0 0 0]; % (M,n) = (XOY) plan
% u = [1 2 0];
% N = [0 0 6];
% [I,rc] = line_plane_intersection(u, N, n, M) % line strictly // plane =>  I = [], rc = 0 expected
%
%% Example #3 : line is included in the plane
%
% n = [1 1 1];
% M = (1/3)*[1 1 1];
% u = [1 1 -2];
% N = [0.5 0.5 0];
% [I,rc] = line_plane_intersection(u, N, n, M) % line belongs to the plane, rc = 2 expected
%% Input parsing
assert(nargin > 3,'Not enough input arguments.');
assert(nargin < 6,'Too many input arguments.');
if nargin < 5    
    verbose = true;    
else    
    assert(islogical(verbose) || isreal(verbose),'verbose must be of type either logical or real numeric.');    
end
assert(isequal(size(u),size(N),size(n),size(M)),'Inputs u, M, n, and M must have the same size.');
assert(isequal(numel(u),numel(N),numel(n),numel(M),3),'Inputs u, M, n, and M must have the same number of elements (3).');
assert(isequal(ndims(u),ndims(N),ndims(n),ndims(M)),'Inputs u, M, n, and M must have the same number of dimensions.');
%% Plane offset parameter
d = -dot(n,M);
%% Specific cases treatment
if ~dot(n,u) % n & u perpendicular vectors
    if dot(n,N) + d == 0 % N in P => line belongs to the plane
        if verbose
            disp('(N,u) line belongs to the (M,n) plane. Their intersection is the whole (N,u) line.');
        end
        I = M;
        rc = 2;
    else % line // to the plane
        if verbose
            disp('(N,u) line is parallel to the (M,n) plane. Their intersection is the empty set.');
        end
        I = [];
        rc = 0;
    end
else
    
    %% Parametric line parameter t
    t = - (d + dot(n,N)) / dot(n,u);
    
    %% Intersection coordinates
    I = N + u*t;
    
    rc = 1;
    
end
end % line_plane_intersection