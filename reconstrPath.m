
function  [visitOrder, LT] = reconstrPath(hamPath, V, sizeV)

w = sizeV(1);
h = sizeV(3);
d = sizeV(2);
    len = length(hamPath);
    visitOrder = zeros(len,3);
    LT = zeros(len,1);
    for i = 1:len
       pt = hamPath(i)-1;
       tz = floor(pt/(w*d));
       ty = floor( (pt - tz*w*d)/w);
       tx = (pt - (tz*d+ty)*w);
       x = tx+1;
       y = ty+1;
       z = tz+1;
       visitOrder(i,:) = [x,y,z];
       LT(i,:) = V(x,y,z);
    end
end