% test intersection between line seg aa->bb and cc->dd
function isIntersect = intersectLL2D(aa, bb, cc, dd)
    
    if max(aa(2),bb(2)) < min(cc(2),dd(2))
        isIntersect = false;
        return;
    end
    
    if max(aa(1),bb(1)) < min(cc(1),dd(1))
        isIntersect = false;
        return;
    end
    
     if max(cc(2),dd(2))<min(aa(2),bb(2))
        isIntersect = false;
        return;
     end

    if max(cc(1),dd(1))<min(aa(1),bb(1))
        isIntersect = false;
        return;
    end
    
    if mult(cc,bb,aa)*mult(bb,dd,aa)<0
        isIntersect = false;
        return;
    end
    
    if mult(aa,bb,cc)*mult(dd,bb,cc)<0
        isIntersect = false;
        return;
    end
    
    isIntersect = true;   
end

function res = mult(a, b, c)
    res = (a(2)-c(2))*(b(1)-c(1))-(b(2)-c(2))*(a(1)-c(1));
    return;
end