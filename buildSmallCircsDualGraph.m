
function [WpXR, WpXL, WpYU, WpYD] = buildSmallCircsDualGraph(I)
    %Wp is the matrix representation of weights for the dual graph of
    %small circuits
    [IdimY, IdimX, IdimZ] = size(I);
    WpXR = zeros(IdimY/2,IdimX/2-1);
    WpXL = zeros(IdimY/2,IdimX/2-1);
    WpYU = zeros(IdimY/2-1,IdimX/2);
    WpYD = zeros(IdimY/2-1,IdimX/2);
    %going right
    for y = 1:IdimY/2
        for x= 1:IdimX/2-1
%             %X direction: |u|+|w|+|x|+|y|+|z|-|e|-|f|
              xx = 2*(x-1)+1; % index of the original grid
              yy = 2*(y-1)+1;
            % going right
            WpXR(y,x) = abs(I(yy,xx+2)-I(yy,xx+1)) + abs(I(yy+1,xx+2)-I(yy+1,xx+1)) +...
                abs(I(yy,xx+3)-I(yy,xx+2)) + abs(I(yy+1,xx+3)-I(yy,xx+3))+abs(I(yy+1,xx+3)-I(yy+1,xx+2))-...
                abs(I(yy+1,xx+2)-I(yy,xx+2)) - abs(I(yy+1,xx+1)-I(yy,xx+1));
            % normalization?
            WpXR(y,x) = WpXR(y,x) ./ 3;
%                             
        end
    end
    
    %going left
     for y = 1:IdimY/2
        for x= IdimX/2-1:-1:1
%             %X direction: |u|+|w|+|x|+|y|+|z|-|e|-|f|
              xx = 2*(x)+1; % index of the original grid
              yy = 2*(y-1)+1;
            % going left
            WpXL(y,x) = abs(I(yy,xx-1)-I(yy,xx)) + abs(I(yy+1,xx-1)-I(yy+1,xx)) +...
                abs(I(yy,xx-2)-I(yy,xx-1)) + abs(I(yy+1,xx-2)-I(yy,xx-2))+abs(I(yy+1,xx-2)-I(yy+1,xx-1))-...
                abs(I(yy+1,xx-1)-I(yy,xx-1)) - abs(I(yy+1,xx)-I(yy,xx));
%                  
            % normalization?
            WpXL(y,x) = WpXL(y,x) ./ 3;
        end
    end
    
    % going down
    %Y direction: |e|+|f| + |x| +|y|+|z|-|u|-|w|
     for y = 1:IdimY/2-1
        for x= 1:IdimX/2
            xx = 2*(x-1)+1; % index of the original grid
            yy = 2*(y-1)+1;
            WpYD(y,x) = abs(I(yy+2,xx)-I(yy+1,xx)) + abs(I(yy+2,xx+1)-I(yy+1,xx+1))+... 
                abs(I(yy+3,xx)-I(yy+2,xx))+abs(I(yy+3,xx+1)-I(yy+2,xx+1))+abs(I(yy+3,xx+1)-I(yy+3,xx))-...
                     abs(I(yy+1,xx+1)-I(yy+1,xx)) - abs(I(yy+2,xx+1)-I(yy+2,xx));     
                             % normalization?
            WpYD(y,x) = WpYD(y,x) ./ 3;
        end
     end
    
    %going up
     for y = IdimY/2-1:-1:1
        for x= 1:IdimX/2
            xx = 2*(x-1)+1; % index of the original grid
            yy = 2*(y)+1;
            WpYU(y,x) = abs(I(yy-1,xx)-I(yy,xx)) + abs(I(yy-1,xx+1)-I(yy,xx+1))+... 
                abs(I(yy-2,xx)-I(yy-1,xx))+abs(I(yy-2,xx+1)-I(yy-1,xx+1))+abs(I(yy-2,xx+1)-I(yy-2,xx))-...
                     abs(I(yy,xx+1)-I(yy,xx)) - abs(I(yy-1,xx+1)-I(yy-1,xx));       
            % normalization?
            WpYU(y,x) = WpYU(y,x) ./ 3;
        end
     end
end
