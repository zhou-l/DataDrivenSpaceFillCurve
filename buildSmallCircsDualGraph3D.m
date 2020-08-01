
function [WpXR, WpXL, WpYU, WpYD, WpZF, WpZB] = buildSmallCircsDualGraph3D(I)
    %Wp is the matrix representation of weights for the dual graph of
    %small circuits
    [IdimY, IdimX, IdimZ] = size(I);
    WpXR = zeros(IdimY/2,IdimX/2-1,IdimZ/2);
    WpXL = zeros(IdimY/2,IdimX/2-1,IdimZ/2);
    WpYU = zeros(IdimY/2-1,IdimX/2,IdimZ/2);
    WpYD = zeros(IdimY/2-1,IdimX/2,IdimZ/2);
    WpZF = zeros(IdimY/2,IdimX/2,IdimZ/2-1);
    WpZB = zeros(IdimY/2,IdimX/2,IdimZ/2-1);
    %going right
    for z = 1:IdimZ/2
        for y = 1:IdimY/2
            for x= 1:IdimX/2-1
    %             %X direction: |u|+|w|+|x|+|y|+|z|-|e|-|f|
                xx = 2*(x-1)+1; % index of the original grid
                yy = 2*(y-1)+1;
                zz = 2*(z-1)+1;
                % going right 
                % first layer
                WpXR(y,x,z) = WpXR(y,x,z) + abs(I(yy,xx+2,zz)-I(yy,xx+1,zz)) + abs(I(yy+1,xx+2,zz)-I(yy+1,xx+1,zz)) +...
                    abs(I(yy,xx+3,zz)-I(yy,xx+2,zz)) + abs(I(yy+1,xx+3,zz)-I(yy,xx+3,zz))+abs(I(yy+1,xx+3,zz)-I(yy+1,xx+2,zz))-...
                    abs(I(yy+1,xx+2,zz)-I(yy,xx+2,zz)) - abs(I(yy+1,xx+1,zz)-I(yy,xx+1,zz));
                zz = zz + 1;
    %            second layer                 
                WpXR(y,x,z) = WpXR(y,x,z) + abs(I(yy,xx+2,zz)-I(yy,xx+1,zz)) + abs(I(yy+1,xx+2,zz)-I(yy+1,xx+1,zz)) +...
                abs(I(yy,xx+3,zz)-I(yy,xx+2,zz)) + abs(I(yy+1,xx+3,zz)-I(yy,xx+3,zz))+abs(I(yy+1,xx+3,zz)-I(yy+1,xx+2,zz))-...
                abs(I(yy+1,xx+2,zz)-I(yy,xx+2,zz)) - abs(I(yy+1,xx+1,zz)-I(yy,xx+1,zz));
                %third layer
                zz = zz - 1;
                WpXR(y,x,z) = WpXR(y,x,z) + abs(I(yy,xx+3,zz+1)-I(yy,xx+3,zz)) + abs(I(yy+1,xx+3,zz+1)-I(yy+1,xx+3,zz)) -...
                    abs(I(yy,xx+2,zz+1)-I(yy,xx+2,zz)) - abs(I(yy+1,xx+2,zz+1)-I(yy+1,xx+2,zz)) -...
                    abs(I(yy,xx+1,zz+1)-I(yy,xx+1,zz)) - abs(I(yy+1,xx+1,zz+1)-I(yy+1,xx+1,zz));
            end
        end
    end
    
    %going left
    for z = 1:IdimZ/2
         for y = 1:IdimY/2
            for x= IdimX/2-1:-1:1
    %             %X direction: |u|+|w|+|x|+|y|+|z|-|e|-|f|
                  xx = 2*(x)+1; % index of the original grid
                  yy = 2*(y-1)+1;
                  zz = 2*(z-1)+1;
                % going left
                % first layer
                WpXL(y,x,z) = abs(I(yy,xx-1,zz)-I(yy,xx,zz)) + abs(I(yy+1,xx-1,zz)-I(yy+1,xx,zz)) +...
                    abs(I(yy,xx-2,zz)-I(yy,xx-1,zz)) + abs(I(yy+1,xx-2,zz)-I(yy,xx-2,zz))+abs(I(yy+1,xx-2,zz)-I(yy+1,xx-1,zz))-...
                    abs(I(yy+1,xx-1,zz)-I(yy,xx-1,zz)) - abs(I(yy+1,xx,zz)-I(yy,xx,zz));
    %             second layer      
                zz = zz + 1;
                WpXL(y,x,z) = WpXL(y,x,z) + abs(I(yy,xx-1,zz)-I(yy,xx,zz)) + abs(I(yy+1,xx-1,zz)-I(yy+1,xx,zz)) +...
                    abs(I(yy,xx-2,zz)-I(yy,xx-1,zz)) + abs(I(yy+1,xx-2,zz)-I(yy,xx-2,zz))+abs(I(yy+1,xx-2,zz)-I(yy+1,xx-1,zz))-...
                    abs(I(yy+1,xx-1,zz)-I(yy,xx-1,zz)) - abs(I(yy+1,xx,zz)-I(yy,xx,zz));       
                %third layer
                zz = zz - 1; 
                WpXL(y,x,z) = WpXL(y,x,z) + abs(I(yy,xx-2,zz+1)-I(yy,xx-2,zz)) + abs(I(yy+1,xx-2,zz+1)-I(yy+1,xx-2,zz)) -...
                    abs(I(yy,xx-1,zz+1)-I(yy,xx-1,zz)) - abs(I(yy+1,xx-1,zz+1)-I(yy+1,xx-1,zz)) -...
                    abs(I(yy,xx,zz+1)-I(yy,xx,zz)) - abs(I(yy+1,xx,zz+1)-I(yy+1,xx,zz));
            end
         end
    end
    
    % going down
    %Y direction: |e|+|f| + |x| +|y|+|z|-|u|-|w|
    for z = 1:IdimZ/2
     for y = 1:IdimY/2-1
        for x= 1:IdimX/2
            xx = 2*(x-1)+1; % index of the original grid
            yy = 2*(y-1)+1;
            zz = 2*(z-1)+1;
            % first layer
            WpYD(y,x,z) = abs(I(yy+2,xx,zz)-I(yy+1,xx,zz)) + abs(I(yy+2,xx+1,zz)-I(yy+1,xx+1,zz))+... 
                abs(I(yy+3,xx,zz)-I(yy+2,xx,zz))+abs(I(yy+3,xx+1,zz)-I(yy+2,xx+1,zz))+abs(I(yy+3,xx+1,zz)-I(yy+3,xx,zz))-...
                     abs(I(yy+1,xx+1,zz)-I(yy+1,xx,zz)) - abs(I(yy+2,xx+1,zz)-I(yy+2,xx,zz));      
            % second layer
            zz = zz + 1;
            WpYD(y,x,z) = WpYD(y,x,z) + abs(I(yy+2,xx,zz)-I(yy+1,xx,zz)) + abs(I(yy+2,xx+1,zz)-I(yy+1,xx+1,zz))+... 
                abs(I(yy+3,xx,zz)-I(yy+2,xx,zz))+abs(I(yy+3,xx+1,zz)-I(yy+2,xx+1,zz))+abs(I(yy+3,xx+1,zz)-I(yy+3,xx,zz))-...
                     abs(I(yy+1,xx+1,zz)-I(yy+1,xx,zz)) - abs(I(yy+2,xx+1,zz)-I(yy+2,xx,zz));   
            %third layer
            zz = zz - 1;
            WpYD(y,x,z) = WpYD(y,x,z) + abs(I(yy+3,xx,zz+1)-I(yy+3,xx,zz)) + abs(I(yy+3,xx+1,zz+1)-I(yy+3,xx+1,zz)) -...
                    abs(I(yy+2,xx,zz+1)-I(yy+2,xx,zz)) - abs(I(yy+2,xx+1,zz+1)-I(yy+2,xx+1,zz)) -...
                    abs(I(yy+1,xx,zz+1)-I(yy+1,xx,zz)) - abs(I(yy+1,xx+1,zz+1)-I(yy+1,xx+1,zz));
        end
     end
    end
    
    %going up
    for z = 1:IdimZ/2
         for y = IdimY/2-1:-1:1
            for x= 1:IdimX/2
                xx = 2*(x-1)+1; % index of the original grid
                yy = 2*(y)+1;
                zz = 2*(z-1)+1;
                % first layer
                WpYU(y,x,z) = abs(I(yy-1,xx,zz)-I(yy,xx,zz)) + abs(I(yy-1,xx+1,zz)-I(yy,xx+1,zz))+... 
                    abs(I(yy-2,xx,zz)-I(yy-1,xx,zz))+abs(I(yy-2,xx+1,zz)-I(yy-1,xx+1,zz))+abs(I(yy-2,xx+1,zz)-I(yy-2,xx,zz))-...
                    abs(I(yy,xx+1,zz)-I(yy,xx,zz)) - abs(I(yy-1,xx+1,zz)-I(yy-1,xx,zz));   
                % second layer
                zz = zz + 1;
                WpYU(y,x,z) = WpYU(y,x,z) + abs(I(yy-1,xx,zz)-I(yy,xx,zz)) + abs(I(yy-1,xx+1,zz)-I(yy,xx+1,zz))+... 
                    abs(I(yy-2,xx,zz)-I(yy-1,xx,zz))+abs(I(yy-2,xx+1,zz)-I(yy-1,xx+1,zz))+abs(I(yy-2,xx+1,zz)-I(yy-2,xx,zz))-...
                    abs(I(yy,xx+1,zz)-I(yy,xx,zz)) - abs(I(yy-1,xx+1,zz)-I(yy-1,xx,zz));
                % third layer
                zz = zz - 1;
                WpYU(y,x,z) = WpYU(y,x,z) + abs(I(yy-2,xx,zz+1)-I(yy-2,xx,zz)) + abs(I(yy-2,xx+1,zz+1)-I(yy-2,xx+1,zz)) -...
                    abs(I(yy-1,xx,zz+1)-I(yy-1,xx,zz)) - abs(I(yy-1,xx+1,zz+1)-I(yy-1,xx+1,zz)) -...
                    abs(I(yy,xx,zz+1)-I(yy,xx,zz)) - abs(I(yy,xx+1,zz+1)-I(yy,xx+1,zz));
            end
         end
    end
    
    %going forward
    for z = IdimZ/2-1:-1:1
         for y = 1:IdimY/2
            for x= 1:IdimX/2
                xx = 2*(x-1)+1; % index of the original grid
                yy = 2*(y-1)+1;
                zz = 2*z + 1;
                % first layer
                WpZF(y,x,z) = abs(I(yy,xx,zz-1)-I(yy,xx,zz)) + abs(I(yy,xx+1,zz-1)-I(yy,xx+1,zz))+... 
                    abs(I(yy,xx,zz-2)-I(yy,xx,zz-1))+abs(I(yy,xx+1,zz-2)-I(yy,xx+1,zz-1))+abs(I(yy,xx+1,zz-2)-I(yy,xx,zz-2))-...
                         abs(I(yy,xx+1,zz-1)-I(yy,xx,zz-1)) - abs(I(yy,xx+1,zz)-I(yy,xx+1,zz));       
                %second layer
                yy = yy + 1;
                WpZF(y,x,z) = WpZF(y,x,z) + abs(I(yy,xx,zz-1)-I(yy,xx,zz)) + abs(I(yy,xx+1,zz-1)-I(yy,xx+1,zz))+... 
                    abs(I(yy,xx,zz-2)-I(yy,xx,zz-1))+abs(I(yy,xx+1,zz-2)-I(yy,xx+1,zz-1))+abs(I(yy,xx+1,zz-2)-I(yy,xx,zz-2))-...
                         abs(I(yy,xx+1,zz-1)-I(yy,xx,zz-1)) - abs(I(yy,xx+1,zz)-I(yy,xx+1,zz));    
                %third layer
                yy = yy - 1;
                WpZF(y,x,z) = WpZF(y,x,z) + abs(I(yy+1,xx,zz-2)-I(yy,xx,zz-2)) + abs(I(yy+1,xx+1,zz-2)-I(yy,xx+1,zz-2)) -...
                    abs(I(yy+1,xx,zz-1)-I(yy,xx,zz-1)) - abs(I(yy+1,xx+1,zz-1)-I(yy,xx+1,zz-1)) -...
                    abs(I(yy+1,xx,zz)-I(yy,xx,zz)) - abs(I(yy+1,xx+1,zz)-I(yy,xx+1,zz));
            end
         end
    end
        
    %going back
    for z = 1:IdimZ/2-1
         for y = 1:IdimY/2
            for x= 1:IdimX/2
                xx = 2*(x-1)+1; % index of the original grid
                yy = 2*(y-1)+1;
                zz = 2*(z-1)+1;
                % first layer
                WpZB(y,x,z) = abs(I(yy,xx,zz+2)-I(yy,xx,zz+1)) + abs(I(yy,xx+1,zz+2)-I(yy,xx+1,zz+1))+... 
                    abs(I(yy,xx,zz+3)-I(yy,xx,zz+2))+abs(I(yy,xx+1,zz+3)-I(yy,xx+1,zz+2))+abs(I(yy,xx+1,zz+3)-I(yy,xx,zz+3))-...
                    abs(I(yy,xx+1,zz+2)-I(yy,xx,zz+2)) - abs(I(yy,xx+1,zz+1)-I(yy,xx+1,zz+1));
                %second layer
                yy = yy + 1;
                WpZB(y,x,z) = WpZB(y,x,z) + abs(I(yy,xx,zz+2)-I(yy,xx,zz+1)) + abs(I(yy,xx+1,zz+2)-I(yy,xx+1,zz+1))+... 
                    abs(I(yy,xx,zz+3)-I(yy,xx,zz+2))+abs(I(yy,xx+1,zz+3)-I(yy,xx+1,zz+2))+abs(I(yy,xx+1,zz+3)-I(yy,xx,zz+3))-...
                    abs(I(yy,xx+1,zz+2)-I(yy,xx,zz+2)) - abs(I(yy,xx+1,zz+1)-I(yy,xx+1,zz+1));
                %third layer
                yy = yy - 1;
                WpZB(y,x,z) = WpZB(y,x,z) + abs(I(yy+1,xx,zz+2)-I(yy,xx,zz+2)) + abs(I(yy+1,xx+1,zz+2)-I(yy,xx+1,zz+2)) -...
                    abs(I(yy+1,xx,zz+1)-I(yy,xx,zz+1)) - abs(I(yy+1,xx+1,zz+1)-I(yy,xx+1,zz+1)) -...
                    abs(I(yy+1,xx,zz)-I(yy,xx,zz)) - abs(I(yy+1,xx+1,zz)-I(yy,xx+1,zz));
            end
         end
    end
    
end
