function [hexAmpConst, hexPhzConst] = generateHexConstants(cenrow, cencol, numRings, apDia, wGap, arrayIn)

 hexFlatDiam = (apDia-numRings*2*wGap)/(2*numRings+1);
    % hexRad = hexFlatDiam/sqrt(3);% center to vertex
    hexSep = hexFlatDiam + wGap;

    [rows,cols]=size(arrayIn);
    

    [X,Y] = meshgrid(-cols/2:cols/2-1,-rows/2:rows/2-1); % Grids with Cartesian (x,y) coordinates
%     disp(X);
%     disp(Y);

    RHOprime = sqrt((X-cencol).^2+(Y-cenrow).^2);
    %disp(size(RHOprime));
    THETA = atan2(Y-cenrow,X-cencol);
    %disp(size(sin(THETA)));

    hexAmpConst = exp(-(RHOprime.*sin(THETA)/(hexFlatDiam/2)).^1000)...
            .*exp(-(RHOprime.*cos(THETA-pi/6)/(hexFlatDiam/2)).^1000)...
            .*exp(-(RHOprime.*cos(THETA+pi/6)/(hexFlatDiam/2)).^1000);
    
    hexPhzConst = and(and(RHOprime.*sin(THETA)<=(hexSep/2),...
             RHOprime.*cos(THETA-pi/6)<=(hexSep/2)),...
             RHOprime.*cos(THETA+pi/6)<=(hexSep/2));
         
end