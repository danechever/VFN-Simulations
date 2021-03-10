% Copyright 2018, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [ OUT ] = hexSegMirror_getField_Modified( hexMirror_struct )
%hexSegMirror_getSupport Returns the support of the pupil function defined
%by a hexagonally segmented mirror 
%   Input: hexMirror_struct - Structure with the following variables 
%   apDia - flat to flat aperture diameter (samples)
%   wGap - width of the gap between segments (samples)
%   numRings - number of rings in the segmented mirror (samples)
%   N - size of NxN computational grid 
%   missingSegments - list of zeros and ones indicating if each segment is present 
%   offset - centering offset vector [N/2+1+offset(1), N/2+1+offset(2)]

apDia = hexMirror_struct.apDia; % flat to flat aperture diameter (samples)
wGap = hexMirror_struct.wGap; % samples
numRings = hexMirror_struct.numRings;% Number of rings in hexagonally segmented mirror 
N = hexMirror_struct.Npad;
charge = hexMirror_struct.charge;

%Tilt (X and Y) and Piston
tiltx = hexMirror_struct.tiltxs;
tilty = hexMirror_struct.tiltys;
piston = hexMirror_struct.pistons;

%Hex Constants
hexAmpConst = hexMirror_struct.hexAmpConst;
hexPhzConst = hexMirror_struct.hexPhzConst;

%Determines number of segments to add to the pupil
if(isfield(hexMirror_struct,'missingSegments'))
    missingSegments = hexMirror_struct.missingSegments;
else
    missingSegments = ones(1,hexSegMirror_numSegments( numRings ));%returns number of segments in the pupil
end

%Grid information for pupil
OUT = zeros(N); %This is the matrix output by this function
[rows,cols] = size(OUT);
[grdX,grdY] = meshgrid(-cols/2:cols/2-1,-rows/2:rows/2-1); % Grids with Cartesian (x,y) coordinates

%Hex information for generating them
hexFlatDiam = (apDia-numRings*2*wGap)/(2*numRings+1);
hexSep = hexFlatDiam + wGap;

%This loop adds individual segments to the pupil
count = 1;
for ringNum = 1:numRings
    cenrow = ringNum*hexSep;
    cencol = 0;
    
    if(missingSegments(count)==1)
        %This part adds the 1st segment in each ring.
        [OUT] = hexSegMirror_addHexSegment_Modified(cenrow,cencol,apDia,...
            charge*piston(count),charge*tiltx(count),charge*tilty(count),OUT,...
            grdX, grdY, hexAmpConst(:,:,count), hexPhzConst(:,:,count) );
    end
    count = count + 1;
    
    for face = 1:6
        step_dir = pi/6*(2*face+5);
        steprow = hexSep*sin(step_dir);
        stepcol = hexSep*cos(step_dir);
        
        stepnum = 1;

        while(stepnum <= ringNum)
            cenrow = cenrow + steprow;
            cencol = cencol + stepcol;
            
            if(face==6 && stepnum==ringNum)
                %disp(['Finished ring ',num2str(ringNum)]);
            else
                if(missingSegments(count)==1)
                    %This part adds the majority of the segments
                    [ OUT ] = hexSegMirror_addHexSegment_Modified(cenrow,cencol,apDia,...
                        charge*piston(count),charge*tiltx(count),charge*tilty(count),OUT,...
                        grdX, grdY, hexAmpConst(:,:,count), hexPhzConst(:,:,count) );
                end
                count = count + 1;
            end
            stepnum = stepnum + 1;
            
        end
    end
end