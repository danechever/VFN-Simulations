function [optField,circum_diam] = makeKeckPupilInputs( vars, seed )
    %Generates a Keck pupil based on a set of inputs in 'vars'
    
    addpath([fileparts(which(mfilename)),'/segMirrorFunctions']);
   
    hexMirror.apDia = 2*vars.apRad; % flat to flat aperture diameter (samples)
    hexMirror.wGap = vars.wGap; % gap width in samples
    hexMirror.numRings = vars.numRings;% Number of rings in hexagonally segmented mirror
    
    hexMirror.Npad = vars.PUP_CRP_SZ; %vars.N;
    central_band_index = ceil(vars.numWavelengths/2);
    hexMirror.charge = vars.charge(central_band_index);
    
    hexMirror.tiltxs = seed(:,1);
    hexMirror.tiltys = seed(:,2);
    hexMirror.pistons = seed(:,3);
    %hexMirror.hexSep = vars.hexSep;
    
    hexMirror.hexAmpConst = vars.hexAmpConst;
    hexMirror.hexPhzConst = vars.hexPhzConst;
    
    [X,Y] = meshgrid(-vars.N/2:vars.N/2-1);
    [THETA,RHO] = cart2pol(X,Y);
    
    %The optimization function will be optimizing angle(optPhz)
    %Thus, the inputs to that function are what will be changed
    %Tiltx,Tilty, and Piston are the 3 variables it should be optimized
    %around

    optField = hexSegMirror_getField_Modified( hexMirror ); %Field in the pupil plane...called 'optField' because it uses optimized inputs
    
    % Assume matrix was shrunk and pad to size if needed
    optField = pad_crop(optField, vars.N); 

    [rowslist,colslist] = find(round(optField));
    apRad = max(sqrt((rowslist-vars.N/2-1).^2 + (colslist-vars.N/2-1).^2));
    circum_diam = hexMirror.apDia;

    secRad = apRad*2600/10949;
    spwidth = 2*apRad*25.4/10949;
    optField = optField.*(1-exp(-(RHO/secRad).^1000));
    optField = optField.*(1-exp(-(Y/(spwidth/2)).^1000));
    optField = optField.*(1-exp(-(RHO.*cos(THETA-30*pi/180)/(spwidth/2)).^1000));
    optField = optField.*(1-exp(-(RHO.*cos(THETA+30*pi/180)/(spwidth/2)).^1000));
end