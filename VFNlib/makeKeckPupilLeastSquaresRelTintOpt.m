function relTints = makeKeckPupilLeastSquaresRelTintOpt( vars, inputs, fiber_props, fibmodes)
    addpath([fileparts(which(mfilename)),'/segMirrorFunctions']);
    addpath(['..' filesep '..' filesep 'VFN-Lab' filesep 'AnalysisCode' filesep 'AnalysisLib']); 
    
%     persistent n_calls;
%     if isempty(n_calls)
%         n_calls=0;
%     end
%     global relTintArray;
     
    %-- Get Ideal Phase Mask
    fiber = fiber_props;
    
    %-- Get inputs
    hexMirror.apDia = 2*vars.apRad; % flat to flat aperture diameter (samples)
    hexMirror.wGap = vars.wGap; % gap width in samples
    hexMirror.numRings = vars.numRings;% Number of rings in hexagonally segmented mirror
    
    hexMirror.Npad = vars.PUP_CRP_SZ;%vars.N;
    central_band_index = ceil(vars.numWavelengths/2);
    hexMirror.charge = vars.charge(central_band_index);
    
    hexMirror.tiltxs = inputs(:,1);
    hexMirror.tiltys = inputs(:,2);
    hexMirror.pistons = inputs(:,3);
    
    hexMirror.hexAmpConst = vars.hexAmpConst;
    hexMirror.hexPhzConst = vars.hexPhzConst;
    
    [X,Y] = meshgrid(-vars.N/2:vars.N/2-1);
    [THETA,RHO] = cart2pol(X,Y);
    
    %-- Generates the segments that will be added to the primary
    Efield_opt = hexSegMirror_getField_Modified( hexMirror );
    
    Efield_opt = pad_crop(Efield_opt, vars.N);

    [rowslist,colslist] = find(round(Efield_opt));
    apRad = max(sqrt((rowslist-vars.N/2-1).^2 + (colslist-vars.N/2-1).^2));
    %circum_diam = hexMirror.apDia;

    %-- Adds apodization around edges of segments
    secRad = apRad*2600/10949;
    spwidth = 2*apRad*25.4/10949;
    Efield_opt = Efield_opt.*(1-exp(-(RHO/secRad).^1000));
    Efield_opt = Efield_opt.*(1-exp(-(Y/(spwidth/2)).^1000));
    Efield_opt = Efield_opt.*(1-exp(-(RHO.*cos(THETA-30*pi/180)/(spwidth/2)).^1000));
    Efield_opt = Efield_opt.*(1-exp(-(RHO.*cos(THETA+30*pi/180)/(spwidth/2)).^1000));
    
    optPhz = angle(Efield_opt);
    
    %--Add wedge effects to the phase here
    %optPhz = optPhz + getWedgeOffset(theta, n);
    
    for ch = 1:vars.numWavelengths
        Epup_opt(:,:,ch) = exp(1i*optPhz*vars.lambda0/vars.lambdas(ch)).*vars.PUPIL;
    end
    
    %-- Get eta_s
    eta_maps = generateCouplingMap_polychromatic( Epup_opt, fiber_props, vars.lambda0, fiber.Fnum, vars.lambdas, vars.totalPower0, vars.lambdaOverD, 3*vars.lambdaOverD, vars.coords, fibmodes);
    
    %-- Average the coupling maps
    map = mean(eta_maps,3);
    
    crp = 2*0.5*vars.lambdaOverD; %The length of one side of the cube to crop the image to

    cmap = map(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
    eta_s = min(cmap,[],'all'); %minimum value in the centroid
    [min_ind(1),min_ind(2)] = find(eta_s==cmap); %indices of minimum value in the centroid
    min_ind = min_ind + (vars.N/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
        
    eta_sX = min_ind(2);
    eta_sY = min_ind(1);

    [radProf2, ~] = VFN_An_radAverage(map,[eta_sX, eta_sY]);
    [radMx, ~] = max(radProf2);
    
    %-- Compute tip/tilt sensitivity for optimizer
    null_border = 1;
    nullgrid = map(eta_sY-null_border:eta_sY+null_border, eta_sX-null_border:eta_sX+null_border);
    %disp(nullReg);
    t_tSens = median(nullgrid, 'all');
    %disp(t_tSens);
    
    %-- Calculate relative integration time in each frame
    relTints = eta_s./(radMx.^2) + t_tSens;
    %disp(relTints);
%     tf = isscalar(relTints);
%     disp(tf);
%     vars.n_calls=vars.n_calls+1;
%     
%     vars.relTintArray(vars.n_calls) = relTints;
end