function [relTints,circum_diam] = makeKeckPupilLeastSquaresRelTintOpt( vars, inputs, fiber_props, fibmodes)

    addpath([fileparts(which(mfilename)),'/segMirrorFunctions']);
    addpath(['..' filesep '..' filesep 'VFN-lab' filesep 'AnalysisCode' filesep 'AnalysisLib']); 
    
%     persistent n_calls;
%     if isempty(n_calls)
%         n_calls=0;
%     end
%     global relTintArray;
%     
    % Get Ideal Phase Mask
    fiber = fiber_props;
    
    % Get small inputs
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
    
    %Generates the segments that will be added to the primary
    Efield_opt = hexSegMirror_getField_Modified( hexMirror );
    
    Efield_opt = pad_crop(Efield_opt, vars.N);

    [rowslist,colslist] = find(round(Efield_opt));
    apRad = max(sqrt((rowslist-vars.N/2-1).^2 + (colslist-vars.N/2-1).^2));
    circum_diam = hexMirror.apDia;

    %Adds apodization around edges of segments
    secRad = apRad*2600/10949;
    spwidth = 2*apRad*25.4/10949;
    Efield_opt = Efield_opt.*(1-exp(-(RHO/secRad).^1000));
    Efield_opt = Efield_opt.*(1-exp(-(Y/(spwidth/2)).^1000));
    Efield_opt = Efield_opt.*(1-exp(-(RHO.*cos(THETA-30*pi/180)/(spwidth/2)).^1000));
    Efield_opt = Efield_opt.*(1-exp(-(RHO.*cos(THETA+30*pi/180)/(spwidth/2)).^1000));
    
    
    optPhz = angle(Efield_opt);
    
    for ch = 1:vars.numWavelengths
        Epup_opt(:,:,ch) = exp(1i*optPhz*vars.lambda0/vars.lambdas(ch)).*vars.PUPIL;
    end
    
    %Get eta_s
    eta_maps = generateCouplingMap_polychromatic( Epup_opt, fiber_props, vars.lambda0, fiber.Fnum, vars.lambdas, vars.totalPower0, vars.lambdaOverD, 3*vars.lambdaOverD, vars.coords, fibmodes);
%     Xshift = zeros(inputs.numWavelengths,1);
%     Yshift = zeros(inputs.numWavelengths,1);
    eta_ss = zeros(vars.numWavelengths, 1);
    radMxs = zeros(vars.numWavelengths,1);
    eta_sX = 0;
    eta_sY = 0;
   
    %etas_offset = zeros(numWavelengths, 1);
    for ch = 1:vars.numWavelengths 
        map = eta_maps(:,:,ch); %one slice of the eta_maps cube
        map_max = max(map,[],'all'); %the maximum value in cmap
        [max_ind(1),max_ind(2)] = find(map_max==map,1); %linear coordinates of max value
        max_rho = sqrt(((vars.N/2+1)-max_ind(1))^2 + ((vars.N/2+1)-max_ind(2))^2);
    
        crp = 2*max_rho; %The length of one side of the cube to crop the image to

        cmap = map(end/2+1-floor(crp/2):end/2+1+floor(crp/2),end/2+1-floor(crp/2):end/2+1+floor(crp/2)); %the centroid
        cmap_min = min(cmap,[],'all'); %minimum value in the centroid
        [min_ind(1),min_ind(2)] = find(cmap_min==cmap); %indices of minimum value in the centroid
        min_ind = min_ind + (vars.N/2-floor(crp/2)); %adjust min values to reflect position in map, not cmap
        
        eta_sX = min_ind(2);
        eta_sY = min_ind(1);
%         Xshift(ch) = inputs.N/2-min_ind(2); %x value is wavelength, y value is offset
%         Yshift(ch) = inputs.N/2-min_ind(1);%^
        eta_ss(ch) = cmap_min;
        
        [radProf2, ~] = VFN_An_radAverage(map,[eta_sX, eta_sY]);
        radMxs(ch) = max(radProf2);
        
%         disp(radMxs);
%         disp(eta_ss);
        %etas_offset(ch) = sqrt((Xshift(ch))^2 + (Yshift(ch)^2));
    end
        
    radMx = mean(radMxs);
    eta_s = mean(eta_ss);
    
%     disp(radMx);
%     disp(eta_s);
    
%     [radProf2, ~] = VFN_An_radAverage(eta_maps,[eta_sX, eta_sY]);
%     radMx = max(radProf2);
%     disp(radMx);
%     disp(eta_s);
    
    %-- Calculate relative integration time in each frame
    relTints = eta_s./(radMx.^2);
    vars.n_calls=vars.n_calls+1;
    
    vars.relTintArray(vars.n_calls) = relTints;
end