function [tilt_valsCopy, clocking, I] = simple_ADCOPT(inpar)

%% Gridscan optimization
    inpar.lam0OverD = inpar.lambdas(ceil(inpar.numWavelengths / 2)) / inpar.keckD;

    dz = inpar.pyyCopy*inpar.lam0OverD;

    % dz = pyyCopy*inpar.lam0OverD_meters; %MFT Value

    dz = dz*890.16;
    dz = py.tuple(dz);
    dz = py.numpy.array(dz);

    samp = linspace(0,90,100);
    tilt_vals = [];

    for i = 1:100
    
        clocking = deg2rad(samp(i));
    %     disp(clocking);
        tilt = py.triple_prism.triple_prism(dz,inpar.n0,inpar.n1,inpar.n2,inpar.n3,inpar.phi1_ADC,inpar.phi2_ADC,inpar.phi3_ADC, clocking,0);
    
        tilt = tilt.tolist();
    
        for j = 1:inpar.numWavelengths
            tilt_out(j) = tilt{j};
        end
    
        tilt_out = tilt_out - tilt_out(ceil(inpar.numWavelengths/2));
        %disp(tilt_out);
        tilt_out = tilt_out/890.16;
    
        % Wavelength dependent dispersion
        tilt_vals(:,i) = tilt_out';
    %     for j = 1:11
    %         tilt_vals(:,i) = tilt_out;
    %     end
    
        % Sum of dispersion at each wavelength squared
        tiltsums(i) = sum(tilt_out.^2);
    
        %Slope of ADC induced tilt
        tiltslopes(i) = tilt_out(end) - tilt_out(1);
    
    end

    slopes = tilt_vals(2,:) - tilt_vals(1,:);

%     figure()
%     hold on
%     title(['Gridsearch Min']);
%     xlabel(['Clocking Angle 0-90 deg']);
%     ylabel(['Sum of Squares']);
%     % plot(samp, tilt1,'Color', 'r');
%     % plot(samp, tilt2,'Color', 'y');
%     % plot(samp, tilt3,'Color', 'g');
%     % plot(samp, tilt4,'Color', 'c');
%     % plot(samp, tilt5,'Color', 'm');
%     plot(samp, tiltsums, 'Color', 'b');
%     hold off

    [M,I] = min(tiltsums);
    clocking = samp(I);
    disp(M);
    disp(I);
    disp(samp(I));
    disp(tilt_vals(:,I));
    
    tilt_valsCopy = tilt_vals;

%     figure()
%     hold on
%     title(['Polychromatic Tilt 2 - 2.4 um']);
%     xlabel(['Clocking Angle 0 - 90 deg']);
%     ylabel(['Tilt (radians)']);
%     plot(samp, tilt_vals,'Color','r');

% for i = 1:inpar.numWavelengths
%     plot(samp, tilt_vals(i,:));
% end

% plot(samp, tilt_vals,'Color', 'r');
% plot(samp, tilt2,'Color', 'y');
% plot(samp, tilt3,'Color', 'g');
% plot(samp, tilt4,'Color', 'c');
% plot(samp, tilt5,'Color', 'm');
% plot(samp, tilt6,'Color', 'r');
% plot(samp, tilt7,'Color', 'y');
% plot(samp, tilt8,'Color', 'g');
% plot(samp, tilt9,'Color', 'c');
% plot(samp, tilt10,'Color', 'm');
% plot(samp, tilt11,'Color', 'r');
% plot(samp, tiltsums, 'Color', 'b');
% hold off
% 
% figure()
% hold on
% title(['Slope of ADC tilt vs increasing clocking angle']);
% xlabel(['Clocking Angle 0 - 90 deg']);
% ylabel(['ADC Dispersion Slope (tilt @ 2.1um - tilt @ 2.0um)']);
% plot(samp, slopes,'Color','b');
% hold off

% figure()
% hold on
% title(['Dispersion across band for clocking angles 0-90deg']);
% xlabel(['Clocking Angle 0 - 90 deg']);
% ylabel(['Dispersion across band (mas) (tilt @ 2.4um - tilt @ 2.0um)']);
% 
% plot(samp, tiltslopes * (3600*1000 *180)/pi,'Color','b');
% 
% hold off

end