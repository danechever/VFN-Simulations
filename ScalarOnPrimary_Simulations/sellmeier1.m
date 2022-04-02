function n = sellmeier1(wvs,temp,K1, K2, K3, L1, L2, L3)

    temp = 273.15 + 3; %Kelvin
    
    Klen = size(K1);
    Llen = size(K1);
    
    wvs2 = wvs^2;
    
    for i = 1:Klen
        K1_tot = sum(K1(i) * temp^(i-1));
        K2_tot = sum(K2(i) * temp^(i-1));
        K3_tot = sum(K3(i) * temp^(i-1));
    end
    
    %disp([K1_tot K2_tot K3_tot]);

    for i = 1:Klen
        L1_tot = sum(L1(i) * temp^(i-1));
        L2_tot = sum(L2(i) * temp^(i-1));
        L3_tot = sum(L3(i) * temp^(i-1));
    end
    
    %disp([K1_tot K2_tot K3_tot]);

    
    arg1 = K1_tot*wvs2/(wvs2 - L1_tot^2);
    arg2 = K2_tot*wvs2/(wvs2 - L2_tot^2);
    arg3 = K3_tot*wvs2/(wvs2 - L3_tot^2);
    
    n = sqrt(arg1 +arg2 + arg3 + 1);
    
end