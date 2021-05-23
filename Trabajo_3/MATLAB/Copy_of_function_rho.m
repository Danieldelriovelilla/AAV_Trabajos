function rhoOut = function_rho(z)

    if z < 150e3
        rho0    = 1.225;
        zs      = 7.524e3;
    else
        rho0    = 3.875e-9;
        zs      = 59.06e3;
    end

     rhoOut = rho0*exp(-z/zs);
     
end
    zs = ones(size(z))*59.06e3;             
    zs(z<=150e3) = 7.524e3;
    
    rho0 = ones(size(z))*3.875e-9;          % [kg/m^3]
    rho0(z<=150e3) = 1.225;                 % [kg/m^3]

    rhoOut = rho0.*exp(-z./zs);
    