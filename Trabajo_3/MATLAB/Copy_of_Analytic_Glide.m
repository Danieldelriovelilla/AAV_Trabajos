function U_Ue = Analytic_Glide(z, beta, E)
    
    RT = 6378e3;                            % [km]
    
    zs = ones(size(z))*59.06e3;             
    zs(z<=150e3) = 7.524e3;
    
    rho0 = ones(size(z))*3.875e-9;          % [kg/m^3]
    rho0(z<=150e3) = 1.225;                 % [kg/m^3]
       
    U_Ue = 1./sqrt( 1 + E.*rho0*RT.*exp(-z./zs)./( 2*beta ));
    
end