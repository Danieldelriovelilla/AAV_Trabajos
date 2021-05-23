function U_Ue = Analytic_Balistic(z, gamma_e, beta)
    
    zs = ones(size(z))*59.06e3;
    zs(z<=150e3) = 7.524e3;
    
    rho0 = ones(size(z))*3.875e-9;
    rho0(z<=150e3) = 1.225;
    
    U_Ue = exp( zs.*rho0.*exp(-z./zs) ./ ( 2*beta*sin(gamma_e) ));
    
end