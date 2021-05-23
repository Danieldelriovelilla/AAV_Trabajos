function U_Ue = Analytic_Balistic(z, gamma_e, beta)
    
    zs = ones(size(z))*59.06e3;             
    zs(z<=150e3) = 7.524e3;
    
    U_Ue = exp( zs.*function_rho(z)./( 2*beta*sin(gamma_e) ));
    
end