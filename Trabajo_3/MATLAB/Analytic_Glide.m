function U_Ue = Analytic_Glide(z, beta, E)
    
    global Rt
       
    U_Ue = 1./sqrt( 1 + E.*Rt*function_rho(z)./( 2*beta ));
    
end