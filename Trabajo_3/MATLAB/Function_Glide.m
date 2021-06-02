function [U_Ue, n] = Function_Glide(z, beta, E)
    
    % Gloval variables
    global Rt
    
    % Function body
    [rho, zs, rho0] = Function_rho(z);
    
    U_Ue = 1./sqrt( 1 + E.*Rt.*rho./( 2*beta ));
    n = 1./( E + 2*beta*exp(z./zs)./(rho0*Rt) );
    
end