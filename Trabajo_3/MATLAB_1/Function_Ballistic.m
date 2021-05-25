function [U_Ue, n] = Function_Ballistic(z, gamma_e, beta, Ue)
    
    % Gloval variables
    global g0
    
    % Function body
    [~, zs, rho0] = Function_rho(z);
    
    B = zs.*rho0./(2*beta*sin(gamma_e));
    C = Ue^2.*rho0./(2*beta*g0);
    
    U_Ue = exp( B.*exp(-z./zs) );
    n = C.*exp( 2*B.*exp(-z./zs) - z./zs );
    
end