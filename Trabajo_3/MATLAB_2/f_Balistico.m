function [u, n] = f_Balistico(z, gamma_e, beta, ue)
    
    g0 = 9.81;
    
    if z < 150e3
        rho0    = 1.225;
        zs      = 7.524e3;
    else
        rho0    = 3.875e-9;
        zs      = 59.06e3;
    end
    
    u = exp( zs*rho0*exp(-z/zs)/( 2*beta*sin(gamma_e) ) );
    
    B = zs*rho0/(sin(gamma_e));
    C = ue^2*rho/(2*beta*g0);
    n = C*exp( 2*B*exp(-z/zs) - z/zs );
    
end