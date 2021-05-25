function [u, n] = f_Sustentado(z, beta, E)
    
    Rt = 6378e3;
    
    if z < 150e3
        rho0    = 1.225;
        zs      = 7.524e3;
    else
        rho0    = 3.875e-9;
        zs      = 59.06e3;
    end
        
    u = ( 1 + E*rho0*Rt*exp(-z/zs)/( 2*beta ))^(-1/2);
    n = ( E + 2*beta*exp(z/zs)/(rho0*Rt) )^(-1);
    
end