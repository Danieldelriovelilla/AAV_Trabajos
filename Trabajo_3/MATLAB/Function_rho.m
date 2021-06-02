function [rho, zs, rho0] = Function_rho(z)
    
    % Function body
    zs = ones(size(z))*59.06e3;             
    zs(z<=150e3) = 7.524e3;
    
    rho0 = ones(size(z))*3.875e-9;          % [kg/m^3]
    rho0(z<=150e3) = 1.225;                 % [kg/m^3]

    rho = rho0.*exp(-z./zs);
    
end
