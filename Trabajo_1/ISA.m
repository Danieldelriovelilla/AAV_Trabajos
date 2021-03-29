function [T, P, rho, a, mu, l] = ISA(z, R, gamma)

    if z > 32e3
        zi = 32e3;
        Ti = 228.65;
        lambdai = 2.8e-3;
        Pi = 0.00868e5;
    else
        zi = 20e3;
        Ti = 216.65;
        lambdai = 1.0e-3;
        Pi = 0.05474e5;
    end
    
    

    
    T = Ti + lambdai*(z - zi);
    P = Pi*(Ti/T)^(9.81/(R*lambdai));
    rho = P/(R*T);
    a = sqrt(gamma*R*T);
    mu = mu0*(T/T0)^(3/2)*(T0+S1)/(T+S1);
    
    l = mu/P*sqrt(Pi*R.*T/2);
end

