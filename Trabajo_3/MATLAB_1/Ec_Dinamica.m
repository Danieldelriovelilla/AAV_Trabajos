function dydt = Ec_Dinamica(y, E, beta)
    
    % Global variables
    global Rt;
    
    % Parameters
    gamma = y(1);
    u = y(2);
    z = y(3);

    % z to r switch
    r = z + Rt;
    
    % Calculo de densidad y gravedad
    [rho, ~, ~] = Function_rho(z);
    g = Function_g(z);

    % Variables
    Lom = rho * u^2 * E /(2*beta);
    Dom = rho * u^2 /(2*beta);
    
    % Differental equations
    gamaP   = 1/u*(Lom  - g*cos(gamma) + u^2/r*cos(gamma));
    uP      = -g*sin(gamma) - Dom;
    zP      = u*sin(gamma);
    dydt = [gamaP; uP; zP];
    
end