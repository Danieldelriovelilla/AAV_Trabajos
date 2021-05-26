function yP = ecuacioneF(yn, E, beta)
    
    % Constantes
    Rt = 6378e3;

    % Valor de variables en un instante t
    gamma = yn(1);
    u = yn(2);
    z = yn(3);
    
    % Atmósfera y gravedad
    r = Rt + z;
    rhoZ = f_rho(z);
    g = f_g(z);

    % Sustentacion y resistencia
    Lom = rhoZ * u^2 * E /(2*beta);
    Dom = rhoZ * u^2 /(2*beta);

    % Ecuaciones diferenciales
    gamaP   = 1/u*(Lom  - g*cos(gamma) + u^2/r*cos(gamma));
    uP      = -g*sin(gamma) - Dom;
    zP      = u*sin(gamma);
    
    % Vector de ecuaciones diferenciales
    yP = [gamaP; uP; zP];

end