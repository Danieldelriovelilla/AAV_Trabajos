function yP = ecuacioneF(yn, E, beta)

    Rt = 6378e3;

    gamma = yn(1);
    u = yn(2);
    z = yn(3);

    r = Rt + z;
    rhoZ = f_rho(z);
    g = f_g(z);


    Lom = rhoZ * u^2 * E /(2*beta);
    Dom = rhoZ * u^2 /(2*beta);

    gamaP   = 1/u*(Lom  - g*cos(gamma) + u^2/r*cos(gamma));
    uP      = -g*sin(gamma) - Dom;
    zP      = u*sin(gamma);


    yP = [gamaP; uP; zP];
end