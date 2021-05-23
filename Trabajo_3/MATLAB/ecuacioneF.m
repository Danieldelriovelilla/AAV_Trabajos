function yP = ecuacioneF(yn, E, beta)

    global Rt;

    gamma = yn(1);
    u = yn(2);
    z = yn(3);

    r = z + Rt;

    rhoZ = function_rho_r(r);
    g = function_g_r(r);


    Lom = rhoZ * u^2 * E /(2*beta);
    Dom = rhoZ * u^2 /(2*beta);

    gamaP   = 1/u*(Lom  - g*cos(gamma) + u^2/r*cos(gamma));
    uP      = -g*sin(gamma) - Dom;
    zP      = u*sin(gamma);


    yP = [gamaP; uP; zP];
end