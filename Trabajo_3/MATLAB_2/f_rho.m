function rho = f_rho(z)

    % Parametros de la atmosfera en funcion de la altura
    if z < 150e3
        rho0    = 1.225;
        zs      = 7.524e3;
    else
        rho0    = 3.875e-9;
        zs      = 59.06e3;
    end

    % Densidad en funcion de z
    rho = rho0*exp(-z/zs);
     
end