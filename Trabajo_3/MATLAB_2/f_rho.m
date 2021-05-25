function rho = f_rho(z)

    if z < 150e3
        rho0    = 1.225;
        zs      = 7.524e3;
    else
        rho0    = 3.875e-9;
        zs      = 59.06e3;
    end

     rho = rho0*exp(-z/zs);
     
end