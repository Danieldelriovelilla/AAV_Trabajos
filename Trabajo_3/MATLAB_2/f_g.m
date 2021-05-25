function g = f_g(z)

    g0 = 9.81;          % [m/s^2]
    Rt = 6378e3;        % [m]
    
    g = g0*( Rt/(z + Rt) );

end