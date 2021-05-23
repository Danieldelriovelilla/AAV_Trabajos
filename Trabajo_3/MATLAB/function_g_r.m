function g = function_g_r(r)

    g0 = 9.81;          % [m/s^2]
    Rt = 6378e3;        % [m]
    z = r - Rt;
    
    g = g0*( Rt/(z + Rt) );

end