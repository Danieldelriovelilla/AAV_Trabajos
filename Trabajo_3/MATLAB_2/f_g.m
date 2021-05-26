function g = f_g(z)
    
    % Constantes
    g0 = 9.81;          
    Rt = 6378e3;        
    
    % Gravedad en funcion de z
    g = g0*( Rt/(z + Rt) );

end