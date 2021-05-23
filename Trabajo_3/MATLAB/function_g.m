function g = function_g(z)
    
    global g0;
    global Rt;    
    
    g = g0*( Rt/(z + Rt) );

end