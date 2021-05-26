function g = Function_g(z)
    
    % Global variables
    global g0;
    global Rt;    
    
    % Function body
    g = g0*( Rt/(z + Rt) );

end