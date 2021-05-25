function dydt = dinamica(t,y)

dydt = [y(3)*sin(y(2)) - function_rho_r(y(3))/(2*beta)*y(1)^2;
        ( y(1)^2*(function_rho_r(y(3))*E)/(2*beta) - cos(y(2))/y(3) )/y(1) - function_g_r(y(3))*cos(y(2));
        y(1)*cos(y(2))];