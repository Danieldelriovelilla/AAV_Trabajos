function numeric = Function_Numeric(beta, E, y0, Dt)
    
    for b = 1:length(beta)
        for e = 1:length(E)
            numeric(b,e) = Runge_Kutta_4(beta(b), E(e), y0, Dt);   
        end
    end

end