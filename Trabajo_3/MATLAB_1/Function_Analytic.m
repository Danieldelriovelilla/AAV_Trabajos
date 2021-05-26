function [ballistic, glide] = Function_Analytic(beta, E, y0)

    % Hight array
    z = linspace(0,y0(1),1e4+1);            % [m]
    
    % Evaluate analytic equations
    for b = 1:length(beta)
        
        % Ballistic
        ballistic(b).beta = beta(b);
        ballistic(b).z = z;
        [ballistic(b).u, ballistic(b).n] = Function_Ballistic(z, y0(3), beta(b), y0(2));
        
        % Glide
        for e = 2:length(E)
            glide(b,e).beta = beta(b);
            glide(b,e).E = E(e);
            glide(b,e).z = z;
            [glide(b,e).u, glide(b,e).n] = Analytic_Glide(z, beta(b), E(e));
        end
        
    end


end