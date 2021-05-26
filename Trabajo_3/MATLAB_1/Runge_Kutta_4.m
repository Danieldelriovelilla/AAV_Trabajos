function [entryDinamicSol] = Runge_Kutta_4(beta, E, yi, Dt)

    % Initial conditions
    y = [yi(3); yi(2); yi(1)];
    mov = [y; 0];

    % Integrate while z is positive
    while y(3) > 0
        k1 = Ec_Dinamica(y, E, beta);
        k2 = Ec_Dinamica(y + k1*Dt/2, E, beta);
        k3 = Ec_Dinamica(y + k2*Dt/2, E, beta);
        k4 = Ec_Dinamica(y + k3*Dt, E, beta);

        yn1 = y + Dt/6*(k1 + 2*k2 + 2*k3 + k4);
        aux = [yn1; k1(2)];
        mov = [mov aux];

        y = yn1;
    end
    mov = mov(:,1:end-1);

    t       = (0:length(mov)-1)*Dt;
    gamat   = mov(1,:);
    ut      = mov(2,:);
    zt      = mov(3,:);
    at      = mov(4,:);

    % Load the solution in a struct
    entryDinamicSol.z0       = yi(1);
    entryDinamicSol.u0       = yi(2);
    entryDinamicSol.gama0    = yi(3);
    entryDinamicSol.beta0    = beta;
    entryDinamicSol.E0       = E;
    entryDinamicSol.t = t';
    entryDinamicSol.tx = max(t);
    entryDinamicSol.z = zt';
    entryDinamicSol.zx = max(zt);
    entryDinamicSol.u = ut';
    entryDinamicSol.ux = max(ut);
    entryDinamicSol.gamma = gamat';
    entryDinamicSol.gammax = max(gamat);
    entryDinamicSol.a = at';
    entryDinamicSol.ax = min(at);

end
