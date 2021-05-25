function [solucion] = RK4(beta, E, yi, Dt)

    yn = [yi(3); yi(2); yi(1)];
    mov = [yn; 0];

    while yn(3) > 0

        k1 = ecuacioneF(yn, E, beta);
        k2 = ecuacioneF(yn + k1*Dt/2, E, beta);
        k3 = ecuacioneF(yn + k2*Dt/2, E, beta);
        k4 = ecuacioneF(yn + k3*Dt, E, beta);

        yn1 = yn + Dt/6*(k1 + 2*k2 + 2*k3 + k4);
        aux = [yn1; k1(2)];
        mov = [mov aux];
        yn = yn1;

    end
    
    mov = mov(:,1:end-1);

    t       = ( (0:length(mov)-1)*Dt )';
    gamat   = ( mov(1,:) )';
    ut      = ( mov(2,:) )';
    zt      = ( mov(3,:) )';
    at      = ( mov(4,:) )';

    solucion = [t zt ut at gamat];

end
