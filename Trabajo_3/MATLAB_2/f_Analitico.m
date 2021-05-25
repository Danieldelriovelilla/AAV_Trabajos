function [ub, nb, us, ns] = f_Analitico(ue, gammae, ze, beta, E)

    g0 = 9.81;
    Rt = 6378e3;

    z = linspace(0,ze, 1e3+1);

    rho0 = zeros(size(z)) + 3.875e-9;
    zs = zeros(size(z)) + 59.06e3;
    idx = find(z < 150e3);
    rho0(idx) = 1.225;
    zs(idx) = 7.524e3;

    for i = 1:length(beta)
        % Balistico
        B = zs.*rho0./(sin(gammae));
        C = ue^2*rho0./(2*beta(i)*g0);
        ub{i} = ( exp( zs.*rho0.*exp(-z./zs)./( 2*beta(i)*sin(gammae) ) ) )';
        nb{i} = ( C.*exp( 2*B.*exp(-z./zs) - z./zs ) )';

        for j = 2:length(E)
            % Sustentado
            us{i,j-1} = ( ( 1 + E(j)*rho0.*Rt.*exp(-z./zs)/( 2*beta(i) )).^(-1/2) )';
            ns{i,j-1} = ( ( E(j) + 2*beta(i)*exp(z./zs)./(rho0*Rt) ).^(-1) )';
        end
    end

end