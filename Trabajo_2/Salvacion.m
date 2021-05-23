%% INTEGRAL HECHA A MANO

R = 2.7;
L = 10;

xi = R/(L)^(1/3); 
alpha = deg2rad([0:1:90]);
a = 10;

Nx = 1e3+1;
Nbeta = 1e2+1;
x = linspace(0,10,Nx);
beta = linspace(0,2*pi,Nbeta);

x_i = linspace(0,10,Nx);
beta_i = linspace(0,2*pi,Nbeta);


[x, beta] = meshgrid(x,beta);

r{1} = x;
r{2} = xi*x.^(1/3).*cos(beta);
r{3} = xi*x.^(1/3).*sin(beta);

n{1} = 1./sqrt( 1 + 9/xi^2*x.^(4/3) );
n{2} = -3./xi.*x.^(2/3).*cos(beta)./sqrt( 1 + 9/xi^2*x.^(4/3) );
n{3} = -3./xi.*x.^(2/3).*sin(beta)./sqrt( 1 + 9/xi^2*x.^(4/3) );

nm{1} = xi*x.^(1/3).*cos(beta).*(-3./xi.*x.^(2/3).*sin(beta)./sqrt( 1 + 9/xi^2*x.^(4/3) )) - ...
        xi*x.^(1/3).*sin(beta).*(-3./xi.*x.^(2/3).*cos(beta)./sqrt( 1 + 9/xi^2*x.^(4/3) ));
nm{2} = (-3./xi.*x.^(2/3).*sin(beta)./sqrt( 1 + 9/xi^2*x.^(4/3) )).*1./sqrt( 1 + 9/xi^2*x.^(4/3) ) - ...
        x.*(-3./xi.*x.^(2/3).*sin(beta)./sqrt( 1 + 9/xi^2*x.^(4/3) ));
nm{3} = x.*(-3./xi.*x.^(2/3).*cos(beta)./sqrt( 1 + 9/xi^2*x.^(4/3) )) - ...
        xi*x.^(1/3).*cos(beta).*(1./sqrt( 1 + 9/xi^2*x.^(4/3) ));


close all
for a = 1:length(alpha)
    cp = 2*( ( cos(alpha(a)) - 3/xi*x.^(2/3).*sin(alpha(a)).*sin(beta) )./sqrt( 1 + 9/xi^2*x.^(4/3) ) ).^2;
    seno = ( cos(alpha(a)) - 3/xi*x.^(2/3).*sin(beta)*sin(alpha(a)) )./sqrt( 1 + 9/xi^2.*x.^(4/3) );
    signal = ones(size(seno));
    signal(seno<0) = 0;
    cp = cp.*signal;
    dA = xi.*sqrt( 1 + 9./xi^2*x.^(4/3) );
    
    
    cpdA = cp.*dA;
    mcpcdA{1} = nm{1}.*cp.*dA;
    mcpcdA{2} = nm{2}.*cp.*dA;
    mcpcdA{3} = nm{3}.*cp.*dA;

    cn(a) = trapz(beta_i,trapz(x_i,cpdA,2))/(pi*R^2);

    cm_1(a) = trapz(beta_i,trapz(x_i,mcpcdA{1},2))/(pi*R^2*L);
    cm_2(a) = trapz(beta_i,trapz(x_i,mcpcdA{2},2))/(pi*R^2*L);
    cm_3(a) = trapz(beta_i,trapz(x_i,mcpcdA{3},2))/(pi*R^2*L);
end

h = figure();
    hold on
    plot(rad2deg(alpha), cn, '-', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$Cn$')
    plot(rad2deg(alpha), cm_1, '-.','Color', 'k', 'LineWidth', 1.5, 'DisplayName', '$Cm_x$')
    plot(rad2deg(alpha), cm_2, '.', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$Cm_y$')
    plot(rad2deg(alpha), cm_3, ':', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$Cm_z$')
    grid on; box on
    legend('Location', 'NorthWest', 'Interpreter', 'Latex')
    xlabel('$\alpha$ [deg]', 'Interpreter', 'Latex')
    xlim([0,90])
    
    set(gca,'TickLabelInterpreter','latex');
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 pos(3), pos(4)])
    print(h, '-dpng', ['Figuras/Coeficientes aerodinámicos.png'],'-r750','-painters')