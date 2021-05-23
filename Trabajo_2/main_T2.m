clc
clear all
close all


%% Ejercicio 2

theta = linspace(0, pi,1e2+1);

cp_N = 2*(sin(theta).^2);
cp_NM = 1.839*(sin(theta).^2);
cp_NB = 2*sin(theta).^2 - 2/3*cos(theta).^2;


h = figure();
    hold on
    plot(rad2deg(theta),cp_N, '-', 'Color','k', 'LineWidth', 2, 'DisplayName', 'N')
    plot(rad2deg(theta),cp_NM, '--', 'Color','k',  'LineWidth', 2, 'DisplayName', 'N-M')
    plot(rad2deg(theta),cp_NB, ':', 'Color','k', 'LineWidth', 2, 'DisplayName', 'N-B')
    grid on; box on;
    legend('Location', 'NorthWest', 'Interpreter', 'Latex')
    xlim([0,180])
    xlabel('$\theta$ [deg]', 'Interpreter', 'Latex')
    ylabel({'$c_p$'}, 'Interpreter', 'Latex')
    Save_as_PDF(h, 'Figuras/Ejercicio_2', 'horizontal');

    
    
%% Ejercicio 3

x = linspace(0,10,1e3+1);
beta_s = -pi/2;
beta_v = pi/2;
alpha = deg2rad(0:15:90);

R = 2.7;
L = 10;
xi = R/(L)^(1/3); 


% Vector tangente:

seno = zeros(length(x),length(alpha));
cps = zeros(length(x),length(alpha));
cpv = zeros(length(x),length(alpha));
for a = 1:length(alpha)
    cps(:,a) = 2*( cos(alpha(a)) - 3/xi*x.^(2/3)*sin(beta_s)*sin(alpha(a)) ).^2./( 1 + 9/xi^2*x.^(4/3) );
    
    cpv(:,a) = 2*( cos(alpha(a)) - 3/xi*x.^(2/3)*sin(beta_v)*sin(alpha(a)) ).^2./( 1 + 9/xi^2*x.^(4/3) );
    seno(:,a) = ( cos(alpha(a)) - 3/xi*x.^(2/3)*sin(beta_v)*sin(alpha(a)) )./sqrt( 1 + 9/xi^2*x.^(4/3) );
    cpv(seno(:,a)<0,a) = 0;
end

figure()
    hold on
    plot(x,cps, 'LineWidth', 2, 'DisplayName', 'N')
    
figure()
    hold on
    plot(x,cpv, 'LineWidth', 2, 'DisplayName', 'N')


  
    
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


[x beta] = meshgrid(x,beta);

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
    
    
%{
    h = figure();
        surf(r{1}, r{2}, r{3}, 'cdata', cp,'edgecolor','none')
        colormap(jet)
        c = colorbar;
        c.TickLabelInterpreter = 'Latex';
        axis equal; box on; grid on;
        view(-60,-15)
        set(gca,'TickLabelInterpreter','latex');
        xlabel('x', 'interpreter', 'latex')
        ylabel({'y'}, 'interpreter', 'latex')
        zlabel({'z'}, 'interpreter', 'latex')
        set(get(gca,'zlabel'),'rotation',0)
        title(string(['Distribucion de cp para $\alpha = ' num2str(rad2deg(alpha(a))) '^o$']), 'interpreter', 'latex')
                
    h = figure();
        surf(r{1}, r{2}, r{3}, 'cdata', signal, 'edgecolor', 'none')
        colormap(jet)
        c = colorbar;
        c.TickLabelInterpreter = 'Latex';
        axis equal; box on; grid on;
        view(-60,-15)
        set(gca,'TickLabelInterpreter','latex');
        xlabel('x', 'interpreter', 'latex')
        ylabel({'y'}, 'interpreter', 'latex')
        zlabel({'z'}, 'interpreter', 'latex')
        set(get(gca,'zlabel'),'rotation',0)
        title(string(['Barlovento - Sotavento para $\alpha = ' num2str(rad2deg(alpha(a))) '^o$']), 'interpreter', 'latex')
%}
        
%     figure()
%         surf(x,beta,cp,'edgecolor','none')
%         colormap(jet)
%         xlabel('x')
%         ylabel('\beta')
%         zlabel('cp')
    
%     figure()
%         surf(r{1}, r{2}, r{3}, 'cdata', cpdA,'edgecolor','none')
%         colormap(jet)
%         colorbar
%         view([-60,10])
%         axis equal
%         title(['cpdA para \alpha = ' num2str(rad2deg(alpha(a)))])
%         view(-60,-15)
%         axis equal

%{
    figure()
        surf(r{1}, r{2}, r{3}, 'cdata', nm{1}.*cp,'edgecolor','none')
        colormap(jet)
        c = colorbar;
        c.TickLabelInterpreter = 'Latex';
        axis equal; box on; grid on;
        view(-60,-15)
        set(gca,'TickLabelInterpreter','latex');
        xlabel('x', 'interpreter', 'latex')
        ylabel({'y'}, 'interpreter', 'latex')
        zlabel({'z'}, 'interpreter', 'latex')
        set(get(gca,'zlabel'),'rotation',0)
        title(string(['\textbf{$\mathbf{r \times u_n|_x \cdot c_p}$ para $\mathbf{\alpha = ' num2str(rad2deg(alpha(a))) '^o}$}']), 'interpreter', 'latex')
        view(-60,-15)
        
    figure()
        surf(r{1}, r{2}, r{3}, 'cdata', nm{2}.*cp,'edgecolor','none')
        colormap(jet)
        c = colorbar;
        c.TickLabelInterpreter = 'Latex';
        axis equal; box on; grid on;
        view(-60,-15)
        set(gca,'TickLabelInterpreter','latex');
        xlabel('x', 'interpreter', 'latex')
        ylabel({'y'}, 'interpreter', 'latex')
        zlabel({'z'}, 'interpreter', 'latex')
        set(get(gca,'zlabel'),'rotation',0)
        title(string(['\textbf{$\mathbf{r \times u_n|_y \cdot c_p}$ para $\mathbf{\alpha = ' num2str(rad2deg(alpha(a))) '^o}$}']), 'interpreter', 'latex')
        view(-60,-15)
    
    figure()
        surf(r{1}, r{2}, r{3}, 'cdata', nm{3}.*cp,'edgecolor','none')
        colormap(jet)
        c = colorbar;
        c.TickLabelInterpreter = 'Latex';
        axis equal; box on; grid on;
        view(-60,-15)
        set(gca,'TickLabelInterpreter','latex');
        xlabel('x', 'interpreter', 'latex')
        ylabel({'y'}, 'interpreter', 'latex')
        zlabel({'z'}, 'interpreter', 'latex')
        set(get(gca,'zlabel'),'rotation',0)
        title(string(['\textbf{$\mathbf{r \times u_n|_z \cdot c_p}$ para $\mathbf{\alpha = ' num2str(rad2deg(alpha(a))) '^o}$}']), 'interpreter', 'latex')
        view(-60,-15)
%}
end

    
%%{
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
%%}





    
%%

%{
syms x beta  v

% Create the symbolic vector
v = sym([1,3]);
v(1) = x;
v(2) = xi*x^(1/3)*cos(beta);
v(3) = xi*x^(1/3)*sin(beta);

surface = double( subs(v, [x, beta], [0 0]) );

% Plot surface
fsurf(v(1), v(2), v(3), [0 10 0 2*pi])
% v(x, beta) = [x, xi*x^(1/3)*cos(beta), xi*x^(1/3)*sin(beta)];

% Tabgent vectors
tx = diff(v,x);
tbeta = diff(v,beta);

% Normal vector
n = cross(tx, tbeta);
n = n/norm(n);

% Velocidad con angulo de ataque

alpha = 0;

u = sym([1,3]);
u(1) = cos(alpha);
u(2) = 0;
u(3) = sin(alpha);

cp = 2*( ( cos(alpha(1)) - 3/xi*x^(2/3)*sin(alpha(1))*sin(beta) )/( 1 + 9/xi^2*x^(4/3) ) )^2;
figure()
    fsurf(cp, [0 10 0 2*pi])

dA = xi*sqrt( 1 + 9/xi^2*x^(4/3) );
figure()
    fsurf(cp*dA, [0 10 0 2*pi])  

%}
    
%%





    %view(0, 0)

% cp = @(x, beta) ((10618401823952946614204613221403.*cos(beta).^2)./(20282409603651670423947251286016.*x.^(1/3)) + ...
%                 (10618401823952946614204613221403.*sin(beta).^2)./(20282409603651670423947251286016.*x.^(1/3)))./((31855205471858839842613839664209.*abs(cos(beta)).^2.*abs(x).^(2/3))./20282409603651670423947251286016 + ...
%                 (31855205471858839842613839664209.*abs(sin(beta)).^2.*abs(x).^(2/3))./20282409603651670423947251286016 +...
%                 abs((10618401823952946614204613221403.*cos(beta).^2)./(20282409603651670423947251286016.*x.^(1/3)) + ...
%                 (10618401823952946614204613221403*sin(beta).^2)./(20282409603651670423947251286016.*x.^(1/3))).^2).^(1/2);
% res = integral2(cp,0,10,0,2*pi);
%cpi = int(cp);
%integral2(cpi,0,10,0,2*pi)


%{
dA = sqrt(  norm(jacobian([v(1), v(2)], [x, beta]))^2 +...
            norm(jacobian([v(2), v(3)], [x, beta]))^2 +...
            norm(jacobian([v(1), v(3)], [x, beta]))^2);
        
        
cp = @(x, beta) 2*( 1 ).^2./( 1 + 9/xi^2*x.^(4/3) ); 

cn = int(cp*dA)
%}