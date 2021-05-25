clc
clear
close all


%% CONSTANTES

E = [0 0.5 1];
beta = [50 280 500];
ze = 100e3;
ue = 7.5e3;
gammae = -5*pi/180;
y0_set = [[150e3; 7.5e3; deg2rad(-5)], [150e3; 8.5e3; deg2rad(-5)],...
          [150e3; 7.5e3; deg2rad(-10)]];

      
%% ANALITICO

[ub, nb, us, ns] = f_Analitico(ue, gammae, ze, beta, E);


%% ECUACIONES NUMERICAS

% Intial condition
ci = [ze; ue; gammae];

% RK4 time step
Dt = 0.1;

for i = 1:length(beta)
    for j = 1:length(E)
        solucion{i,j} = RK4(beta(i), E(j), ci, Dt);   
    end
end
            

%% PLOTS

colors = [0, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.4660, 0.6740, 0.1880;
          0.4940, 0.1840, 0.5560];  

line_style = {'-','--',':'};

for i = 1:length(beta)
    for j = 1:length(E)
        
%{
        figure(1)
            hold on
            plot(solucion(b,e).u/solucion(b,e).u0,...
                 solucion(b,e).z/1000, line_style{b},...
                 'Color', colors(e,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str))
            axis([0, 1, 0, 122])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
            xlabel('$U/U_0$', 'Interpreter', 'Latex')
            ylabel('$z$ [km]', 'Interpreter', 'Latex')
             
        figure(2)
            hold on
            plot(rad2deg(solucion(b,e).gamma),...
                 solucion(b,e).z/1000,line_style{b},...
                 'Color', colors(e,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str))
            axis([-90, 4.75, 0, 122])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
            xlabel('$\gamma$', 'Interpreter', 'Latex')
            ylabel('$z$ [km]', 'Interpreter', 'Latex')
        
        figure(3)
            hold on
            plot(-solucion(b,e).a/9.81,...
                 solucion(b,e).z/1000,line_style{b},...
                 'Color', colors(e,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str))
            axis([0, 15, 0, 121])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
            xlabel('$a$', 'Interpreter', 'Latex')
            ylabel('$z$ [km]', 'Interpreter', 'Latex')
%}
        h{1} = figure(1);
            hold on
            plot(solucion{i,j}(:,1), rad2deg(solucion{i,j}(:,end)),...
                 line_style{i}, 'Color', colors(j,:), 'LineWidth', 1.5,...
                 'DisplayName', ['$\beta = ' num2str(beta(i)) ' $ ; $E = ' num2str(E(j)) '$'])
            axis([0, 1996, -90, 4.2])
            grid on; box on;
            h{1}.Position = [680   558   1.4*560   420];
            legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)
            xlabel('$t$ [s]', 'Interpreter', 'Latex')
            ylabel({'$\gamma$';'[$\mathrm{^o}$]'}, 'Interpreter', 'Latex')
            
        h{2} = figure(2);
            hold on
            plot(solucion{i,j}(:,1), solucion{i,j}(:,3)/ue,...
                line_style{i}, 'Color', colors(j,:), 'LineWidth', 1.5,...
                'DisplayName', ['$\beta = ' num2str(beta(i)) ' $ ; $E = ' num2str(E(j)) '$'])
            axis([0, 1996, 0, 1])
            grid on; box on;
            h{2}.Position = [680   558   1.4*560   420];
            legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)
            xlabel('$t$ [s]', 'Interpreter', 'Latex')
            ylabel({'$U/U_0$'}, 'Interpreter', 'Latex')

        h{3} = figure(3);
            hold on
            plot(solucion{i,j}(:,1),solucion{i,j}(:,2)/1000,...
                line_style{i}, 'Color', colors(j,:), 'LineWidth', 1.5,...
                'DisplayName', ['$\beta = ' num2str(beta(i)) ' $ ; $E = ' num2str(E(j)) '$'])
            axis([0, 1996, 0, 121])
            grid on; box on;
            h{3}.Position = [680   558   1.4*560   420];
            legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)
            xlabel('$t$ [s]', 'Interpreter', 'Latex')
            ylabel({'$z$';'[km]'}, 'Interpreter', 'Latex')
            
        h{4} = figure(4);
            hold on
            plot(solucion{i,j}(:,1),-solucion{i,j}(:,4)/9.81,...
                line_style{i}, 'Color', colors(j,:), 'LineWidth', 1.5,...
                'DisplayName', ['$\beta = ' num2str(beta(i)) ' $ ; $E = ' num2str(E(j)) '$'])
            axis([0, 1996, 0, 14.9])
            grid on; box on;
            h{4}.Position = [680   558   1.4*560   420];
            legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)
            xlabel('$t$ [s]', 'Interpreter', 'Latex')
            ylabel({'$a$';'[g]'}, 'Interpreter', 'Latex')
            

    end
end

figure(1);
Save_as_PDF(h{1}, 'Figures\test1', 'horizontal')

figure(2);
Save_as_PDF(h{2}, 'Figures\test2', 'horizontal')

figure(3);
Save_as_PDF(h{3}, 'Figures\test3', 'horizontal')

figure(4);
Save_as_PDF(h{4}, 'Figures\test4', 'horizontal')

close all

%% Balistico - Analitico
z = linspace(0,200e3,1e3+1);
gamma_e = deg2rad(-1);
beta = 150;                                 % [kg/m^2]

U_Balistic = Analytic_Balistic(z, gamma_e, beta);

h = figure();
    plot(U_Balistic, z);
    xlabel('$U/U_e$', 'Interpreter', 'Latex')
    ylabel('$z$ [m]', 'Interpreter', 'Latex')
    title('Analitica: balistica', 'Interpreter', 'Latex')
    

%% Planeo - Analaitico
z = linspace(0,100e3,1e3+1);
beta = 500;                                 % [kg/m^2]
E = 0.25;

U_Glide = Analytic_Glide(z, beta, E);

h = figure();
    plot(U_Glide, z);
    xlabel('$U/U_e$', 'Interpreter', 'Latex')
    ylabel('$z$ [m]', 'Interpreter', 'Latex')
    title('Analitica: planeo', 'Interpreter', 'Latex')