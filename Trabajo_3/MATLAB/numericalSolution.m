clc
clear 
close all


%% Global variables

global g0; g0 = 9.81;
global Rt; Rt = 6378e3;



%% PARAMETERS & IC

% Parameters
E = [0 0.5 1];
beta = [50 280 500];

% Intial condition
yi = [100e3; 7.5e3; -5*pi/180];

% RK4 time step
Dt = 0.1;



%% INTEGRATE EQUATIONS

for b = 1:length(beta)
    for e = 1:length(E)
        entryDinamicSol(b,e) = Runge_Kutta_4(beta(b), E(e), yi, Dt);   
    end
end
            


%% PLOTS

colors = [0, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.4660, 0.6740, 0.1880;
          0.4940, 0.1840, 0.5560];  

line_style = {'-','--',':'};
leg_str = "['$E = ' num2str(E(e)) ' $ ; $\beta = ' num2str(beta(b)) '$']";

for b = 1:length(beta)
    for e = 1:length(E)
        figure(1)
            hold on
            plot(entryDinamicSol(b,e).u/entryDinamicSol(b,e).u0,...
                 entryDinamicSol(b,e).z/1000, line_style{b},...
                 'Color', colors(e,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str))
            axis([0, 1, 0, 122])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
            xlabel('$U/U_0$', 'Interpreter', 'Latex')
            ylabel('$z$ [km]', 'Interpreter', 'Latex')
             
        figure(2)
            hold on
            plot(rad2deg(entryDinamicSol(b,e).gamma),...
                 entryDinamicSol(b,e).z/1000,line_style{b},...
                 'Color', colors(e,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str))
            axis([-90, 4.75, 0, 122])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
            xlabel('$\gamma$', 'Interpreter', 'Latex')
            ylabel('$z$ [km]', 'Interpreter', 'Latex')
        
        figure(3)
            hold on
            plot(-entryDinamicSol(b,e).a/9.81,...
                 entryDinamicSol(b,e).z/1000,line_style{b},...
                 'Color', colors(e,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str))
            axis([0, 15, 0, 121])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
            xlabel('$a$', 'Interpreter', 'Latex')
            ylabel('$z$ [km]', 'Interpreter', 'Latex')
        
        figure(4)
            hold on
            plot(entryDinamicSol(b,e).t/60, rad2deg(entryDinamicSol(b,e).gamma),...
                 line_style{b}, 'Color', colors(e,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str))
            axis([0, 33.17, -90, 4.2])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
            xlabel('$t$ [s]', 'Interpreter', 'Latex')
            ylabel('$\gamma$ [$\mathrm{^o}$]', 'Interpreter', 'Latex')
            
        figure(5)
            hold on
            plot(entryDinamicSol(b,e).t/60,...
                entryDinamicSol(b,e).u/entryDinamicSol(b,e).u0,...
                line_style{b}, 'Color', colors(e,:), 'LineWidth', 1.5,...
                'DisplayName', eval(leg_str))
            axis([0, 33.23, 0, 1])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
            xlabel('$t$ [s]', 'Interpreter', 'Latex')
            ylabel('$U/U_0$', 'Interpreter', 'Latex')

        figure(6)
            hold on
            plot(entryDinamicSol(b,e).t/60,entryDinamicSol(b,e).z/1000,...
                line_style{b}, 'Color', colors(e,:), 'LineWidth', 1.5,...
                'DisplayName', eval(leg_str))
            axis([0, 33.17, 0, 122])
            grid on; box on;
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)%, 'Location', 'EastOutside')
            xlabel('$t$ [s]', 'Interpreter', 'Latex')
            ylabel('$z$ [km]', 'Interpreter', 'Latex')

    end
end
