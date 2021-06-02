clc
clear
close all


%% GLOBAL VARIABLES

global g0; g0 = 9.81;               % [m/s^2]
global Rt; Rt = 6378e3;             % [m]


%% PARAMETERS & IC

% Capsule parameters
E = [0 0.5 1];
beta = [50 250 500];
% Intial condition
y0_set = [[100e3, 7.5e3, deg2rad(-5)]',...
          [100e3, 8.5e3, deg2rad(-5)]',...
          [100e3, 7.5e3, deg2rad(-10)]'];

% RK4 time step
Dt = 0.1;                           % [s]


%% PLOT DEFINITION

colors = [0, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980;
          0.4660, 0.6740, 0.1880;
          0.4940, 0.1840, 0.5560]; 

line_style = {'-','-.',':'};
leg_str = "['$E = ' num2str(E(e)) ' $ ; $\beta = ' num2str(beta(b)) '$']";
PDF_str = "['Figures/Num_Fig_' num2str(f) '_ue_',  num2str(y0(2)/100), '_gammae_', num2str(rad2deg(y0(3)))]";      

      
for y = 1:size(y0_set,2)
    
    y0 = y0_set(:,y);

    % Numeric solution
    numeric = Function_Numeric(beta, E, y0, Dt);
    numerics{y} = numeric;
    
    % Analytic solution
    [ballistic, glide] = Function_Analytic(beta, E, y0);
    analytics(y).ballistic = ballistic;
    analytics(y).glide = glide;
    

    % Numeric
    for b = 1:length(beta)
        for e = 1:length(E)

            h{1} = figure(1);
                hold on
                plot(numeric(b,e).t/60, rad2deg(numeric(b,e).gamma),...
                     line_style{b}, 'Color', colors(e,:), 'LineWidth', 1.5,...
                     'DisplayName', eval(leg_str))
                axis([0, max([numeric(:,:).tx])/60, -90, rad2deg(max([numeric(:,:).gammax]))])
                grid on; box on;
                xlabel('$t$ [min]', 'Interpreter', 'Latex')
                ylabel({'$\gamma$';'[$\mathrm{^o}$]'}, 'Interpreter', 'Latex')

            h{2} = figure(2);
                hold on
                plot(numeric(b,e).t/60,...
                    numeric(b,e).u/numeric(b,e).u0,...
                    line_style{b}, 'Color', colors(e,:), 'LineWidth', 1.5,...
                    'DisplayName', eval(leg_str))
                axis([0, max([numeric(:,:).tx])/60, 0, max([numeric(:,:).ux])/numeric(b,e).u0])
                grid on; box on;
                xlabel('$t$ [min]', 'Interpreter', 'Latex')
                ylabel({'$\frac{U}{U_0}$'}, 'Interpreter', 'Latex')

            h{3} = figure(3);
                hold on
                plot(numeric(b,e).t/60,numeric(b,e).z/1000,...
                    line_style{b}, 'Color', colors(e,:), 'LineWidth', 1.5,...
                    'DisplayName', eval(leg_str))
                axis([0, max([numeric(:,:).tx])/60, 0, max([numeric(:,:).zx])/1000])
                grid on; box on;
                xlabel('$t$ [min]', 'Interpreter', 'Latex')
                ylabel({'$z$'; '[km]'}, 'Interpreter', 'Latex')
                
            h{4} = figure(4);
                hold on
                plot(numeric(b,e).t/60,-numeric(b,e).a/9.81,line_style{b},...
                     'Color', colors(e,:), 'LineWidth', 1.5,...
                     'DisplayName', eval(leg_str))
                axis([0, max([numeric(:,:).tx])/60,-.1, -min([numeric(:,:).ax])/9.81])
                grid on; box on;
                xlabel('$t$ [min]', 'Interpreter', 'Latex')
                ylabel({'$n$';'[g]'}, 'Interpreter', 'Latex')
                
        end
    end

    for f = 1:4
        Save_as_PDF(h{f}, eval(PDF_str), 'vertical')
        figure(f);
            legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 3)
        h{f} = figure(f);
        h{f}.Position = [7.0729, 5.8021, 1.15*5.8333, 1.15*4.3750];
        Save_as_PDF(h{f}, eval(PDF_str), 'horizontal')
    end
    
    close all
    clear h

end


%% ANALYTIC AND NUMERIC COMPARISON

close all
line_style = {'-',':','--'};

% BALLISTIC
e = 1;
leg_str_ba = "['An. $\gamma_e =' num2str(rad2deg(y0_set(3,y))) '$, $\beta =' num2str(beta(b)) '$']";
leg_str_bn = "['Nu. $\gamma_e =' num2str(rad2deg(y0_set(3,y))) '$, $\beta =' num2str(beta(b)) '$']";

% Speed
h = figure(1);
    hold on
    for y = 1:2:3
        for b = 1:2:3
            plot(analytics(y).ballistic(b).u, analytics(y).ballistic(b).z/1000,...
            line_style{b}, 'Color', colors(y,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str_ba));
            plot(numerics{y}(b,e).u/numerics{y}(b,e).u0,...
                 numerics{y}(b,e).z/1000, line_style{b},...
                 'Color', colors(y+1,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str_bn))
        end
    end
    axis([0, max([numerics{y}(:,:).ux])/numerics{y}(b,e).u0, 0, max([numerics{y}(:,:).zx])/1000])
    grid on; box on;
    legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 2)
    xlabel('$U/U_e$', 'Interpreter', 'Latex')
    ylabel({'$z$';'[km]'}, 'Interpreter', 'Latex')
    
    Save_as_PDF(h, 'Figures/Comparison_Ballistic_U', 'vertical')
    figure(1);
       legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns',2)
    h = figure(1);
    h.Position = [7.0729, 5.8021, 1.15*5.8333, 1.15*4.3750];
    Save_as_PDF(h, 'Figures/Comparison_Ballistic_U', 'horizontal')

    
% Accleration   
h = figure(2);
    hold on
    for y = 1:2:3
        for b = 1:2:3
            plot(analytics(y).ballistic(b).n, analytics(y).ballistic(b).z/1000,...
            line_style{b}, 'Color', colors(y,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str_ba));
            plot(-numerics{y}(b,e).a/9.81,numerics{y}(b,e).z/1000,...
                 line_style{b}, 'Color', colors(y+1,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str_bn))
        end
    end
    axis([-0.1, -min([numeric(:,:).ax])/9.81, 0, max([numerics{y}(:,:).zx])/1000])
    grid on; box on;
    legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 2)
    xlabel('$n$ [g]', 'Interpreter', 'Latex')
    ylabel({'$z$';'[km]'}, 'Interpreter', 'Latex')
    
    Save_as_PDF(h, 'Figures/Comparison_Ballistic_a', 'vertical',0,10)
    figure(2);
       legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns',2)
    h = figure(2);
    h.Position = [7.0729, 5.8021, 1.15*5.8333, 1.15*4.3750];
    Save_as_PDF(h, 'Figures/Comparison_Ballistic_a', 'horizontal')
    

% GLIDE

close all

leg_str_ga = "['An. $E =' num2str(E(e)) '$, $\beta =' num2str(beta(b)) '$']";
leg_str_gn = "['Nu. $E =' num2str(E(e)) '$, $\beta =' num2str(beta(b)) '$']";

% Speed
h = figure(1);
    hold on
y = 1;
    for e = 2:3
        for b = 1:2:3
            plot(analytics(y).glide(b,e).u, analytics(y).glide(b,e).z/1000,...
            line_style{b}, 'Color', colors(2*e-3,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str_ga));
            plot(numerics{y}(b,e).u/numerics{y}(b,e).u0,...
                 numerics{y}(b,e).z/1000, line_style{b},...
                 'Color', colors(2*e-2,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str_gn))
        end
    end
    axis([0, max([numerics{y}(:,:).ux])/numerics{y}(b,e).u0, 0, max([numerics{y}(:,:).zx])/1000])
    grid on; box on;
    legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 2)
    xlabel('$U/U_e$', 'Interpreter', 'Latex')
    ylabel({'$z$';'[km]'}, 'Interpreter', 'Latex')
    
    Save_as_PDF(h, 'Figures/Comparison_Glide_U', 'vertical')
    figure(1);
       legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns',2)
    h = figure(1);
    h.Position = [7.0729, 5.8021, 1.15*5.8333, 1.15*4.3750];
    Save_as_PDF(h, 'Figures/Comparison_Glide_U', 'horizontal')
  
% Accleration    
h = figure(2);
    hold on
    for e = 2:3
        for b = 1:2:3
            plot(analytics(y).glide(b,e).n, analytics(y).glide(b,e).z/1000,...
            line_style{b}, 'Color', colors(2*e-3,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str_ga));
            plot(-numerics{y}(b,e).a/9.81,numerics{y}(b,e).z/1000,...
                 line_style{b}, 'Color', colors(2*e-2,:), 'LineWidth', 1.5,...
                 'DisplayName', eval(leg_str_gn))
        end
    end
    axis([-0.1, 6, 0, max([numerics{y}(:,:).zx])/1000])
    grid on; box on;
    legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns', 2)
    xlabel('$n$ [g]', 'Interpreter', 'Latex')
    ylabel({'$z$';'[km]'}, 'Interpreter', 'Latex')
    
    Save_as_PDF(h, 'Figures/Comparison_Glide_a', 'vertical',0,10)
    figure(2);
       legend('Location', 'NorthOutside', 'Interpreter', 'Latex', 'NumColumns',2)
    h = figure(2);
    h.Position = [7.0729, 5.8021, 1.15*5.8333, 1.15*4.3750];
    Save_as_PDF(h, 'Figures/Comparison_Glide_a', 'horizontal')

    
%% LIFT

L = {};
D = {};
h = figure();
    hold on
for e = 1:length(E)
    [rho, ~, ~] = Function_rho(numerics{3}(1,e).z);
    L = cat(1, L, rho.*numerics{3}(1,e).u.^2*E(e)/(2*beta(1)));
    D = cat(1, D, rho.*numerics{3}(1,e).u.^2/(2*beta(1)));
    plot(numerics{3}(1,e).t/60, L{e}(:), 'Color', colors(e,:),...
         'LineWidth', 1.5, 'DisplayName', ['$L$:',' $E = ' num2str(E(e)) ' $'])
    plot(numerics{3}(1,e).t/60, D{e}(:), line_style{3}, 'Color', colors(e,:),...
         'LineWidth', 1.5, 'DisplayName', ['$D$:',' $E = ' num2str(E(e)) ' $'])
end
    axis([0, max([numeric(:,:).tx])/60, 0, max(D{1}(:))])
    Save_as_PDF(h, 'Figures/Lift', 'horizontal')
    h.Position = [7.0729, 5.8021, 1.15*5.8333, 1.15*4.3750];
    grid on; box on;
    xlabel('$t$ [min]', 'Interpreter', 'Latex')
    ylabel({'[N]'}, 'Interpreter', 'Latex')
    legend('Interpreter', 'Latex', 'Location', 'NorthEast', 'NumColumns',3)
    Save_as_PDF(h, 'Figures/Lift', 'horizontal')
    
%% FUGOIDE

T = 2*pi/( sqrt( (1 - .825^2)/( (7.524e3) ) ) )/ 60