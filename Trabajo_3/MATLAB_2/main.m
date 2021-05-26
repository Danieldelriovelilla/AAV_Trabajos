clc
clear
close all


%% PARAMETROS DE ESTUDIO

E = [0 0.5 1];
beta = [50 280 500];
ze = [120e3, 120e3, 120e3];
ue = [7e3, 9e3, 7e3];
gammae = [-5*pi/180, -5*pi/180, -7.5*pi/180];
ci = [ze; ue; gammae];

% Paso de tiempo del integrador
Dt = 0.1;

 

%% DEFINIR PLOT
colors = [0, 64, 255;
          0, 191, 255;
          191, 0, 255;
          144, 12, 63]/255;

line_style = {'-','--',':'};


%% ESTUDIO CON VARIAS CONDICIONES INICIALES

for c = 1:length(ci)
    
    
    % ANALITICO
    [ub{c}, nb{c}, us{c}, ns{c}, z] = f_Analitico(ue(c), gammae(c), ze(c), beta, E);
    
    % NUMERICO
    limsx = [];
    limsm = [];
    for i = 1:length(beta)
        for j = 1:length(E)
            sol_ = RK4(beta(i), E(j), ci(:,c), Dt);
            limsx = cat(1,limsx,max(sol_));
            limsm = cat(1,limsm,min(sol_));
            sol{i,j} = sol_; 
        end
    end
    solucion{c} = sol;
    
    % PLOTS
    limx(c,:) = max(limsx);
    limm(c,:) = min(limsm);
    for i = 1:length(beta)
        for j = 1:length(E)
            h{1} = figure(1);
                hold on
                plot(solucion{c}{i,j}(:,1), rad2deg(solucion{c}{i,j}(:,end)),...
                     line_style{i}, 'Color', colors(j,:), 'LineWidth', 1.5,...
                     'DisplayName', ['$\beta = ' num2str(beta(i)) ' $ ; $E = ' num2str(E(j)) '$'])
                axis([0, limx(c,1), -90, rad2deg(limx(c,end))])
                grid on; box on;
                h{1}.Position = [680   558   1.4*560   420];
                legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)
                xlabel('$t$ [s]', 'Interpreter', 'Latex')
                ylabel({'$\gamma$';'[$\mathrm{^o}$]'}, 'Interpreter', 'Latex')

            h{2} = figure(2);
                hold on
                plot(solucion{c}{i,j}(:,1), solucion{c}{i,j}(:,3)/ue(c),...
                    line_style{i}, 'Color', colors(j,:), 'LineWidth', 1.5,...
                    'DisplayName', ['$\beta = ' num2str(beta(i)) ' $ ; $E = ' num2str(E(j)) '$'])
                axis([0, limx(c,1), 0, limx(c,3)/ue(c)])
                grid on; box on;
                h{2}.Position = [680   558   1.4*560   420];
                legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)
                xlabel('$t$ [s]', 'Interpreter', 'Latex')
                ylabel({'$\frac{U}{U_0}$'}, 'Interpreter', 'Latex')

            h{3} = figure(3);
                hold on
                plot(solucion{c}{i,j}(:,1),solucion{c}{i,j}(:,2)/1000,...
                    line_style{i}, 'Color', colors(j,:), 'LineWidth', 1.5,...
                    'DisplayName', ['$\beta = ' num2str(beta(i)) ' $ ; $E = ' num2str(E(j)) '$'])
                axis([0, limx(c,1), 0, limx(c,2)/1000])
                grid on; box on;
                h{3}.Position = [680   558   1.4*560   420];
                legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)
                xlabel('$t$ [s]', 'Interpreter', 'Latex')
                ylabel({'$z$';'[km]'}, 'Interpreter', 'Latex')

            h{4} = figure(4);
                hold on
                plot(solucion{c}{i,j}(:,1),-solucion{c}{i,j}(:,4)/9.81,...
                    line_style{i}, 'Color', colors(j,:), 'LineWidth', 1.5,...
                    'DisplayName', ['$\beta = ' num2str(beta(i)) ' $ ; $E = ' num2str(E(j)) '$'])
                axis([0, limx(c,1), 0, -limm(c,4)/9.81])
                grid on; box on;
                h{4}.Position = [680   558   1.4*560   420];
                legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)
                xlabel('$t$ [s]', 'Interpreter', 'Latex')
                ylabel({'$a$';'[g]'}, 'Interpreter', 'Latex')


        end
    end
    
    % Gurardar figuras en PDF
    figure(1);
    Save_as_PDF(h{1}, ['Figures\G_ci_' num2str(c)], 'horizontal')

    figure(2);
    Save_as_PDF(h{2}, ['Figures\U_ci_' num2str(c)], 'horizontal')

    figure(3);
    Save_as_PDF(h{3},  ['Figures\Z_ci_' num2str(c)], 'horizontal')

    figure(4);
    Save_as_PDF(h{4},  ['Figures\A_ci_' num2str(c)], 'horizontal')

    close all
    
end




%% VALIDACION RESULTADOS ANALITICOS

close all

for c = 1:2:length(ci)
    
    for i = 1:2:length(beta)

        h{1} = figure(1);
            hold on
            plot(ub{c}{i}, z/1000,line_style{i}, 'Color', colors(c,:),...
                 'LineWidth', 1.5, 'DisplayName',...
                 ['$A - \gamma_e = ' num2str(rad2deg(gammae(c))) ' $ ; $\beta = ' num2str(beta(i)) '$'])
            plot(solucion{c}{i,1}(:,3)/ue(c), solucion{c}{i,1}(:,2)/1000,line_style{i}, 'Color', colors(c+1,:),...
                 'LineWidth', 1.5, 'DisplayName',...
                 ['$N - \gamma_e = ' num2str(rad2deg(gammae(c))) ' $ ; $\beta = ' num2str(beta(i)) '$'])
             axis([0, limx(1,3)/ue(1), 0, limx(1,2)/1000])
            grid on; box on;
            h{1}.Position = [680   558   1.4*560   420];
            legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)    


        h{2} = figure(2);
            hold on
            plot(nb{c}{i}/9.81, z/1000,line_style{i}, 'Color', colors(c,:),...
                 'LineWidth', 1.5, 'DisplayName',...
                 ['$A - \gamma_e = ' num2str(rad2deg(gammae(c))) ' $ ; $\beta = ' num2str(beta(i)) '$'])
            plot(-solucion{c}{i,1}(:,4)/9.81, solucion{c}{i,1}(:,2)/1000,line_style{i}, 'Color', colors(c+1,:),...
                 'LineWidth', 1.5, 'DisplayName',...
                 ['$N - \gamma_e = ' num2str(rad2deg(gammae(c))) ' $ ; $\beta = ' num2str(beta(i)) '$'])
             axis([-0.1, 19.8, 0, limx(1,2)/1000])
            grid on; box on;
            h{2}.Position = [680   558   1.4*560   420];
            legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)     
    end
    
end

for i = 1:length(beta)
    for j = 2:length(E)
        h{3} = figure(3);
            hold on
            plot(us{1}{i,j}, z/1000,line_style{j}, 'Color', colors(i,:),...
                 'LineWidth', 1.5, 'DisplayName',...
                 ['$A - E = ' num2str(E(j)) ' $ ; $\beta = ' num2str(beta(i)) '$'])
            plot(solucion{1}{i,j}(:,3)/ue(c), solucion{1}{i,j}(:,2)/1000,line_style{i}, 'Color', colors(i+1,:),...
                 'LineWidth', 1.5, 'DisplayName',...
                 ['$N - E = ' num2str(E(j)) ' $ ; $\beta = ' num2str(beta(i)) '$'])
            axis([0, limx(1,3)/ue(1), 0, limx(1,2)/1000])
            grid on; box on;
            h{3}.Position = [680   558   1.4*560   420];
            legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)    

        h{4} = figure(4);
            hold on
            plot(ns{1}{i,j}/9.81, z/1000,line_style{j}, 'Color', colors(i,:),...
                 'LineWidth', 1.5, 'DisplayName',...
                 ['$A - E = ' num2str(E(j)) ' $ ; $\beta = ' num2str(beta(i)) '$'])
            plot(-solucion{1}{i,j}(:,4)/9.81, solucion{1}{i,j}(:,2)/1000,line_style{i}, 'Color', colors(i+1,:),...
                 'LineWidth', 1.5, 'DisplayName',...
                 ['$N - E = ' num2str(E(j)) ' $ ; $\beta = ' num2str(beta(i)) '$'])
            axis([-0.1, 6.2, 0, limx(1,2)/1000])
            grid on; box on;
            h{4}.Position = [680   558   1.4*560   420];
            legend('Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns', 1)    
    end
end