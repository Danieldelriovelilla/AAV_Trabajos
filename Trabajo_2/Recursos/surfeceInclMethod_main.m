clc
close all
clear all


filename    = 'CAD_3/stl_3.stl';
[p,t,tnorm] = import_stl_fast(filename,1);

figure(2)
patch('vertices',p,'faces',t,'FaceColor','none','edgecolor','k')

alpha = deg2rad(60);
U_infty = zeros(size(tnorm));
U_infty(:,1) = 0;
U_infty(:,2) = cos(alpha);
U_infty(:,3) = sin(alpha);

% U_infty(:,2) = cos(alpha);
% U_infty(:,1) = sin(alpha);

Stheta = zeros(length(tnorm),1);

for I=1:length(tnorm)

    Stheta(I,1) = dot(U_infty(I,:),tnorm(I,:))./norm(U_infty(I,:));
    
end

% Cálculo del cp. Método de Newton
cps = 2*Stheta.^2;
cps = [cps;0];

% Cálculo del cp en las zonas de sotavento
I_cpO           = find(Stheta>0);
cps(I_cpO,1)    = 0;


figure(1)
colormap(jet)
set(gca,'CLim',[0 2])
pat = patch('vertices',p,'faces',t,'FaceColor','b','edgecolor','none');
xlabel('x'),ylabel('y'),zlabel('z'),


set(gca,'CLim',[0 2])
cdata = cps;
set(pat,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled');
colorbar
title('Numetico \alpha = 45')
view([15,-15])
axis equal

