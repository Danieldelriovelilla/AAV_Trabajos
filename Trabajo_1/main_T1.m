clc
clear all
close all


%% DATOS

M1 = 6.9;
Mm = 29;
gamma = 1.4;
P = 101325;         % Pa
T = 80;             % K 

gas = gas_obj();



%% EJERCICIO 1

resultados_1 = struct();

% Apartado a)
gas.gas_TP(Mm, gamma, P, T);
resultados_1.a = [ gas.rho; gas.e ];

% Apartado b)
[T0, P0, rho0, h0] = gas.remanso(M1, gas.T, gas.P, gas.rho);
resultados_1.b = [ T0; P0 ];

% Apartado c)
resultados_1.c = [ gas.h; h0 ];

% Resultados
resultados_1 = struct2table(resultados_1);


%% EJERCICIO 2

M1 = 8;
resultados_2 = struct();

M2 = sqrt( ( 1 + ((gamma-1)/2)*M1^2 )/...
    ( gamma*M1^2 - (gamma-1)/2 ) );

% Apartado a)
[T01, P01, rho01, h01] = gas.remanso(M1, gas.T, gas.P, gas.rho);

[M2, T2, P2, rho2, h2] = gas.NOCH(M1, gas.T, gas.P, gas.rho, gas.h);
[T02, P02, rho02, h02] = gas.remanso(M2, T2, P2, rho2);

resultados_2.a = T2/T01;

% Apartado b)
resultados_2.b = T02;

% Apartado c)
resultados_2.c = P2/(gas.rho*(gas.a*M1)^2);

% Resultados
resultados_2 = struct2table(resultados_2);


%% EJERCICIO 3

resultados_3 = struct;
M1 = 8;
delta = deg2rad(20);

syms beta
eqn = tan(delta) == 2*cot(beta)*( ( M1^2*(sin(beta))^2 - 1 )/...
                    ( M1^2*(gas.gamma+cos(2*beta)) + 2 ) );
                
beta1 = double(vpasolve(eqn,beta,[0., pi/2]));
beta2 = double(vpasolve(eqn,beta,[pi/4, pi/2]));

% beta debil
Mn1 = M1*sin(beta1);

T2 =  gas.T*( ( 1 + 2*gas.gamma/(gas.gamma+1)*(Mn1^2-1) )*...
        ( (2 + (gas.gamma-1)*Mn1^2)/((gas.gamma+1)*Mn1^2) ) );

resultados_3.a = [rad2deg(beta1); Mn1; T2];

% beta fuerte
Mn1 = M1*sin(beta2);

T2 =  gas.T*( ( 1 + 2*gas.gamma/(gas.gamma+1)*(Mn1^2-1) )*...
        ( (2 + (gas.gamma-1)*Mn1^2)/((gas.gamma+1)*Mn1^2) ) );

resultados_3.a2 = [rad2deg(beta2); Mn1; T2];

% Resultados
resultados_3 = struct2table(resultados_3);


%% EJERCICIO 4

resultados_4 = struct();

L = 2.2;

z = [90, 60, 20]*1e3;
v = [8, 6, 0.5]*1e3;

% for i = 1:3
%     [T(i), P(i), rho(i), a(i), mu(i), l(i)] = ISA(z(i), gas.R, gas.gamma);
% end

T = [ 186.9 , 245.5 , 216.7 ];
a = 340.294*[ 0 , 9.229e-1 , 8.671e-1 ];
P = 101325*[ 1.812e-6 , 2.005e-4 , 5.403e-2 ];
rho = 1.225*[ 2.789e-6 , 2.354e-4 , 7.187e-2 ];
mu = 1.7894e-5*[ 0 , 8.805e-1, 7.945e-1 ];


% 90 km
a(1) = sqrt(gas.gamma*T(1)*8314/28.73);
S1 = 120;
T0 = 291.15;
mu0 = 18.27e-6;
mu(1) = 1.458e-6*T(1)^(3/2)/(T(1)+110.4);
l = mu./P.*sqrt(pi*gas.R.*T/2);

% Resultados: alturas
resultados_4.z = z';
resultados_4.v = v';

% Apartado a)
Re = rho.*v*L./mu;
resultados_4.a = Re';

% Apartado b)
M = v./a;
resultados_4.b = M';

% Apartado c)
Ku = l/L;
resultados_4.c = Ku';

% Resultados
resultados_4 = struct2table(resultados_4);
