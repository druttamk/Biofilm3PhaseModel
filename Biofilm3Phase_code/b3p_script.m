%% 3 phase cell ecm and water with osmotic influx
%  in cylindrical coordinates
clc; clear; close all;

code = "Biofilm3Phase_code.m";

T = 100; dr = 0.05; R = 12;  r = 0:dr:R; n = length(r);
% Dimensionless parameters and their numerical value
psid = 0; gmma = 1; k1 = 0.5;
D = 0.9917;  kmb = 1.21000; kms = 0.30250; pe = 4.033; %pe = 4.03;%2.52
psin = 0.3; psim = 0.7; kb = 6; Q = 0.24;
phic0 = psin - 0.05; phie0 = psim - 0.05; phia0 = 0.1;

% Initial conditions
Rin = 1; b = 5*10^-3; H0 = 0.2; hstar = 2*b; %Rin initial radius of the biofilm
h = (H0-b).*(1 - tanh(5.*(r - Rin)))./2 + b; 
phic = phic0.*(1 - tanh(5.*(r - Rin)))./2; 
phie = phie0.*(1 - tanh(5.*(r - Rin)))./2;
phia = phia0.*(1 - tanh(5.*(r - Rin)))./2; 

cs = ones(1,n); cb = zeros(1,n);
y0 = [h,phic,phie,phia,cs,cb];

tf = [0 T];
solq = ode15s(@(t,c) thinfm(t,c,dr,n,psin,psim,Q,gmma,kb,D, kmb, kms, pe, hstar), tf, y0);
 
%% ploting the solution profiles
dt = 1; t = 0:dt:T; ssol = deval(solq,t);
nt = length(t);

for i = 1:15*dt:nt
% 
hv = @(h) heaviside(h - hstar);

h = ssol(1:n,i); phic = ssol(n+1:2*n,i);
phie = ssol(2*n+1:3*n,i); phia = ssol(3*n+1:4*n,i);
cs = ssol(4*n+1:5*n,i); cb = ssol(5*n+1:6*n,i);

phit = phic + phie + phia; v0 = Q.*(phie.^3 - 0.04^3);

gc = psin.*phic.*cb./(k1 + cb)  - psid.*phic; ge = psim.*phic.*cb./(k1 + cb)  + psid.*phic;
gh = (gc + ge).*hv(h); G = gh./max(gh);

figure(1), %clf
subplot(2, 2, 1),plot(r,h), hold on,xlabel('Radial distance (r)'), ylabel('Height of biofilm (h)')
subplot(2, 2, 2),plot(r,phic), hold on,xlabel('Radial distance (r)'), ylabel('\phi_c')
subplot(2, 2, 3),plot(r,phie), hold on,xlabel('Radial distance (r)'), ylabel('\phi_e')
subplot(2, 2, 4),plot(r,phia), hold on,xlabel('Radial distance (r)'), ylabel('\phi_a')

figure(2), %clf
subplot(2, 2, 1),plot(r,cs), hold on, xlabel('Radial distance (r)'), ylabel('Nutrient conc in sub. n_s')
subplot(2, 2, 2),plot(r,cb), hold on, xlabel('Radial distance (r)'), ylabel('Nutrient conc. in bio. n_b')
subplot(2, 2, 3),plot(r,v0), hold on,xlabel('Radial distance (r)'), ylabel('Volumetric flow (v_0)')
subplot(2, 2, 4),plot(r,G), hold on,xlabel('Time (t)'), ylabel('Growth rate (G/max(G))')

end













 
