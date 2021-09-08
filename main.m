clc; clear all; close all

%% 1D elastic model parameters
h=[10 10 10]*1e-3;
vp=[200 400 600 1000]*1e-3;
vs=[100 200 300 500]*1e-3;
rho=[1.9 1.9 1.9 1.9];

%% Frequency and phase velocity axis
f=1:3;
c=(80:0.25e-2:500)*1e-3;

%% Computation of Rayleigh waves modes
cR = dispersion_modes(vs,vp,rho,h,f,c);

%% Figure
figure
for i=1:size(cR,1)
plot(f,cR,'Linewidth',2), hold on, grid on
end
set(gca,'fontsize',18,'TickLabelInterpreter','latex'), 
title('Rayleigh dispersion curves','FontSize',22,'Interpreter','Latex')
ylabel('Rayleigh phase velocity (m/s)','FontSize',22,'Interpreter','Latex'), 
xlabel('Frequency (Hz)','FontSize',22,'Interpreter','Latex')
%export_fig('Rayleigh_dispersion.pdf','-transparent') 