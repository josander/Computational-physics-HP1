%% Plot energy, temperature and pressure

% load the data file
clf
clear all
clc
%%
data = importdata('energy.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
Size = size(data);
%%
figure(1);
plot(data(:,1),data(:,2:end-2));

% labels
title('Energy','interpreter','latex','fontsize',14);
y = ylabel('Energy [eV]','interpreter','latex','fontsize',10);
xlabel('Time [ps]','interpreter','latex','fontsize',10);
plotTickLatex2D
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
e = legend('Total Energy','Kinetic Energy','Potential energy');
set(e,'Interpreter','latex', 'Location', 'east')
print(gcf,'-depsc2','energy.eps') 

figure(2);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
plot(data(:,1),data(:,end-1));
meanTemp = mean(data(1500:Size(1),end-1))
title('Temperature','interpreter','latex','fontsize',14);
y = ylabel('Temperature [K]','interpreter','latex','fontsize',10);
xlabel('Time [ps]','interpreter','latex','fontsize',10);
plotTickLatex2D
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
print(gcf,'-depsc2','temperature.eps')

figure(3);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
plot(data(:,1),data(:,end));
meanPress = mean(data(1000:Size(1),end))
title('Pressure','interpreter','latex','fontsize',14);
y = ylabel('Pressure [eV/\AA$^3$]','interpreter','latex','fontsize',10);
xlabel('Time [ps]','interpreter','latex','fontsize',10);
%axis([0 Size(1)*0.01 -0.004 0.004]);
plotTickLatex2D
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
print(gcf,'-depsc2','pressure.eps')
%%
sT= 5.028783; 	
sP= 4.884246; 

%% Plot the correlation data
corrSamp = 100;  % The maximum value for k to be plotted, [0.01ps]
figure(4);
corrData = importdata('correlation.data');
%%
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
Size = size(corrData);

figure(4);
plot((0:0.01:(corrSamp-0.1)/100)',corrData(1:corrSamp,1));
hold on
plot([0 (corrSamp-0.01)/100], [exp(-2) exp(-2)],'g-');
%axis([0 100 -0.4 0.4])
title('Correlation function for temperature','interpreter','latex','fontsize',14);
y = ylabel('$\Phi _T (k)$ [-]','interpreter','latex','fontsize',10);
x = xlabel('Time lag $k$ [ps]','interpreter','latex','fontsize',10);
%s4et(x, 'Units', 'Normalized', 'Position', [0.5, 0, 0]);
hold off
plotTickLatex2D
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
l = legend('Correlation function','$ e^{-2}$');
set(l,'Interpreter','latex');
print(gcf,'-depsc2','correlationT.eps')

figure(5);
plot((0:0.01:(corrSamp-0.1)/100)',corrData(1:corrSamp,2));
hold on
plot([0 (corrSamp-0.01)/100], [exp(-2) exp(-2)],'g-');
%axis([0 100 0 0.4])
title('Correlation function for pressure','interpreter','latex','fontsize',14);
y = ylabel('$\Phi _P (k)$ [-]','interpreter','latex','fontsize',10);
x = xlabel('Time lag $k$ [ps]','interpreter','latex','fontsize',10);
%set(x, 'Units', 'Normalized', 'Position', [0.5, -0.01, 0]);
hold off
plotTickLatex2D
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
l = legend('Correlation function','$ e^{-2}$');
set(l,'Interpreter','latex');

print(gcf,'-depsc2','correlationP.eps')
%% Plot the displacement in 3D

figure(5);

dispData = importdata('displacement.data');
%%
set(gcf,'renderer','painters','PaperPosition',[0 0 6 8]);
%Size = size(corrData);

steps = 10000;

plot3(dispData(1:steps,1), dispData(1:steps,2), dispData(1:steps,3),'b  --', 'LineWidth', 0.001);
t = title('Displacement of 1 atom at 500 C$^\circ$','fontsize',14);
set(t,'interpreter','latex');
ylabel('Y [\AA]','Interpreter','latex');
xlabel('X [\AA]','Interpreter','latex');
zlabel('Z [\AA]','Interpreter','latex');
axis equal
grid on

print(gcf,'-depsc2','diffusionLiquid500.eps')

%% Plot the MSD

figure(6);
MSDdata = importdata('MSD.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);

plot(MSDdata(:,1),MSDdata(:,2));
ylabel('Displacement');
xlabel('Time');

print(gcf,'-depsc2','MSD.eps');

%% Plot the Velocity correlation function

figure(7);

MSDdata = importdata('MSD.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);

plot(MSDdata(:,1),MSDdata(:,3));
ylabel('Velocity');
xlabel('Time');

print(gcf,'-depsc2','VCF.eps');

%% Plot the Spectral function

figure(8);
MSDdata = importdata('MSD.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);

plot(MSDdata(:,4),MSDdata(:,5));
ylabel('Dont know');
xlabel('Omega');

print(gcf,'-depsc2','spectral.eps');
