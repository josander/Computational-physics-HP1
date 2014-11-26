%% Plot energy, temperature and pressure

% load the data file
clf
clc
data = importdata('energy.data');

set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
Size = size(data);

figure(1);
clf
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

%% Plot the means of the temperature

figure(12);
clf

data = importdata('energy.data');
Size = size(data);
startCut = 0;


for i = startCut+250:500:Size(1)-250
   subplot(2,1,1)
   plot(i, mean(data(i-250+1:i+250,end-1)), 'o')
   hold on 
   
   subplot(2,1,2)
   plot(i, mean(data(i-250+1:i+250,end)), 'o')
   hold on 
end

%% Plot the correlation data
corrSamp = 100;  % The maximum value for k to be plotted, [0.01ps]
corrData = importdata('correlation.data');

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

figure(6);
clf
dispData = importdata('displacement.data');

set(gcf,'renderer','painters','PaperPosition',[0 0 6 8]);
Size = size(dispData);
start_cut = 8000;
steps = Size(1);

plot3(dispData(start_cut:steps,1), dispData(start_cut:steps,2), dispData(start_cut:steps,3),'g  -', 'LineWidth', 0.001);
t = title('Displacement of 1 atom at 500 C$^\circ$','fontsize',14);
set(t,'interpreter','latex');
ylabel('Y [\AA]','Interpreter','latex');
xlabel('X [\AA]','Interpreter','latex');
zlabel('Z [\AA]','Interpreter','latex');
axis equal
grid on

print(gcf,'-depsc2','diffusionLiquid500.eps')

%% Plot the cell size as a function of time
figure(7);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
cellData = importdata('cellSize.data');
plot(cellData(2:end,1), cellData(2:end,2));

title('Length of supercell during equilibration','interpreter','latex','fontsize',14);
xlabel('Time [ps]','Interpreter','latex','fontsize',10);
ylabel('Cell Size [\AA]','Interpreter','latex','fontsize',10);

print(gcf,'-depsc2','cellSize.eps');
%% Import MSD-data for the solid 
sMSDdata = importdata('MSD.data');

% Import MSD-data for the liquid 
lMSDdata = importdata('MSD.data');

% Plot the MSD

figure(8);
clf

set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

plot(sMSDdata(:,1),sMSDdata(:,2), lMSDdata(:,1),lMSDdata(:,2));
title('Mean Square Displacement (MSD)','interpreter','latex','fontsize',14);
y = ylabel('$\Delta_{MSD} ($k$) [$\AA$^2]$','interpreter','latex','fontsize',10);
x = xlabel('Time lag k [ps]','interpreter','latex','fontsize',10);
l = legend('$\Delta_{MSD}$, T = 500C$^\circ$','$\Delta_{MSD}$, T = 900C$^\circ$')
%axis([0 1 0 6])
set(l,'Interpreter','latex');

plotTickLatex2D
print(gcf,'-depsc2','MSD1.eps');

%% Plot convergence of MSD

MSDdata = importdata('MSD.data');
Size = size(MSDdata);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

figure(9);
clf

for(i = 1:Size(1))
   plot(MSDdata(i,1),MSDdata(i,2)/(6*MSDdata(i,1)),'-');
   hold on
end

title('Self Diffusion Coefficient from MSD','interpreter','latex','fontsize',14);
y = ylabel('$\Delta_{MSD}/6k [$\AA$^2/ps]$','interpreter','latex','fontsize',10);
x = xlabel('Time lag k [ps]','interpreter','latex','fontsize',10);
plotTickLatex2D
print(gcf,'-depsc2','selfDiffusion.eps');


%% MSD for 500C 
plot(sMSDdata(:,1),sMSDdata(:,2), lMSDdata(:,1),lMSDdata(:,2));
title('Mean Square Displacement(MSD)','interpreter','latex','fontsize',14);
y = ylabel('$\Delta_{MSD} ($k$) [$\AA$^2]$','interpreter','latex','fontsize',10);
x = xlabel('Time lag k [ps]','interpreter','latex','fontsize',10);
l = legend('$\Delta_{MSD}$, T = 500C$^\circ$','$\Delta_{MSD}$, T = 900C$^\circ$')


axis([0 1 0 0.5]);
set(l,'Interpreter','latex');

plotTickLatex2D
print(gcf,'-depsc2','MSD2.eps');

%% Plot the Velocity correlation function

figure(10);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
plot(sMSDdata(:,1),sMSDdata(:,3)/sMSDdata(1,3), lMSDdata(:,1),lMSDdata(:,3)/lMSDdata(1,3));
title('Velocity correlation function','interpreter','latex','fontsize',14);

ylabel(' $\Phi_v$(k) [-]','interpreter','latex','fontsize',10);
xlabel('Time lag $k$ [ps]','interpreter','latex','fontsize',10);
l = legend('$\Phi_v$, T = 500C$^\circ$','$\Phi_v$, T = 900C$^\circ$')
axis([0 1 -0.5 1])
set(l,'Interpreter','latex');
plotTickLatex2D
print(gcf,'-depsc2','VCF.eps');

%% Plot the Spectral function

omfact=10; % to convert to ps^-1
figure(8);

clf

set(gcf,'renderer','painters','PaperPosition',[0 0 5 4]);
plot(omfact*sMSDdata(:,4),sMSDdata(:,5),omfact*lMSDdata(:,4),lMSDdata(:,5));

title('Spectral function','interpreter','latex','fontsize',14);

ylabel(' $\hat{\Phi}_v(\omega$) [ps/\AA]','interpreter','latex','fontsize',10);
xlabel(' $\omega$ [ps$^{-1}$]','interpreter','latex','fontsize',10);
l = legend('$\Phi_v$, T = 500C$^\circ$','$\Phi_v$, T = 900C$^\circ$');

set(l,'Interpreter','latex');
plotTickLatex2D
print(gcf,'-depsc2','spectral.eps');
