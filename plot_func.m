%% Plot energy, temperature and pressure

% load the data file
clf
clear all
clc
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
y = ylabel('Pressure [$eV/\AA^3$]','interpreter','latex','fontsize',10);
xlabel('Time [ps]','interpreter','latex','fontsize',10);
plotTickLatex2D
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
print(gcf,'-depsc2','pressure.eps')
%%
sT= 5.028783; 	
sP= 4.884246; 

%% Plot the correlation data

figure(4);
corrData = importdata('correlation.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);
Size = size(corrData);

subplot(2,1,1);
plot(corrData(:,1));
hold on
plot([0 Size(1)], [exp(-2) exp(-2)],'g-');

subplot(2,1,2);
plot(corrData(:,2));
hold on
plot([0 Size(1)], [exp(-2) exp(-2)],'g-');

%% Plot the displacement in 3D

figure(5);

dispData = importdata('displacement.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);
Size = size(corrData);

plot3(dispData(:,1), dispData(:,2), dispData(:,3));
ylabel('Y');
xlabel('X');
zlabel('Z');

print(gcf,'-depsc2','diffusionLiquid.eps')

%% Plot the MSD

figure(6);

MSDdata = importdata('MSD.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);

plot(MSDdata(:,1));
ylabel('Displacement');
xlabel('Timestep');

print(gcf,'-depsc2','MSD.eps');

