%% Plot energy, temperature and pressure

lengthEq = 5000;

% load the data file
clc
data = importdata('energy.data');

set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
Size = size(data);

figure(1);
clf
plot(data(:,1),data(:,2:end-2));

figure(2);
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
plot(data(:,1),data(:,end-1));
meanTemp = mean(data(lengthEq:Size(1),end-1))

figure(3);
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
plot(data(:,1),data(:,end));
meanPress = mean(data(lengthEq:Size(1),end))

% Plot the means of the temperature

figure(4);
clf

Size = size(data);
startCut = 1; 

for i = startCut+250:500:Size(1)-250
   subplot(2,1,1)
   plot(i, mean(data(i-250:i+250,end-1)), 'o')
   hold on 

   subplot(2,1,2)
   plot(i, mean(data(i-250:i+250,end)), 'o')
   hold on 
end

% Plot the correlation data
corrSamp = 100;  % The maximum value for k to be plotted, [0.01ps]
corrData = importdata('correlation.data');

set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
Size = size(corrData);

figure(5);
subplot(2,1,1)
plot((0:0.01:(corrSamp-0.1)/100)',corrData(1:corrSamp,1));
hold on
plot([0 (corrSamp-0.01)/100], [exp(-2) exp(-2)],'g-');
hold off

subplot(2,1,2)
plot((0:0.01:(corrSamp-0.1)/100)',corrData(1:corrSamp,2));
hold on
plot([0 (corrSamp-0.01)/100], [exp(-2) exp(-2)],'g-');
hold off

% Plot the displacement in 3D

figure(6);
clf
dispData = importdata('displacement.data');
title('Displacement in 3D','interpreter','latex','fontsize',14);


set(gcf,'renderer','painters','PaperPosition',[0 0 6 8]);
Size = size(dispData);
start_cut = 1;
steps = Size(1);

plot3(dispData(start_cut:steps,1), dispData(start_cut:steps,2), dispData(start_cut:steps,3),'g  -', 'LineWidth', 0.001);
axis equal
grid on

% Plot the cell size as a function of time
figure(7);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
cellData = importdata('cellSize.data');
plot(cellData(2:end,1), cellData(2:end,2));

title('Change in lattice parameter','interpreter','latex','fontsize',14);


% Import MSD-data for the solid 
sMSDdata = importdata('MSD.data');

% Import MSD-data for the liquid 
lMSDdata = importdata('MSD.data');

% Plot the MSD

figure(8);
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
plot(sMSDdata(:,1),sMSDdata(:,2));

% Plot convergence of MSD

MSDdata = importdata('MSD.data');
Size = size(MSDdata);
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

figure(9);
clf
title('Self diffusion coefficient from MSD','interpreter','latex','fontsize',14);

for(i = 1:Size(1))
   plot(MSDdata(i,1),MSDdata(i,2)/(6*MSDdata(i,1)),'o');
   hold on
end

% Plot the Velocity correlation function

figure(10);
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
plot(sMSDdata(:,1),sMSDdata(:,3)/sMSDdata(1,3));
title('Velocity correlation function','interpreter','latex','fontsize',14);

mean(sMSDdata(:,3))


% Plot the Spectral function

omfact=10; % to convert to ps^-1
figure(11);
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 5 4]);
plot(omfact*sMSDdata(:,4),sMSDdata(:,5));

title('Spectral function','interpreter','latex','fontsize',14);

