%% Plot energy, temperature and pressure

% load the data file
clf
clear all
clc
data = importdata('energy.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);
Size = size(data);

figure(1);
plot(data(:,1),data(:,2:end-2));

% labels
ylabel('Energy');
xlabel('Time');
print(gcf,'-depsc2','energy.eps') 

figure(2);
plot(data(:,1),data(:,end-1));
meanTemp = mean(data(1500:Size(1),end-1))

ylabel('Temerature');
xlabel('Time');
print(gcf,'-depsc2','temperature.eps')

figure(3);
plot(data(:,1),data(:,end));
meanPress = mean(data(1000:Size(1),end))

ylabel('Pressure');
xlabel('Time');
print(gcf,'-depsc2','pressure.eps')

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
