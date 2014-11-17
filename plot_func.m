% Plot energy, temperature and pressure

% load the data file
clf
data = importdata('energy.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);

%%

figure(1);
plot(data(:,1),data(:,2:end-2));

% labels
ylabel('Energy');
xlabel('Time');
print(gcf,'-depsc2','energy.eps')
%% 

figure(2);
plot(data(:,1),data(:,end-1));

ylabel('Temerature');
xlabel('Time');
print(gcf,'-depsc2','temperature.eps')
%%
plot(data(:,1),data(:,end), [0 100],[6.32*10^-7 6.32*10^-7]);



ylabel('Pressure');
xlabel('Time');
print(gcf,'-depsc2','pressure.eps')

%%

corrData = importdata('correlation.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);

%%
subplot(2,1,1);
plot(corrData(:,1));

subplot(2,1,2);
plot(corrData(:,2));
