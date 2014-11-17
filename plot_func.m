% plot

% load the data file
clf
data = importdata('energy.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);

%%


%plot

plot(data(:,1),data(:,2:end-2));


% labels
ylabel('Energy');
xlabel('Time');
print(gcf,'-depsc2','energy.eps')
%% 
plot(data(:,1),data(:,end-1));

ylabel('Temerature');
xlabel('Time');
print(gcf,'-depsc2','temperature.eps')

%% 
plot(data(:,1),data(:,end), [0 100],[6.32*10^-7 6.32*10^-7]);

ylabel('Pressure');
xlabel('Time');
print(gcf,'-depsc2','pressure.eps')
