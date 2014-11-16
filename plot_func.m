% plot

% load the data file
clf
data = importdata('energy.data');
%%

set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);

%plot
hold on
plot(data(:,1),data(:,2:end));


% labels
xlabel('Energy');
ylabel('Time');
print(gcf,'-depsc2','samp2.eps')

