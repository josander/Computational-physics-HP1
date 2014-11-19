
clf
%data = importdata('energy.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 8]);
%%
data1 = importdata('energy.data');
%% Plot first time-step

subplot(3,1,1)
plot(data1(:,1),data1(:,2:end-2));

% labels
title('Timestep = 0.02 ps','interpreter','latex','fontsize',14);
ylabel('Energy [eV]','interpreter','latex','fontsize',10);
xlabel('Time [ps]','interpreter','latex','fontsize',10);
 	
axis([0 50 0 880]);
breakyaxis([50,830]);
%plotTickLatex2D

%%
data2 = importdata('energy.data');
%%
subplot(3,1,2)
plot(data2(:,1),data2(:,2:end-2));

% labels
title('Timestep = 0.01 ps','interpreter','latex','fontsize',14);
ylabel('Energy [eV]','interpreter','latex','fontsize',10);
xlabel('Time [ps]','interpreter','latex','fontsize',10);
 	
axis([0 50 0 870]);
breakyaxis([25,833]);
%plotTickLatex2D
 %%
data3 = importdata('energy.data');
%%
subplot(3,1,3)
plot(data3(:,1),data3(:,2:end-2));

% labels
title('Timestep = 0.001 ps','interpreter','latex','fontsize',14);
ylabel('Energy [eV]','interpreter','latex','fontsize',10);
xlabel('Time [ps]','interpreter','latex','fontsize',10);
 	
axis([0 50 0 870]);
breakyaxis([25,833]); 
%plotTickLatex2D

print(gcf,'-depsc2','timestep.eps') 