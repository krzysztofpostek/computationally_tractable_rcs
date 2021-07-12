% Histograms plotting


figure(2);
subplot(1,2,1);
hold on;

[f,xi] = ksdensity(-EVaR_data(:,1));
plot(xi,f/length(EVaR_data(:,1)),'k','LineWidth',1);
[f,xi] = ksdensity(-EVaR_data(:,2));
plot(xi,f/length(EVaR_data(:,1)),'--k','LineWidth',1);
xlabel('Worst-case EVaR','FontSize',10);
ylabel('Density','FontSize',10);

legend('Robust','Nominal','Location','North');

hold off;

subplot(1,2,2);
hold on;

[f,xi] = ksdensity(return_data(:,1));
plot(xi,f/length(EVaR_data(:,1)),'k','LineWidth',1);
[f,xi] = ksdensity(return_data(:,2));
plot(xi,f/length(EVaR_data(:,1)),'--k','LineWidth',1);
xlabel('Worst-case mean return','FontSize',10);
ylabel('Density','FontSize',10);

hold off;