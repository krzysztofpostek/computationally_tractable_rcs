figure(1)

plot(EVaR_range,Optimal_mu_for_EVaR(:,1),'LineWidth',2);

hold on;

scatter(Risky_assets_wc_EVaR_corresponding_return(:,1),Risky_assets_wc_return_corresponding_EVaR(:,2),'k','filled');
legend('Efficient frontier','Individual risky assets','Location','NorthWest');
ylabel('Worst-case mean monthly return');
xlabel('Worst-case monthly EVaR');

hold off;