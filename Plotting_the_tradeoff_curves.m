%%%% This file plots the tradeoff curves for the nominal and robust
%%%% portfolios
plot(EVaR_range,ResultsReversed)
hold on
plot(EVaR_range,Mysterious_curve(:,2))
plot(Mysterious_curve(:,1),ResultsReversed(:,2),'--')
plot(Mysterious_curve(:,1),Mysterious_curve(:,2),'--')
legend('Robust portfolio','Nominal portfolio','Nominal portfolio - nominal EVaR vs worst case return','Nominal portfolio - worst case EVaR vs nominal return','Nominal portfolio - worst case EVaR vs worst case return')
%%%% End of file