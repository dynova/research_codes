Stot = Fig4aScript;
Qvec = logspace(-2,4,100);
lw = 3;
hold all

plot(Qvec,Stot,'k','Linewidth',lw);
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XMinorTick','on',...
    'YMinorTick','on','XScale','log','YScale','log');
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
xlim([1e-2 1e4]);
ylim([1e1 1e9]);
xlabel('Synthesis rate of Substrate (Q) [nM s^{-1}]','FontWeight','bold','FontSize',20);
ylabel('Total Substrate [S]_T [nM]','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.08, 0.5, 0]);
box off
hold off
export_fig Fig4a -pdf -painters -transparent