load Full.mat;
load PDS.mat;
lw = 3;
Q2vec = logspace(-2,4,100);
figure(1)
hold all
g1 = plot(Q2vec,A,'Color','k','LineWidth',lw);
g2 = plot(Q2vec,B,'Color','k','LineStyle',':','LineWidth',lw);
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XMinorTick','on',...
    'YMinorTick','on','XScale','log','YScale','log');
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
xlabel('Synthesis Rate of Substrate 2 (Q_2) [nM s^{-1}]','FontWeight','bold','FontSize',20);
ylabel('Total Substrate [S_1]_T [nM]','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.08, 0.5, 0]);
s = {'Full Model', 'Proc. E3, Dist./Seq. DUB Model'};
t = legend(s,'Location','West','FontSize',16,'box','off');
box off
xlim([min(Q2vec) max(Q2vec)]);
ylim([1e3 1e4]);
hold off
export_fig Fig4C_revised -pdf -painters -transparent