load DDS1.mat
load DDS2.mat
rvec = logspace(-2,4,100);
figure(1)
lw = 3;
hold all
g1 = plot(rvec,DDS11./DDS12,'Color','k','LineStyle',':','LineWidth',lw);
g2 = plot(rvec,DDS21./DDS22,'Color','k','LineStyle','-','LineWidth',lw);
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XMinorTick','on','YTick',...
    0.1:0.1:1,'XScale','log');
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
xlabel('r','FontWeight','bold','FontSize',20);
ylabel('Normalized Fraction S^*','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.08, 0.5, 0]);
s={' Unsaturated', ' Saturated'};
f = legend([g1 g2],s);
set(f,'Location','East',...
    'FontSize',16,'box','off');
box off
hold off
export_fig DDSa -pdf -painters -transparent

figure(2)
lw = 3;
hold all
k1 = plot(rvec,DDS12,'Color','k','LineStyle',':','LineWidth',lw);
k2 = plot(rvec,DDS22,'Color','k','LineStyle','-','LineWidth',lw);
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XMinorTick','on',...
    'YMinorTick','on','XScale','log','YScale','log','YTick',...
    [1e1, 1e2, 1e3, 1e4]);
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
xlabel('r','FontWeight','bold','FontSize',20);
ylabel('Total Substrate [S]_T [nM]','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.08, 0.5, 0]);
s={' Unsaturated', ' Saturated'};
f = legend([k1 k2],s);
set(f,'Location','East',...
    'FontSize',16,'box','off');
box off
hold off
export_fig DDSb -pdf -painters -transparent