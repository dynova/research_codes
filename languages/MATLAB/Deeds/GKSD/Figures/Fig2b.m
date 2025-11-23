load 2b
lw = 3;
hold all
h1 = semilogx(rvec,KMI1./KMI2,'r:','LineWidth',lw);
h2 = semilogx(rvec,KMI3./KMI4,'r','LineWidth',lw);
h3 = semilogx(rvec,KMGK1./KMGK2,'b:','LineWidth',lw);
h4 = semilogx(rvec,KMGK3./KMGK4,'b','LineWidth',lw);
h5 = semilogx(rvec,KMSD1./KMSD2,'k:','LineWidth',lw);
h6 = semilogx(rvec,KMSD3./KMSD4,'k','LineWidth',lw);
h7 = semilogx(NaN, NaN,'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',lw);
h8 = semilogx(NaN, NaN,'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',lw);
h9 = semilogx(NaN, NaN,'Color',[1 1 1],'LineStyle',':','LineWidth',lw);
h10 = semilogx(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
hold off
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XMinorTick','on','YTick',...
    0.1:0.1:1,'XScale','log');
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
xlabel('r','FontWeight','bold','FontSize',20);
ylabel('Normalized Fraction S^*','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.09, 0.5, 0]);
% title('Varying K_M as the measure of saturation','FontWeight','bold','FontSize',22)
s={' GK',' Intermediate',' Full','','',' Unsaturated', ' Saturated'};
f = legend([h4 h2 h6 h9 h10 h7 h8],s);
set(f,'Location','SouthEast',...
    'FontSize',16,'box','off');
box off
export_fig Fig2b -pdf -painters -transparent