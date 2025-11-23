load DDS3.mat
load DDS4.mat
DDS3 = load('~/Desktop/Programming/cpp/Deeds/SingleSubstrate/DDS_Unsat.txt');
DDS4 = load('~/Desktop/Programming/cpp/Deeds/SingleSubstrate/DDS_Sat.txt');
rvec = logspace(-2,1,100);
cvec = logspace(-2,1,4);
figure(1)
lw = 3;
hold all
g1 = plot(rvec,DDS31./DDS32,'Color','k','LineStyle',':','LineWidth',lw);
g2 = plot(rvec,DDS41./DDS42,'Color','k','LineStyle','-','LineWidth',lw);
g3 = plot(cvec,DDS3(:,1),'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');
g4 = plot(cvec,DDS4(:,1),'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');            
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
set(f,'Location','West',...
    'FontSize',16,'box','off');
box off
hold off
export_fig DDSc -pdf -painters -transparent

figure(2)
lw = 3;
hold all
k1 = plot(rvec,DDS32,'Color','k','LineStyle',':','LineWidth',lw);
k2 = plot(rvec,DDS42,'Color','k','LineStyle','-','LineWidth',lw);
k3 = plot(cvec,DDS3(:,2)./10,'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');
k4 = plot(cvec,DDS4(:,2)./10,'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');            
hold off
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
set(f,'Location','West',...
    'FontSize',16,'box','off');
box off
hold off
export_fig DDSd -pdf -painters -transparent