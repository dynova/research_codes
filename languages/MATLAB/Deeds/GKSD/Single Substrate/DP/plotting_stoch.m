load DP3.mat
load DP4.mat
DP3 = load('~/Desktop/Programming/cpp/Deeds/SingleSubstrate/DP_Unsat.txt');
DP4 = load('~/Desktop/Programming/cpp/Deeds/SingleSubstrate/DP_Sat.txt');
rvec = logspace(-2,1,100);
cvec = logspace(-2,1,4);
figure(1)
lw = 3;
hold all
g1 = plot(rvec,DP31./DP32,'Color','k','LineStyle',':','LineWidth',lw);
g2 = plot(rvec,DP41./DP42,'Color','k','LineStyle','-','LineWidth',lw);
g3 = plot(cvec,DP3(:,1),'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');
g4 = plot(cvec,DP4(:,1),'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');            
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
export_fig DPc -pdf -painters -transparent

figure(2)
lw = 3;
hold all
k1 = plot(rvec,DP32,'Color','k','LineStyle',':','LineWidth',lw);
k2 = plot(rvec,DP42,'Color','k','LineStyle','-','LineWidth',lw);
k3 = plot(cvec,DP3(:,2)./10,'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');
k4 = plot(cvec,DP4(:,2)./10,'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');            
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
export_fig DPd -pdf -painters -transparent