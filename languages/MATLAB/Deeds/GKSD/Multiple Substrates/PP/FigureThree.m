load PPM.mat
load PPC.mat
Q2vec = 2*logspace(-4,0,100);
Pvec = Q2vec(1:10:100);
figure(1)
lw = 3;
hold all
g1 = plot(Q2vec,PP1./PP2,'Color','b','LineStyle','-','LineWidth',lw);
g2 = plot(Pvec,PP3(:,1),'bo','MarkerSize',10);
hold off
set(gcf,'units','inches','PaperPosition',[0 0 2.26 1.59])
ax = gca;
set(ax,'FontSize',14,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XMinorTick','on','YTick',...
    0.1:0.1:1,'XScale','log');
ylim([0 1]);
ax.XAxis.TickLength = [0.025 0.0375];
xlabel('[Q_2]','FontWeight','bold','FontSize',18);
ylabel('Normalized Fraction S_1^*','FontWeight','bold','FontSize',18,...
    'Units','Normalized','Position',[-0.09, 0.5, 0]);
title('Processive-Processive Model','FontWeight','bold','FontSize',20)
box off
export_fig PPA -pdf -painters -transparent

figure(2)
lw = 3;
hold all
h1 = plot(Q2vec,PP2,'Color','b','LineStyle','-','LineWidth',lw);
h2 = plot(Pvec,PP3(:,2)./10,'bo','MarkerSize',10);
hold off
set(gcf,'units','inches','PaperPosition',[0 0 2.26 1.59])
ax = gca;
set(ax,'FontSize',14,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XMinorTick','on',...
    'YMinorTick','on','XScale','log','YScale','log');
ax.TickLength = [0.025 0.0375];
ylim([1e1 1e2]);
xlabel('[Q_2]','FontWeight','bold','FontSize',18);
ylabel('Total Substrate [S_1]_T [nM]','FontWeight','bold','FontSize',18,...
    'Units','Normalized','Position',[-0.06, 0.5, 0]);
title('Processive-Processive Model','FontWeight','bold','FontSize',20)
box off
export_fig PPB -pdf -painters -transparent