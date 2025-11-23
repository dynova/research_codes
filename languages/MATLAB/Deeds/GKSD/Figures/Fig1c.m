x = [0 0 0 8 7 25 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1 3 2 1 1 1 1 5 8 5 1 8 0 2 1 0 0 6 38 0 0 0 5 35 1 1 0 2 0 0 1 zeros(1,17) 2 7 3 1 1 0 0 0 0 0 0 4 1 10 1 1 1 1 4 1 1 1 10 2 0 0 2 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 zeros(1,40) 0 1 1 0 0 2 1 0 0 1 0 0 1 0 1 0 2 1 1 19 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 2 0 1 0 2 5 1 20 2 0 1 4 1 1 18 16 0 1 23 2 0 0 5 3 5 0 0 0 1 1 1 2 20 27 0 0 0 0 3 1 0 1 2 0 0 1 0 2 1 0 2 14 13 0 1 0 2 0 3 2 0 0 0 0 1 0 1 2 1 1 0 3 1 2 3 0 1 3 0 0 0 0 0 0 1 0 1 3 1 1 0 2 0 4 0 0 0 1 1 4 2 3 5 2 3 7 25 2 2 1 1 1 2 1 1 1 1 1 1 2 7 1 1 12 1 24 0 3 0 1 1 20 10 92 11 0 0 0 0 0 0 0 0 12 0 0 0 1 2 2 1 1 0 12 0 1 0 1 1 1 1 3 1 5 1 0 1 0 1 0 0 0 1 2 0 0 0 1 0 0 0 0 1 3 0 0 0 5 10 0 1 0 3 0 0 0 0 2 0 0 0 1 0 3 2 0 0 12 3 8 0 0 1 6 0 0 0 1 0 0]';
num_zeros = length(x)-length(find(x)); % 225
y = x(x>0);
pmf = hist(y,max(y))'./numel(y);
figure(1)
set(gcf,'PaperPosition', [0,0,2.26,1.59])
plot(pmf,'bo','LineSmoothing','on','LineWidth',3)
ax = gca;
set(ax,'FontSize',15,'FontWeight','bold','LineWidth',4,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'YScale','log',...
    'YTick',[1e-2,1e-1,1e0],'XScale','log','XTick',[1e1,1e2]);
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
xlabel('Number of Substrates per E3 Ligase','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[0.5, -0.09, 0]);
ylabel('Probability','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.09, 0.5, 0]);
box off
% export_fig Fig1c -pdf -painters -transparent