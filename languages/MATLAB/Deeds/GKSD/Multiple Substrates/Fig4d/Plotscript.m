files = dir('*.mat');
for i = 1:numel(files)
    load(files(i).name);
end
Q2vec = 2*logspace(-4,4,100);
rvec = logspace(-2,4,1000);
w = zeros(1,length(Q2vec));
x = zeros(1,length(Q2vec));
y = zeros(1,length(Q2vec));
z = zeros(1,length(Q2vec));
r50w = zeros(1,length(Q2vec));
r50x = zeros(1,length(Q2vec));

for i = 1:length(Q2vec)
    [~,w(i)] = ...
    min(abs(eval(['F' num2str(i) '1'])-(0.5*min(eval(['F' num2str(i) '1']))+25000)));
    r50w(i) = rvec(w(i));
end
for i = 1:length(Q2vec)
    [~,x(i)] = ...
    min(abs(eval(['R' num2str(i) '1'])-(0.5*min(eval(['R' num2str(i) '1']))+0.5*max(eval(['R' num2str(i) '1'])))));
    r50x(i) = rvec(x(i));
end
for i = 1:length(Q2vec)
    y(i) = eval(['F' num2str(i) '2(x(i))']);
end
for i = 1:length(Q2vec)
    z(i) = 0.1001E1.*(0.25E4+(0.2E-4+0.18E-3.*y(i)).^(-1).*y(i).*Q2vec(i)).*...
        (1+1000.*(0.25E5+(1+(-1).*y(i)).*(0.2E-4+0.18E-3.*y(i)).^...
        (-1).*Q2vec(i)).^(-1)).*(0.2E-3+0.9992E0.*(0.35E4+(0.2E-4+0.18E-3.*...
        y(i)).^(-1).*y(i).*Q2vec(i)).^(-1));
end
hold all
p1 = plot(Q2vec(1:3:100),r50w(1:3:100),'ko','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','k');
p2 = plot(Q2vec,z,'Color','k','LineStyle','-','LineWidth',3);
p3 = plot(Q2vec,r50x,'Color','b','LineStyle',':','LineWidth',3);
p4 = semilogx(NaN, NaN,'Color',[0.5 0.5 0.5],'LineStyle',':','LineWidth',lw);
p5 = semilogx(NaN, NaN,'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',lw);
p6 = semilogx(NaN, NaN,'Color',[1 1 1],'LineStyle',':','LineWidth',lw);
p7 = semilogx(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
xlim([1e-4 1e2]);
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','XScale','log','YScale','log','LineWidth',3);
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
xlabel('Synthesis rate of Substrate 2 (Q_2) [nM s^{-1}]','FontWeight','bold','FontSize',20);
ylabel('r_{50}([S_1]_T)','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.08, 0.5, 0]);
s = {'Simulation','Semi-Analytical','','','Full Model','Proc. E3, Seq. DUB Model'};
f = legend([p1 p2 p3 p6 p7 p4 p5],s);
set(f,'FontWeight','bold','FontSize',16,'Location','Best','box','off');
box off
hold off
export_fig Fig4d -pdf -painters -transparent