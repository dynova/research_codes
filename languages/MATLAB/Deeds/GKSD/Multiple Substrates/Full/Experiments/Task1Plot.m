syms a2 k2 Q1 Q2 d1 d2 KM1 KM2 DT
load Task1
n = 9;
Q2vec = logspace(-4,4,n);
rvec = logspace(-2,2,100);
% Empirical extraction of r50 %%
[~,x1] = min(abs(T12-mean([5018.5,50000])));
[~,x2] = min(abs(T22-mean([5018.5,50000])));
[~,x3] = min(abs(T32-mean([5018.4,50000])));
[~,x4] = min(abs(T42-mean([5017.8,50000])));
[~,x5] = min(abs(T52-mean([5014.8,50000])));
[~,x6] = min(abs(T62-mean([5012.1,50000])));
[~,x7] = min(abs(T72-mean([12587,50000])));
[~,x8] = min(abs(T82-mean([46222,50000])));
[~,x9] = min(abs(T92-mean([49622,50000])));
a2vec(1) = T13(x1)./T14(x1);
a2vec(2) = T23(x2)./T24(x2);
a2vec(3) = T33(x3)./T34(x3);
a2vec(4) = T43(x4)./T44(x4);
a2vec(5) = T53(x5)./T54(x5);
a2vec(6) = T63(x6)./T64(x6);
a2vec(7) = T73(x7)./T74(x7);
a2vec(8) = T83(x8)./T84(x8);
a2vec(9) = T93(x9)./T94(x9);
% a2vec = linspace(0.01,1,n);
%% Analytical expression with a2 taking values from simulation
r50expr = k2.^(-1).*(Q2.*a2.*(d1+a2.*((-1).*d1+d2)).^(-1)+d1.*(d1+d2).^(-1).*(d1+d1.*((-1).*d1+d2).*(d1+d2).^(-1)).^(-1).*Q1).*(1+KM1.*(Q2.*(1+(-1).*a2).*(d1+a2.*((-1).*d1+d2)).^(-1)+(1/2).*d1.^(-1).*Q1).^(-1)).*(d2.*DT.^(-1)+(d2+k2).*(Q2.*a2.*(d1+a2.*((-1).*d1+d2)).^(-1)+KM2+(1/2).*d2.^(-1).*Q1).^(-1));
for i = 1:length(Q2vec)
%     for j = 1:length(a2vec)
        r50analytical(i) = double(subs(r50expr,{a2,k2,Q1,Q2,d1,d2,KM1,KM2,DT},{a2vec(i),0.999,1,Q2vec(i),2e-5,2e-4,1000.02,1000.2,1}));
%     end
end
% Simulation r50
r50simulation(1) = rvec(x1);
r50simulation(2) = rvec(x2);
r50simulation(3) = rvec(x3);
r50simulation(4) = rvec(x4);
r50simulation(5) = rvec(x5);
r50simulation(6) = rvec(x6);
r50simulation(7) = rvec(x7);
r50simulation(8) = rvec(x8);
r50simulation(9) = rvec(x9);

hold all
g1 = plot(Q2vec,r50simulation,'bo','MarkerSize',10,'MarkerEdgeColor','none','MarkerFaceColor','b');
g2 = plot(Q2vec,r50analytical,'r','LineWidth',3);            
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','XScale','log','LineWidth',3);
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
xlabel('Synthesis rate of Substrate 2 (Q_2) [nM s^{-1}]','FontWeight','bold','FontSize',20);
ylabel('r_{50}([S_1]_T)','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.08, 0.5, 0]);
s = {'Simulation','Semi-Analytical'};
legend(s,'FontSize',14,'FontWeight','bold','Location','Northwest');
box off
hold off
export_fig r50S1T -pdf -painters -transparent