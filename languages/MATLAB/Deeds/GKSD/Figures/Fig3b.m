syms Q r1 r2 d a k2 Pt S00
x = 100;
Qvec = logspace(-1,5,x);
HillF = log(81)/(log(81)+log(((Q+10*((d+k2)/a)*r1)*(Q+10*((d+k2)/a)*r2))/((9*Q+10*((d+k2)/a)*r1)*(9*Q+10*((d+k2)/a)*r2))) + log(((10*k2*Pt+9*Q+10*(((d+k2)/a)+Pt)*r2)/(10*k2*Pt+Q+10*(((d+k2)/a)+Pt)*r2))));
for kk = 1:x
    HillF_Q(kk) = subs(HillF,{Q,r1,r2,d,a,k2,Pt},...
        {Qvec(kk),2e-5,2e-4,1e-3,1e-4,1-1e-3,0.1});
end 
figure(1)
hold all
plot(Qvec,HillF_Q,'k','LineSmoothing','on','Linewidth',3);
xlim([1e-1 1e5])
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XScale','log',...
    'XMinorTick','on','XTick',[1e-1 1e2 1e5]);
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = [linspace(1e-1,1e2,10) ...
    linspace(11200,1e5,9)];
ylabel('Hill coefficient (n_{eff})','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.08, 0.5, 0]);
xlabel('Synthesis Rate (Q) [nM s^{-1}]','FontWeight','bold','FontSize',20);
box off
export_fig Fig3b -pdf -painters -transparent