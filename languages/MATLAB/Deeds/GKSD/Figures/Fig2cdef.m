syms Q r1 r2 d a k2 Pt S00
x = 100;
Qvec = logspace(-1,5,x);
avec = logspace(-1,5,x);
KMvec = logspace(1,-5,x);
InvSat = logspace(-2,-8,x);
EC50GK = ones(1,x);
EC50I = 1 + Q/(2*k2*Pt)+r1*(((d+k2)/a)+Pt)/(k2*Pt);
EC50F = 1 + Q*r2/((r1+r2)*k2*Pt)+r2*(((d+k2)/a)+Pt)/(k2*Pt);
HillI = log(81)/(log(81)+2*log((((d+k2)/a)/S00+0.1)/(((d+k2)/a)/S00+0.9))+...
    log(((k2+r1)*Pt/Q+((d+k2)/a)/S00+0.9)/((k2+r1)*Pt/Q+((d+k2)/a)/S00+0.1)));
HillF = log(81)/(log(81)+log(((9*r1+r2)*(Q+((d+k2)/a)*(9*r1+r2))*...
    (Q+((d+k2)/a)*(9*r2+r1))*(r2*(9*Q+(9*r2+r1)*(((d+k2)/a)+Pt))+Pt*k2*(9*r2+r1)))/...
    ((9*r2+r1)*(9*Q+((d+k2)/a)*(9*r1+r2))*(9*Q+((d+k2)/a)*(9*r2+r1))*...
    (r2*(Q+(9*r1+r2)*(((d+k2)/a)+Pt))+Pt*k2*(9*r1+r2)))));

old_a = 1e-2;
new_a = 1e-4;
for kk = 1:x
    HillGK(kk) = log(81)./(log(81)+2*log((InvSat(kk)+0.1)./(InvSat(kk)+0.9)));
    EC50I_Q(kk) = subs(EC50I,{Q,r1,r2,d,a,k2,Pt},...
        {Qvec(kk),2e-5,2e-5,1e-3,new_a,1-1e-3,.1});
    EC50F_Q(kk) = subs(EC50F,{Q,r1,r2,d,a,k2,Pt},...
        {Qvec(kk),2e-5,2e-4,1e-3,new_a,1-1e-3,.1});
    HillI_Q(kk) = subs(HillI,{Q,r1,r2,d,a,k2,Pt,S00},...
        {Qvec(kk),2e-5,2e-5,1e-3,new_a,1-1e-3,.1,Qvec(kk)./(2e-5)});
    HillF_Q(kk) = subs(HillF,{Q,r1,r2,d,a,k2,Pt},...
        {Qvec(kk),2e-5,2e-4,1e-3,new_a,1-1e-3,.1});
    EC50I_KM(kk) = subs(EC50I,{Q,r1,r2,d,a,k2,Pt},...
        {2e-2,2e-5,2e-5,1e-3,avec(kk),1-1e-3,.1});
    EC50F_KM(kk) = subs(EC50F,{Q,r1,r2,d,a,k2,Pt},...
        {2e-2,2e-5,2e-4,1e-3,avec(kk),1-1e-3,.1});
    HillI_KM(kk) = subs(HillI,{Q,r1,r2,d,a,k2,Pt,S00},...
        {2e-2,2e-5,2e-5,1e-3,avec(kk),1-1e-3,.1,1000});
    HillF_KM(kk) = subs(HillF,{Q,r1,r2,d,a,k2,Pt},...
        {2e-2,2e-5,2e-4,1e-3,avec(kk),1-1e-3,.1});
end

figure(1)
hold all
plot(Qvec,EC50GK,'b','LineSmoothing','on','Linewidth',3);
plot(Qvec,EC50I_Q,'r','LineSmoothing','on','Linewidth',3);
plot(Qvec,EC50F_Q,'k','LineSmoothing','on','Linewidth',3);
xlim([1e-1 1e5]);
hold off
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'XScale','log',...
    'XMinorTick','on','YScale','log','XTick',[1e-1 1e2 1e5],...
    'YTick',[1e0 1e3 1e6]);
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = [linspace(1e-1,1e2,10) ...
    linspace(11200,1e5,9)];
ax.YAxis.MinorTickValues = [linspace(1e0,1e3,10) ...
    linspace(112e3,1e6,9)];
ylabel('r_{50}','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.09, 0.5, 0]);
xlabel('Total Substrate ([S]_T) [nM]','FontWeight','bold','FontSize',20);
s= {' GK', ' Intermediate',' Full'};
legend(s,'FontWeight','bold','FontSize',16,...
    'Location','Best',...
    'EdgeColor',[1 1 1]);
box off
% export_fig Fig2c -pdf -painters -transparent

% figure(2)
% hold all
% plot(Qvec,HillGK,'b','LineSmoothing','on','Linewidth',3);
% plot(Qvec,HillI_Q,'r','LineSmoothing','on','Linewidth',3);
% plot(Qvec,HillF_Q,'k','LineSmoothing','on','Linewidth',3);
% hold off
% ax = gca;
% set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
%     'PlotBoxAspectRatio',[1.79/1.25 1 1],'YScale','log',...
%     'XScale','log','XMinorTick','on','YMinorTick','on',...
%     'YTick',[1e2 1e4 1e6 1e8],'XTick',[1e-1 1e2 1e5]);
% xlim([1e-1 1e5]);
% ax.XAxis.TickLength = [0.025 0.0375];
% ax.YAxis.TickLength = [0.025 0.0375];
% ax.XAxis.MinorTick = 'on';
% ax.YAxis.MinorTick = 'on';
% ax.XAxis.MinorTickValues = [linspace(1e-1,1e2,10) ...
%     linspace(112e2,1e5,9)];
% ax.YAxis.MinorTickValues = [linspace(1e0,1e2,10) ...
%     linspace(12e2,1e4,9) linspace(12e4,1e6,9) linspace(12e6,1e8,9)];
% xlabel('Synthesis Rate (Q) [nM s^{-1}]','FontWeight','bold','FontSize',20);
% ylabel('Hill coefficient (n_{eff})','FontWeight','bold','FontSize',20,...
%     'Units','Normalized','Position',[-0.09, 0.5, 0]);
% s= {' GK', ' Intermediate',' Full'};
% legend(s,'Position',[0.201785714285714 0.490562649640863 0.228571428571429 0.266666666666667],...
%     'FontSize',16,'box','off');
% box off
% export_fig Fig2f -pdf -painters -transparent

% figure(3)
% hold all
% plot(KMvec,EC50GK,'b','LineSmoothing','on','Linewidth',3);
% plot(KMvec,EC50I_KM,'r','LineSmoothing','on','Linewidth',3);
% plot(KMvec,EC50F_KM,'k','LineSmoothing','on','Linewidth',3);
% hold off
% ax = gca;
% set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
%     'PlotBoxAspectRatio',[1.79/1.25 1 1],'YScale','log',...
%     'XScale','log','XMinorTick','on',...
%     'XTick',[1e-5 1e-2 1e1]);
% xlim([1e-5 1e1]);
% ax.XAxis.TickLength = [0.025 0.0375];
% ax.YAxis.TickLength = [0.025 0.0375];
% ax.XAxis.MinorTick = 'on';
% ax.YAxis.MinorTick = 'on';
% ax.XAxis.MinorTickValues = [linspace(1e-5,1e-2,10) ...
%     linspace(1.12,1e1,9)];
% xlabel('K_M [nM]','FontWeight','bold','FontSize',20);
% ylabel('r_{50}','FontWeight','bold','FontSize',20,...
%     'Units','Normalized','Position',[-0.09, 0.5, 0]);
% s= {' GK',' Intermediate',' Full'};
% legend(s,'FontWeight','bold','FontSize',16,...
%     'Location','Best',...
%     'EdgeColor',[1 1 1]);
% box off
% export_fig Fig2e -pdf -painters -transparent
% 
figure(4)
hold all
plot(KMvec,HillGK,'b','LineSmoothing','on','Linewidth',3);
plot(KMvec,HillI_KM,'r','LineSmoothing','on','Linewidth',3);
plot(KMvec,HillF_KM,'k','LineSmoothing','on','Linewidth',3);
hold off
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold','LineWidth',3,...
    'PlotBoxAspectRatio',[1.79/1.25 1 1],'YScale','log',...
    'XScale','log','XMinorTick','on','YMinorTick','on',...
    'YTick',[1e2 1e4 1e6 1e8],'XTick',[1e-5 1e-2 1e1]);
xlim([1e-5 1e1]);
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = [linspace(1e-5,1e-2,10) ...
    linspace(1.12,1e1,9)];
ax.YAxis.MinorTickValues = [linspace(1e0,1e2,10) ...
    linspace(12e2,1e4,9) linspace(12e4,1e6,9) linspace(12e6,1e8,9)];
xlabel('K_M [nM]','FontWeight','bold','FontSize',20);
ylabel('Hill coefficient (n_{eff})','FontWeight','bold','FontSize',20,...
    'Units','Normalized','Position',[-0.09, 0.5, 0]);
s= {' GK',' Intermediate',' Full'};
legend(s,'FontWeight','bold','FontSize',16,...
    'Location','Best',...
    'EdgeColor',[1 1 1]);
box off
export_fig Fig2d -pdf -painters -transparent