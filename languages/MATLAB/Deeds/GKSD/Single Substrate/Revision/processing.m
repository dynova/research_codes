load DDS.mat
x1 = squeeze(out_r50(:,1));
x2 = squeeze(out_r50(:,2));
y1 = squeeze(out_neff(:,1));
y2 = squeeze(out_neff(:,2));
load PDS.mat
x3 = squeeze(out_r50(:,1));
x4 = squeeze(out_r50(:,2));
y3 = squeeze(out_neff(:,1));
y4 = squeeze(out_neff(:,2));
load DDT.mat
x5 = squeeze(out_r50(:,1));
x6 = squeeze(out_r50(:,2));
y5 = squeeze(out_neff(:,1));
y6 = squeeze(out_neff(:,2));
load PDT.mat
x7 = squeeze(out_r50(:,1));
x8 = squeeze(out_r50(:,2));
y7 = squeeze(out_neff(:,1));
y8 = squeeze(out_neff(:,2));
load DP.mat
x9 = squeeze(out_r50(:,1));
x10 = squeeze(out_r50(:,2));
y9 = squeeze(out_neff(:,1));
y10 = squeeze(out_neff(:,2));
load PP.mat
x11 = squeeze(out_r50(:,1));
x12 = squeeze(out_r50(:,2));
y11 = squeeze(out_neff(:,1));
y12 = squeeze(out_neff(:,2));

temp1 = horzcat([x1(:); x3(:); x5(:); x7(:); x9(:); x11(:)]);
temp2 = horzcat([x2(:); x4(:); x6(:); x8(:); x10(:); x12(:)]);
temp3 = horzcat([y1(:); y3(:); y5(:); y7(:); y9(:); y11(:)]);
temp4 = horzcat([y2(:); y4(:); y6(:); y8(:); y10(:); y12(:)]);

x = [temp1;temp2];
y = [temp3;temp4];
figure(1)
hold all
histogram(x,'normalization','probability');
title('r_{50} ratios','FontSize',16)
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold');
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
xlabel('Ratio','FontSize',14)
ylabel('Probability','FontSize',14)
xlim([0 4]);
ylim([0 0.2]);
hold off
export_fig r50_ratios -pdf -painters -transparent

figure(2)
hold all
histogram(y,'normalization','probability');
title('n_{eff} ratios','FontSize',16)
ax = gca;
set(ax,'FontSize',16,'FontWeight','bold');
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
xlabel('Ratio','FontSize',14)
ylabel('Probability','FontSize',14)
xlim([0 6]);
ylim([0 0.2]);
hold off
export_fig neff_ratios -pdf -painters -transparent