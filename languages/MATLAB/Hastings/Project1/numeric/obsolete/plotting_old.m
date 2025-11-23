load running.mat;
N = 25;
arr_b1 = linspace(0.01,0.99,3);
arr_b2 = linspace(0.01,0.99,3);
arr_d = logspace(-2,2,1000);
arr_h = linspace(2,N,N-1);
len_b1 = length(arr_b1);
len_b2 = length(arr_b2);
len_d = length(arr_d);
len_h = length(arr_h);

figure('WindowState','fullscreen');
lw = 2;

ind_cases = 1;
ind_b1 = 1;
ind_b2 = 1;
arr_hh = squeeze(hh(:,ind_b1,ind_b2,2:11:24))';
subplot(3,3,ind_cases);
hold all
plot(arr_d,arr_hh,'Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ylabel('Prob. of system recovery');
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
hold off

ind_cases = 2;
ind_b1 = 1;
ind_b2 = 2;
arr_hh = squeeze(hh(:,ind_b1,ind_b2,2:11:24))';
subplot(3,3,ind_cases);
hold all
plot(arr_d,arr_hh,'Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
hold off

ind_cases = 3;
ind_b1 = 1;
ind_b2 = 3;
arr_hh = squeeze(hh(:,ind_b1,ind_b2,2:11:24))';
subplot(3,3,ind_cases);
hold all
plot(arr_d,arr_hh,'Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
hold off

ind_cases = 4;
ind_b1 = 2;
ind_b2 = 1;
arr_hh = squeeze(hh(:,ind_b1,ind_b2,2:11:24))';
subplot(3,3,ind_cases);
hold all
plot(arr_d,arr_hh,'Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ylabel('Prob. of system recovery');
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
hold off

ind_cases = 5;
ind_b1 = 2;
ind_b2 = 2;
arr_hh = squeeze(hh(:,ind_b1,ind_b2,2:11:24))';
subplot(3,3,ind_cases);
hold all
plot(arr_d,arr_hh,'Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
hold off

ind_cases = 6;
ind_b1 = 2;
ind_b2 = 3;
subplot(3,3,ind_cases);
hold all
h1 = plot(arr_d,squeeze(hh(:,ind_b1,ind_b2,2))','Linewidth',lw);
h2 = plot(arr_d,squeeze(hh(:,ind_b1,ind_b2,13))','Linewidth',lw);
h3 = plot(arr_d,squeeze(hh(:,ind_b1,ind_b2,24))','Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
s = {' H = 3',' H = 14',' H = 25'}; 
f = legend([h1,h2,h3],s);
set(f,'Location','Best','box','off');
hold off

ind_cases = 7;
ind_b1 = 3;
ind_b2 = 1;
arr_hh = squeeze(hh(:,ind_b1,ind_b2,2:11:24))';
subplot(3,3,ind_cases);
hold all
plot(arr_d,arr_hh,'Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlabel('Dispersal rate (D)');
ylabel('Prob. of system recovery');
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
hold off

ind_cases = 8;
ind_b1 = 3;
ind_b2 = 2;
arr_hh = squeeze(hh(:,ind_b1,ind_b2,2:11:24))';
subplot(3,3,ind_cases);
hold all
plot(arr_d,arr_hh,'Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlabel('Dispersal rate (D)');
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
hold off

ind_cases = 9;
ind_b1 = 3;
ind_b2 = 3;
arr_hh = squeeze(hh(:,ind_b1,ind_b2,2:11:24))';
subplot(3,3,ind_cases);
hold all
plot(arr_d,arr_hh,'Linewidth',lw);
title(['\beta_1 = ' num2str(arr_b1(ind_b1)) ...
    ', \beta_2 = ' num2str(arr_b2(ind_b2))],'FontWeight','bold',...
    'FontSize',18);
ax = gca;
set(ax,'Xscale','Log','FontWeight','bold','LineWidth',lw,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlabel('Dispersal rate (D)');
xlim([0 100]);
ylim([0 1]);
yline(0.5,'k--','Linewidth',lw);
hold off

print(gcf,'Fig3.png','-dpng','-r300');