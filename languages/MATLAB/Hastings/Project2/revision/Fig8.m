time_max = 1.1e3;
num = 502;
x = linspace(time_max/2,time_max,num);

figure('WindowState','fullscreen');
lw = 3;
ms = 1;

ind_cases = 1;
subplot(3,2,ind_cases);
hold all
load Fig8_12.mat;
h1 = plot(x, acf1_lb,'bo','MarkerSize',ms);
h2 = plot(x, acf1_mdn,'b-','LineWidth',lw);
h3 = plot(x, acf1_ub,'bo','MarkerSize',ms);
h4 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t1 = acf1_t;
load Fig8_22.mat;
h5 = plot(x, acf1_lb,'ko','MarkerSize',ms);
h6 = plot(x, acf1_mdn,'k-','LineWidth',lw);
h7 = plot(x, acf1_ub,'ko','MarkerSize',ms);
h8 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t2 = acf1_t;
title('Low dispersal','FontWeight','bold','FontSize',20);
ylabel('Lag-1 Autocorrelation','FontSize',18);
s = {' x_1', [' \tau = ' num2str(t1)],' x_2',[' \tau = ' num2str(t2)]};
f = legend([h2,h4,h6,h8],s);
set(f,'Location','NorthWest','box','off','FontSize',18);
xlim([0 time_max]);
ylim([0 1]);
xticks(0:200:1000);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'XMinorTick','off','FontSize',18);
ax.XAxis.TickLength = [0.020 0.020];
ax.YAxis.TickLength = [0.020 0.020];

ind_cases = 2;
subplot(3,2,ind_cases);
hold all
load Fig8_13.mat;
h1 = plot(x, acf1_lb,'bo','MarkerSize',ms);
h2 = plot(x, acf1_mdn,'b-','LineWidth',lw);
h3 = plot(x, acf1_ub,'bo','MarkerSize',ms);
h4 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t1 = acf1_t;
load Fig8_23.mat;
h5 = plot(x, acf1_lb,'ko','MarkerSize',ms);
h6 = plot(x, acf1_mdn,'k-','LineWidth',lw);
h7 = plot(x, acf1_ub,'ko','MarkerSize',ms);
h8 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t2 = acf1_t;
title('High dispersal','FontWeight','bold','FontSize',20);
s = {' x_1', [' \tau = ' num2str(t1)],' x_2',[' \tau = ' num2str(t2)]};
f = legend([h2,h4,h6,h8],s);
set(f,'Location','NorthWest','box','off','FontSize',18);
xlim([0 time_max]);
ylim([0 1]);
xticks(0:200:1000);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'XMinorTick','off','FontSize',18);
ax.XAxis.TickLength = [0.020 0.020];
ax.YAxis.TickLength = [0.020 0.020];

ind_cases = 3;
subplot(3,2,ind_cases);
hold all
load Fig8_12.mat;
h1 = plot(x, vr_lb,'bo','MarkerSize',ms);
h2 = plot(x, vr_mdn,'b-','LineWidth',lw);
h3 = plot(x, vr_ub,'bo','MarkerSize',ms);
h4 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t1 = vr_t;
load Fig8_22.mat;
h5 = plot(x, vr_lb,'ko','MarkerSize',ms);
h6 = plot(x, vr_mdn,'k-','LineWidth',lw);
h7 = plot(x, vr_ub,'ko','MarkerSize',ms);
h8 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t2 = vr_t;
ylabel('Variance','FontSize',18);
s = {' x_1', [' \tau = ' num2str(t1)],' x_2',[' \tau = ' num2str(t2)]};
f = legend([h2,h4,h6,h8],s);
set(f,'Location','NorthWest','box','off','FontSize',18);
xlim([0 time_max]);
ylim([0 2e-2]);
xticks(0:200:1000);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'XMinorTick','off','FontSize',18);
ax.XAxis.TickLength = [0.020 0.020];
ax.YAxis.TickLength = [0.020 0.020];

ind_cases = 4;
subplot(3,2,ind_cases);
hold all
load Fig8_13.mat;
h1 = plot(x, vr_lb,'bo','MarkerSize',ms);
h2 = plot(x, vr_mdn,'b-','LineWidth',lw);
h3 = plot(x, vr_ub,'bo','MarkerSize',ms);
h4 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t1 = vr_t;
load Fig8_23.mat;
h5 = plot(x, vr_lb,'ko','MarkerSize',ms);
h6 = plot(x, vr_mdn,'k-','LineWidth',lw);
h7 = plot(x, vr_ub,'ko','MarkerSize',ms);
h8 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t2 = vr_t;
s = {' x_1', [' \tau = ' num2str(t1)],' x_2',[' \tau = ' num2str(t2)]};
f = legend([h2,h4,h6,h8],s);
set(f,'Location','NorthWest','box','off','FontSize',18);
xlim([0 time_max]);
ylim([0 2e-2]);
xticks(0:200:1000);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'XMinorTick','off','FontSize',18);
ax.XAxis.TickLength = [0.020 0.020];
ax.YAxis.TickLength = [0.020 0.020];

ind_cases = 5;
subplot(3,2,ind_cases);
hold all
load Fig8_12.mat;
h1 = plot(x, cv_lb,'bo','MarkerSize',ms);
h2 = plot(x, cv_mdn,'b-','LineWidth',lw);
h3 = plot(x, cv_ub,'bo','MarkerSize',ms);
h4 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t1 = cv_t;
load Fig8_22.mat;
h5 = plot(x, cv_lb,'ko','MarkerSize',ms);
h6 = plot(x, cv_mdn,'k-','LineWidth',lw);
h7 = plot(x, cv_ub,'ko','MarkerSize',ms);
h8 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t2 = cv_t;
xlabel('Time','FontSize',18);
ylabel('CV','FontSize',18);
s = {' x_1', [' \tau = ' num2str(t1)],' x_2',[' \tau = ' num2str(t2)]};
f = legend([h2,h4,h6,h8],s);
set(f,'Location','NorthWest','box','off','FontSize',18);
xlim([0 time_max]);
ylim([0 0.15]);
xticks(0:200:1000);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'XMinorTick','off','FontSize',18);
ax.XAxis.TickLength = [0.020 0.020];
ax.YAxis.TickLength = [0.020 0.020];

ind_cases = 6;
subplot(3,2,ind_cases);
hold all
load Fig8_13.mat;
h1 = plot(x, cv_lb,'bo','MarkerSize',ms);
h2 = plot(x, cv_mdn,'b-','LineWidth',lw);
h3 = plot(x, cv_ub,'bo','MarkerSize',ms);
h4 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t1 = cv_t;
load Fig8_23.mat;
h5 = plot(x, cv_lb,'ko','MarkerSize',ms);
h6 = plot(x, cv_mdn,'k-','LineWidth',lw);
h7 = plot(x, cv_ub,'ko','MarkerSize',ms);
h8 = plot(NaN, NaN,'Color',[1 1 1],'LineStyle','-','LineWidth',lw);
t2 = cv_t;
xlabel('Time','FontSize',18);
s = {' x_1', [' \tau = ' num2str(t1)],' x_2',[' \tau = ' num2str(t2)]};
f = legend([h2,h4,h6,h8],s);
set(f,'Location','NorthWest','box','off','FontSize',18);
xlim([0 time_max]);
ylim([0 0.15]);
xticks(0:200:1000);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'XMinorTick','off','FontSize',18);
ax.XAxis.TickLength = [0.020 0.020];
ax.YAxis.TickLength = [0.020 0.020];

set(gcf,'PaperSize',[20 12]);
export_fig Fig8 -pdf -painters -transparent