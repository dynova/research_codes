sigma_mu = 5e-2;
b_vec = linspace(1e-2, 1-1e-2, 1e2-1);
len_b = length(b_vec);
d_vec = [0 0.01 1];
len_d = length(d_vec);
[acf1_mu, v_mu, cv_mu] = deal(zeros(len_b,len_d));

for ind_b = 1:len_b
    b = b_vec(ind_b);
    
    for ind_d = 1:len_d
        d = d_vec(ind_d);
        
        acf1_mu(ind_b,ind_d) = ...
            ((d * exp(2*d) + sqrt(b) + exp(2*d) * sqrt(b) + b + exp(2*d) * b) ...
            / (d + 2 * sqrt(b)*(sqrt(b) + 1)) * exp(-2*(d + sqrt(b)*(sqrt(b)+1))));
        
        v_mu(ind_b,ind_d) = ...
            (sigma_mu^2 * (sqrt(b) + 1) * (2 * sqrt(b) * (sqrt(b) + 1) + d) ...
            / (8 * sqrt(b) * (d + sqrt(b) * (sqrt(b) + 1))));
        
        cv_mu(ind_b,ind_d) = ...
            (1+sqrt(b)) * (sigma_mu * sqrt(2) / (4 * (1 + sqrt(b))) ...
            * sqrt((d + 2 * sqrt(b) * (sqrt(b) + 1)) ...
            / (sqrt(b) * (sqrt(b) + 1) * (d + sqrt(b) * (sqrt(b) + 1)))));
    end
end

cases_mu = repelem({acf1_mu, v_mu, cv_mu}, 1, len_d);

lw = 2;

ind_cases = 1;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases),'k-','LineWidth',lw);
title('Isolated','FontWeight','bold','FontSize',14);
ylabel('Lag-1 Autocorrelation','FontSize',10);
ytickformat('%.1f')
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];

ind_cases = 2;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases),'k-','LineWidth',lw);
ytickformat('%.1f')
title('Low dispersal','FontWeight','bold','FontSize',14);
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];

ind_cases = 3;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases),'k-','LineWidth',lw);
ytickformat('%.1f')
title('High dispersal','FontWeight','bold','FontSize',14);
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 4;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases-3),'k-','LineWidth',lw);
ytickformat('%.1f')
ylabel('Variance','FontSize',10);
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 5;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases-3),'k-','LineWidth',lw);
ytickformat('%.1f')
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 6;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases-3),'k-','LineWidth',lw);
ytickformat('%.1f')
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 7;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases-6),'k-','LineWidth',lw);
xlabel('{\boldmath$\beta$}','FontWeight','bold','interpreter','latex');
ylabel('CV','FontSize',10);
ytickformat('%.2f')
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 8;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases-6),'k-','LineWidth',lw);
xlabel('{\boldmath$\beta$}','FontWeight','bold','interpreter','latex');
ytickformat('%.2f')
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 9;
subplot(3,3,ind_cases);
hold all
plot(b_vec, cases_mu{ind_cases}(:,ind_cases-6),'k-','LineWidth',lw);
xlabel('{\boldmath$\beta$}','FontWeight','bold','interpreter','latex');
ytickformat('%.2f')
ylim([min(cases_mu{ind_cases}(:)) max(cases_mu{ind_cases}(:))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

export_fig Fig2 -pdf -painters -transparent