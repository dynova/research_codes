sigma_mu = 5e-2;
b1 = 2e-1;
b2_vec = linspace(1e-2, 1-1e-2, 1e2-1);
len_b2 = length(b2_vec);
d_vec = [0.01 1];
len_d = length(d_vec);
[acf1_x1, acf1_x2, v_x1, v_x2, cv_x1, cv_x2] = deal(zeros(len_b2,len_d));

for ind_b2 = 1:len_b2
    b2 = b2_vec(ind_b2);
    
    for ind_d = 1:len_d
        d = d_vec(ind_d);
        
        m1 = 1 + sqrt(b1);
        m2 = 1 + sqrt(b2);
        a_11 = b1 - 1 + 4*m1 - 3*m1^2 - d;
        a_22 = b2 - 1 + 4*m2 - 3*m2^2 - d;
        a_12 = d;
        a_21 = d;
        d_11 = sigma_mu^2*m1^2;
        d_22 = sigma_mu^2*m2^2;
        tr = a_11 + a_22;
        dt = a_11*a_22 - a_12*a_21;
        s1 = @(omega)(d_11*a_22^2+d_22*a_12^2+d_11*omega.^2)./((omega.^2-dt).^2+tr^2*omega.^2);
        s2 = @(omega)(d_22*a_11^2+d_11*a_21^2+d_22*omega.^2)./((omega.^2-dt).^2+tr^2*omega.^2);
        t1 = @(omega)cos(omega).*(d_11*a_22^2+d_22*a_12^2+d_11*omega.^2)./((omega.^2-dt).^2+tr^2*omega.^2);
        t2 = @(omega)cos(omega).*(d_22*a_11^2+d_11*a_21^2+d_22*omega.^2)./((omega.^2-dt).^2+tr^2*omega.^2);
        
        v_x1(ind_b2,ind_d) = (1/pi)*integral(s1,0,Inf,'AbsTol',1e-6,'RelTol',1e-3);
        
        v_x2(ind_b2,ind_d) = (1/pi)*integral(s2,0,Inf,'AbsTol',1e-6,'RelTol',1e-3);
        
        acf1_x1(ind_b2,ind_d) = (1/pi)*integral(t1,0,Inf,'AbsTol',1e-6,'RelTol',1e-3)./v_x1(ind_b2,ind_d);
        
        acf1_x2(ind_b2,ind_d) = (1/pi)*integral(t2,0,Inf,'AbsTol',1e-6,'RelTol',1e-3)./v_x2(ind_b2,ind_d);
        
        cv_x1(ind_b2,ind_d) = sqrt(v_x1(ind_b2,ind_d))/m1;
        
        cv_x2(ind_b2,ind_d) = sqrt(v_x2(ind_b2,ind_d))/m2;
    end
end

cases_x1 = repelem({acf1_x1, v_x1, cv_x1}, 1, len_d);
cases_x2 = repelem({acf1_x2, v_x2, cv_x2}, 1, len_d);

% figure('WindowState','fullscreen');
lw = 2;

ind_cases = 1;
subplot(3,2,ind_cases);
hold all
plot(b2_vec, cases_x1{ind_cases}(:,ind_cases),'b-','LineWidth',lw);
plot(b2_vec, cases_x2{ind_cases}(:,ind_cases),'k-','LineWidth',lw);
title('Low dispersal','FontWeight','bold','FontSize',14);
ylabel('Lag-1 Autocorrelation');
ylim([min([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)]) max([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)])]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 2;
subplot(3,2,ind_cases);
hold all
plot(b2_vec, cases_x1{ind_cases}(:,ind_cases),'b-','LineWidth',lw);
plot(b2_vec, cases_x2{ind_cases}(:,ind_cases),'k-','LineWidth',lw);
title('High dispersal','FontWeight','bold','FontSize',14);
ylim([min([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)]) max([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)])]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 3;
subplot(3,2,ind_cases);
hold all
plot(b2_vec, cases_x1{ind_cases}(:,ind_cases-2),'b-','LineWidth',lw);
plot(b2_vec, cases_x2{ind_cases}(:,ind_cases-2),'k-','LineWidth',lw);
ylabel('Variance');
ylim([min([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)]) max([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)])]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','off');
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 4;
subplot(3,2,ind_cases);
hold all
h1 = plot(b2_vec, cases_x1{ind_cases}(:,ind_cases-2),'b-','LineWidth',lw);
h2 = plot(b2_vec, cases_x2{ind_cases}(:,ind_cases-2),'k-','LineWidth',lw);
s = {' x_1',' x_2'};
f = legend([h1,h2],s);
set(f,'Location','NorthEast','box','off');
ylim([min([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)]) max([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)])]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','off');
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 5;
subplot(3,2,ind_cases);
hold all
plot(b2_vec, cases_x1{ind_cases}(:,ind_cases-4),'b-','LineWidth',lw);
plot(b2_vec, cases_x2{ind_cases}(:,ind_cases-4),'k-','LineWidth',lw);
ylabel('CV');
xlabel('{\boldmath$\beta_2$}','FontWeight','bold','interpreter','latex');
ylim([min([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)]) max([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)])]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','off');
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 6;
subplot(3,2,ind_cases);
hold all
plot(b2_vec, cases_x1{ind_cases}(:,ind_cases-4),'b-','LineWidth',lw);
plot(b2_vec, cases_x2{ind_cases}(:,ind_cases-4),'k-','LineWidth',lw);
xlabel('{\boldmath$\beta_2$}','FontWeight','bold','interpreter','latex');
ylim([min([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)]) max([cases_x1{ind_cases}(:); cases_x2{ind_cases}(:)])]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YScale','linear','YMinorTick','off');
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

export_fig Fig7 -pdf -painters -transparent