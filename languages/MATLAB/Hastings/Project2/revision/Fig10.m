sigma_mu = 5e-2;
b2_vec = linspace(0.01,0.99,3);
b1_vec = linspace(0.01,0.99,99);
d_vec = [0.01 0.1 1];
len_b1 = length(b1_vec);
len_b2 = length(b2_vec);
len_d = length(d_vec);
[acf1_x1, acf1_x2, v_x1, v_x2, cv_x1, cv_x2] = deal(zeros(len_b1,len_b2,len_d));

for ind_b1 = 1:len_b1
    b1 = b1_vec(ind_b1);
    
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
            
            v_x1(ind_b1,ind_b2,ind_d) = (1/pi)*integral(s1,0,Inf,'AbsTol',1e-6,'RelTol',1e-3);
            
            v_x2(ind_b1,ind_b2,ind_d) = (1/pi)*integral(s2,0,Inf,'AbsTol',1e-6,'RelTol',1e-3);
            
            acf1_x1(ind_b1,ind_b2,ind_d) = (1/pi)*integral(t1,0,Inf,'AbsTol',1e-6,'RelTol',1e-3)./v_x1(ind_b1,ind_b2,ind_d);
            
            acf1_x2(ind_b1,ind_b2,ind_d) = (1/pi)*integral(t2,0,Inf,'AbsTol',1e-6,'RelTol',1e-3)./v_x2(ind_b1,ind_b2,ind_d);
            
            cv_x1(ind_b1,ind_b2,ind_d) = sqrt(v_x1(ind_b1,ind_b2,ind_d))/m1;
            
            cv_x2(ind_b1,ind_b2,ind_d) = sqrt(v_x2(ind_b1,ind_b2,ind_d))/m2;
        end
        
    end
    
end

% figure('WindowState','fullscreen');
lw = 2;

subplot(3,3,1);
hold all
plot(b1_vec, acf1_x2(:,1,1),'k-','LineWidth',lw);
plot(b1_vec, acf1_x2(:,2,1),'b-','LineWidth',lw);
plot(b1_vec, acf1_x2(:,3,1),'r-','LineWidth',lw);
t = title('Low dispersal','FontWeight','bold','FontSize',14);
ylabel('Lag-1 Autocorrelation','FontSize',10);
ylim([min(squeeze(acf1_x2(:))) max(squeeze(acf1_x2(:)))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];

subplot(3,3,2);
hold all
plot(b1_vec, acf1_x2(:,1,2),'k-','LineWidth',lw);
plot(b1_vec, acf1_x2(:,2,2),'b-','LineWidth',lw);
plot(b1_vec, acf1_x2(:,3,2),'r-','LineWidth',lw);
title('Moderate dispersal','FontWeight','bold','FontSize',14);
ylim([min(squeeze(acf1_x2(:))) max(squeeze(acf1_x2(:)))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];

subplot(3,3,3);
hold all
plot(b1_vec, acf1_x2(:,1,3),'k-','LineWidth',lw);
plot(b1_vec, acf1_x2(:,2,3),'b-','LineWidth',lw);
plot(b1_vec, acf1_x2(:,3,3),'r-','LineWidth',lw);
title('High dispersal','FontWeight','bold','FontSize',14);
ylim([min(squeeze(acf1_x2(:))) max(squeeze(acf1_x2(:)))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];

subplot(3,3,4);
hold all
plot(b1_vec, v_x2(:,1,1),'k-','LineWidth',lw);
plot(b1_vec, v_x2(:,2,1),'b-','LineWidth',lw);
plot(b1_vec, v_x2(:,3,1),'r-','LineWidth',lw);
ylabel('Variance','FontSize',10);
ylim([min(squeeze(v_x2(:))) max(squeeze(v_x2(:)))]);
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.3f');
ax.YTick = [0.002 0.004 0.006];

set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
hold off

subplot(3,3,5);
hold all
plot(b1_vec, v_x2(:,1,2),'k-','LineWidth',lw);
plot(b1_vec, v_x2(:,2,2),'b-','LineWidth',lw);
plot(b1_vec, v_x2(:,3,2),'r-','LineWidth',lw);
ylim([min(squeeze(v_x2(:))) max(squeeze(v_x2(:)))]);
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.3f');
ax.YTick = [0.002 0.004 0.006];

set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
hold off

subplot(3,3,6);
hold all
h1 = plot(b1_vec, v_x2(:,1,3),'k-','LineWidth',lw);
h2 = plot(b1_vec, v_x2(:,2,3),'b-','LineWidth',lw);
h3 = plot(b1_vec, v_x2(:,3,3),'r-','LineWidth',lw);
s = {' \beta_2 = 0.01',' \beta_2 = 0.50',' \beta_2 = 0.99'};
f = legend([h1,h2,h3],s,'FontSize',7);
set(f,'Location','Best','box','off','FontSize',10);
ylim([min(squeeze(v_x2(:))) max(squeeze(v_x2(:)))]);
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.3f');
ax.YTick = [0.002 0.004 0.006];

set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];
hold off

subplot(3,3,7);
hold all
plot(b1_vec, cv_x2(:,1,1),'k-','LineWidth',lw);
plot(b1_vec, cv_x2(:,2,1),'b-','LineWidth',lw);
plot(b1_vec, cv_x2(:,3,1),'r-','LineWidth',lw);
xlabel('{\boldmath$\beta_1$}','FontWeight','bold','interpreter','latex');
ylabel('CV','FontSize',10);
ylim([min(squeeze(cv_x2(:))) max(squeeze(cv_x2(:)))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];

subplot(3,3,8);
hold all
plot(b1_vec, cv_x2(:,1,2),'k-','LineWidth',lw);
plot(b1_vec, cv_x2(:,2,2),'b-','LineWidth',lw);
plot(b1_vec, cv_x2(:,3,2),'r-','LineWidth',lw);
xlabel('{\boldmath$\beta_1$}','FontWeight','bold','interpreter','latex');
ylim([min(squeeze(cv_x2(:))) max(squeeze(cv_x2(:)))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];

subplot(3,3,9);
hold all
plot(b1_vec, cv_x2(:,1,3),'k-','LineWidth',lw);
plot(b1_vec, cv_x2(:,2,3),'b-','LineWidth',lw);
plot(b1_vec, cv_x2(:,3,3),'r-','LineWidth',lw);
xlabel('{\boldmath$\beta_1$}','FontWeight','bold','interpreter','latex');
ylim([min(squeeze(cv_x2(:))) max(squeeze(cv_x2(:)))]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw,'YMinorTick','Off'); 
ax.XAxis.TickLength = [0.025 0.0375];
ax.YAxis.TickLength = [0.025 0.0375];

export_fig Fig10 -pdf -painters -transparent