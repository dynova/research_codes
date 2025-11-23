% rng shuffle;
% num_sim = 5e2;
% time_max = 1.1e3;
% time_crit = 1e3;
% num_steps = 1e4;
% b_init = 0.99;
% b_0 = 1e-3;
% sigma_mu = 5e-2;
% num_cases = 3;
% [x1_init, x2_init] = deal(1 + sqrt(b_init));
% [x1, x2] = deal(zeros(num_cases,num_sim,num_steps+1));
% x1(:,:,1) = x1_init;
% x2(:,:,1) = x2_init;
% h = time_max/num_steps;
% hs = sqrt(h);
% z1 = randn(num_cases,num_sim,num_steps);
% z2 = randn(num_cases,num_sim,num_steps);
% d_vec = [0 0.01 1];
% time_vec = 0:h:time_max;
% b = arrayfun(@(t) b_init - b_0*t, time_vec);
% 
% for ii = 1:num_cases
%     d = d_vec(ii);
%     
%     for jj = 1:num_sim
%         
%         for kk = 1:num_steps
%             
%             if (x1(ii,jj,kk) <= 0 && x2(ii,jj,kk) > 0)
%                 V11 = [sigma_mu^2*x1(ii,jj,kk)^2 sigma_mu^2*x1(ii,jj,kk)^2 sigma_mu^2*x1(ii,jj,kk)^2];
%                 V12 = 0;
%                 V21 = 0;
%                 V22 = [sigma_mu^2*x2(ii,jj,kk)^2 sigma_mu^2*x2(ii,jj,kk)^2 sigma_mu^2*x2(ii,jj,kk)^2];
%                 B11 = sqrt(V11);
%                 B12 = 0;
%                 B21 = 0;
%                 B22 = sqrt(V22);
%                 x1_mean = d*x2(ii,jj,kk);
%                 x2_mean = b(kk)*x2(ii,jj,kk) - x2(ii,jj,kk) - d*x2(ii,jj,kk)...
%                     + 2*x2(ii,jj,kk)*x2(ii,jj,kk) - x2(ii,jj,kk)*x2(ii,jj,kk)*x2(ii,jj,kk);
%                 tmp = x1_mean*h + B11(ii)*hs*z1(ii,jj,kk) + B12*hs*z2(ii,jj,kk);
%                 x1(ii,jj,kk+1) = max(0,tmp);
%                 x2(ii,jj,kk+1) = x2(ii,jj,kk) + x2_mean*h + B21*hs*z1(ii,jj,kk) + B22(ii)*hs*z2(ii,jj,kk);
%                 
%             elseif (x1(ii,jj,kk) <= 0 && x2(ii,jj,kk) <= 0)
%                 x1(ii,jj,kk+1) = 0;
%                 x2(ii,jj,kk+1) = 0;
%                 
%             elseif (x2(ii,jj,kk) <= 0 && x1(ii,jj,kk) > 0)
%                 V11 = [sigma_mu^2*x1(ii,jj,kk)^2 sigma_mu^2*x1(ii,jj,kk)^2 sigma_mu^2*x1(ii,jj,kk)^2];
%                 V12 = 0;
%                 V21 = 0;
%                 V22 = [sigma_mu^2*x2(ii,jj,kk)^2 sigma_mu^2*x2(ii,jj,kk)^2 sigma_mu^2*x2(ii,jj,kk)^2];
%                 B11 = sqrt(V11);
%                 B12 = 0;
%                 B21 = 0;
%                 B22 = sqrt(V22);
%                 x1_mean = b(kk)*x1(ii,jj,kk) - x1(ii,jj,kk) - d*x1(ii,jj,kk)...
%                     + 2*x1(ii,jj,kk)*x1(ii,jj,kk) - x1(ii,jj,kk)*x1(ii,jj,kk)*x1(ii,jj,kk);
%                 x2_mean = d*x1(ii,jj,kk);
%                 tmp = x2_mean*h + B21*hs*z1(ii,jj,kk) + B22(ii)*hs*z2(ii,jj,kk);
%                 x1(ii,jj,kk+1) = x1(ii,jj,kk) + x1_mean*h + B11(ii)*hs*z1(ii,jj,kk) + B12*hs*z2(ii,jj,kk);
%                 x2(ii,jj,kk+1) = max(0,tmp);
%                 
%             elseif (x1(ii,jj,kk) > 0 && x2(ii,jj,kk) > 0)
%                 V11 = [sigma_mu^2*x1(ii,jj,kk)^2 sigma_mu^2*x1(ii,jj,kk)^2 sigma_mu^2*x1(ii,jj,kk)^2];
%                 V12 = 0;
%                 V21 = 0;
%                 V22 = [sigma_mu^2*x2(ii,jj,kk)^2 sigma_mu^2*x2(ii,jj,kk)^2 sigma_mu^2*x2(ii,jj,kk)^2];
%                 B11 = sqrt(V11);
%                 B12 = 0;
%                 B21 = 0;
%                 B22 = sqrt(V22);
%                 x1_mean = b(kk)*x1(ii,jj,kk) - x1(ii,jj,kk) - d*x1(ii,jj,kk) + d*x2(ii,jj,kk)...
%                     + 2*x1(ii,jj,kk)*x1(ii,jj,kk) - x1(ii,jj,kk)*x1(ii,jj,kk)*x1(ii,jj,kk);
%                 x2_mean = b(kk)*x2(ii,jj,kk) - x2(ii,jj,kk) - d*x2(ii,jj,kk) + d*x1(ii,jj,kk)...
%                     + 2*x2(ii,jj,kk)*x2(ii,jj,kk) - x2(ii,jj,kk)*x2(ii,jj,kk)*x2(ii,jj,kk);
%                 tmp_1 = x1(ii,jj,kk) + x1_mean*h + B11(ii)*hs*z1(ii,jj,kk) + B12*hs*z2(ii,jj,kk);
%                 tmp_2 = x2(ii,jj,kk) + x2_mean*h + B21*hs*z1(ii,jj,kk) + B22(ii)*hs*z2(ii,jj,kk);
%                 x1(ii,jj,kk+1) = max(0,tmp_1);
%                 x2(ii,jj,kk+1) = max(0,tmp_2);
%                 
%             end
%             
%         end
%         
%     end
%     
% end
% 
% save Fig1.mat -v7.3;

load Fig1.mat;
lw = 2;

ind_cases = 1;
subplot(2,3,ind_cases);
hold all
plot(time_vec, squeeze(x1(ind_cases,1:ceil(num_sim/10),:)),'Color',[0.5 0.5 0.5]);
plot(time_vec, squeeze(x1(ind_cases,1,:)),'k','LineWidth',lw);
plot(time_vec, mean(squeeze(x1(ind_cases,:,:))),'r','LineWidth',lw);
axis([0 time_max 0 1.2*x1_init]);
xline(time_crit,'k--');
xlabel('Time');
title('Isolated','FontWeight','bold','FontSize',14);
str = {'Population x_1,', 'Multiplicative noise'};
ylabel([str{1} newline str{2}]);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 2;
subplot(2,3,ind_cases);
hold all
plot(time_vec, squeeze(x1(ind_cases,1:ceil(num_sim/10),:)),'Color',[0.5 0.5 0.5]);
plot(time_vec, squeeze(x1(ind_cases,1,:)),'k','LineWidth',lw);
plot(time_vec, mean(squeeze(x1(ind_cases,:,:))),'r','LineWidth',lw);
axis([0 time_max 0 1.2*x1_init]);
xline(time_crit,'k--');
xlabel('Time');
title('Low dispersal','FontWeight','bold','FontSize',14);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

ind_cases = 3;
subplot(2,3,ind_cases);
hold all
plot(time_vec, squeeze(x1(ind_cases,1:ceil(num_sim/10),:)),'Color',[0.5 0.5 0.5]);
plot(time_vec, squeeze(x1(ind_cases,1,:)),'k','LineWidth',lw);
plot(time_vec, mean(squeeze(x1(ind_cases,:,:))),'r','LineWidth',lw);
axis([0 time_max 0 1.2*x1_init]);
xline(time_crit,'k--');
xlabel('Time');
title('High dispersal','FontWeight','bold','FontSize',14);
hold off

ax = gca;
set(ax,'FontWeight','bold','LineWidth',lw);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];

print -dpng -r1200 Fig1.png
RemoveWhiteSpace([], 'file', 'Fig1.png');