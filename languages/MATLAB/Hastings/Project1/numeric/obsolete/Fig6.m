load main.mat;
N = 10;
len_d = 3;
len_b1 = 3;
len_b2 = 3;
len_H = N-1;
len_L = N-1;
arr_d = linspace(0.01,0.99,len_d);
arr_b1 = linspace(0.01,0.99,len_b1);
arr_b2 = linspace(0.01,0.99,len_b2);
arr_H = linspace(2,N,len_H);
arr_L = linspace(1,N-1,len_L);
arr_thresh = linspace(0.500,1.000,501*20);
len_thresh = length(arr_thresh);
proportion = zeros(1,len_thresh);
% f = zeros(1,len_thresh);
% metric = zeros(1,len_thresh);
% tbl = zeros(len_thresh,5);

format bank;

[a,b] = min(prob(:));
[ind_d,ind_b1,ind_b2,ind_h,ind_l] = ind2sub(size(prob),b);
d = arr_d(ind_d); b1 = arr_b1(ind_b1); b2 = arr_b2(ind_b2); h = arr_H(ind_h); l = arr_L(ind_l);
disp([d,b1,b2,h,l]);
% 0.5, 0.01, 0.01, 10, 9, prob = 0.274924915494288

[c,d] = max(prob(:));
[ind_d,ind_b1,ind_b2,ind_h,ind_l] = ind2sub(size(prob),d);
d = arr_d(ind_d); b1 = arr_b1(ind_b1); b2 = arr_b2(ind_b2); h = arr_H(ind_h); l = arr_L(ind_l);
disp([d,b1,b2,h,l]);
% 0.99, 0.99, 0.99, 9, 1, prob = 0.999977685408689

format long;
for ind_thresh = 1:len_thresh
    thresh = arr_thresh(ind_thresh);
    idx = prob(:) > thresh;
    arr = find(idx);
    proportion(ind_thresh) = numel(arr)/numel(prob(~isnan(prob)));
%     [ind_d,ind_b1,ind_b2,ind_h] = ind2sub(size(prob),arr);
%     f(1,ind_thresh) = norm((arr_d(ind_d)+arr_b1(ind_b1)+arr_b2(ind_b2))*N + arr_h(ind_h),1);
%     [metric(ind_thresh),b] = min((arr_d(ind_d)+arr_b1(ind_b1) + ...
%         arr_b2(ind_b2))*N + arr_h(ind_h));
%     tbl(ind_thresh,:) = [arr_d(ind_d(b)), arr_b1(ind_b1(b)), ...
%         arr_b2(ind_b2(b)), arr_h(ind_h(b)), thresh];
end

figure(1)
hold all
plot(arr_thresh,proportion,'k','LineWidth',3)
xlabel('Probability of system recovery (r)','FontSize',18);
ylabel('Proportion of parameter space (\nu)','FontSize',18);
ax = gca;
set(ax,'FontWeight','bold','LineWidth',3,'FontSize',16);
ax.XAxis.TickLength = [0.025 0.025];
ax.YAxis.TickLength = [0.025 0.025];
xlim([0.5 1]);
ylim([0 1]);
hold off
print(gcf,'Fig6.png','-dpng','-r600');

% stat = arr_thresh+proportion-f/norm(f,1);
% 
% figure(2)
% hold all
% plot(arr_thresh,stat,'k','LineWidth',3)
% xlabel('\eta','FontSize',18);
% ylabel('\sigma','FontSize',18);
% ax = gca;
% set(ax,'FontWeight','bold','LineWidth',3,'FontSize',16);
% ax.XAxis.TickLength = [0.025 0.025];
% ax.YAxis.TickLength = [0.025 0.025];
% xlim([0.5 1]);
% hold off
% print(gcf,'Fig5.png','-dpng','-r300');