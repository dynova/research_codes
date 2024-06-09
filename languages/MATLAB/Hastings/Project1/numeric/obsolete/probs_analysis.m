load everything.mat;
[X,Y] = meshgrid(arr_b1/N, arr_b2/N);

figure(1)
for ind_d = 1:len_d
    hold all
    Z = squeeze(hh(ind_d,:,:));
    surf(X,Y,Z);
    zlim([0 1]);
    view(3)
    hold off
end

% counter = 1;
% s={' \beta_1 = 0.01, \beta_2 = 0.01';...
%     ' \beta_1 = 0.01, \beta_2 = 0.50';...
%     ' \beta_1 = 0.01, \beta_2 = 0.99';...
%     ' \beta_1 = 0.50, \beta_2 = 0.01';...
%     ' \beta_1 = 0.50, \beta_2 = 0.50';...
%     ' \beta_1 = 0.50, \beta_2 = 0.99';...
%     ' \beta_1 = 0.99, \beta_2 = 0.01';...
%     ' \beta_1 = 0.99, \beta_2 = 0.50';...
%     ' \beta_1 = 0.99, \beta_2 = 0.99'};
% for i = 1:len_b1
%     for j = 1:len_b2
%         figure(counter);
%         counter = counter + 1;
%         hold all
%         plot(arr_d/N,squeeze(hh(:,i,j)),'b','linewidth',3);
%         plot(arr_d/N,squeeze(ll(:,i,j)),'r','linewidth',3);
%         hold off
%         title(s{counter-1},'fontsize',18);        
%     end
% end