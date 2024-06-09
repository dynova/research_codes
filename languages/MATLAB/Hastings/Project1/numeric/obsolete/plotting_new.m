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
arr_thresh = [0.9 0.95 0.99];
len_thresh = length(arr_thresh);
proportion = zeros(len_d,len_thresh,len_b1,len_b2);
cell_d = {'low', 'medium', 'high'};
fs = 11;

for ind_thresh = 1:len_thresh
    for ind_d = 1:len_d
        for ind_b1 = 1:len_b1
            for ind_b2 = 1:len_b2
                thresh = arr_thresh(ind_thresh);
                idx = squeeze(prob(ind_d,ind_b1,ind_b2,:,:)) > thresh;
                proportion(ind_d,ind_thresh,ind_b1,ind_b2) = sum(idx(:))/45;
            end
        end
    end
end

% ind_cases = 0;
% figure(3);
% set(gcf, 'WindowState', 'maximized');
% 
% for ind_d = 1:len_d
%     d = arr_d(ind_d);
%     sgtitle('Proportion of parameter space (\nu) such that r > \eta','FontSize',18,'FontWeight','bold');
%     for ind_thresh = 1:len_thresh
%         thresh = arr_thresh(ind_thresh);
%         ind_cases = ind_cases + 1;
%         subplot(len_d,len_thresh,ind_cases);
%         xvalues = {'D + 0.01', 'D + 0.5', 'D + 0.99'};
%         yvalues = {'D + 0.01', 'D + 0.5', 'D + 0.99'};
%         h = heatmap(xvalues,yvalues,squeeze(proportion(ind_d,ind_thresh,:,:)),...
%             'MissingDataLabel','','MissingDataColor','white',...
%             'CellLabelFormat','%0.2f','CellLabelColor','white',...
%             'FontSize',fs,'Colormap',autumn,'GridVisible','off',...
%             'ColorLimits',[0 1],'ColorbarVisible','on');
%         h.Title = strcat(['D = ', num2str(d), ', \eta = ', num2str(thresh)]);
%         h.XLabel = '\beta_1';
%         h.YLabel = '\beta_2';
%     end
% end
% print(figure(3),'Fig3.png','-dpng','-r300');
% 
% for ind_d = len_d:len_d
%     d = arr_d(ind_d);
%     ind_cases = 0;
%     figure(4);
%     set(gcf, 'WindowState', 'maximized');
%     sgtitle(['Probability of system recovery (r) with ' cell_d{ind_d} ' dispersal (D = ' num2str(d) ')'],'FontSize',18,'FontWeight','bold');
%     for ind_b1 = 1:len_b1
%         b1 = d + arr_b1(ind_b1);
%         for ind_b2 = 1:len_b2
%             b2 = d + arr_b2(ind_b2);
%             ind_cases = ind_cases + 1;
%             subplot(len_b1,len_b2,ind_cases);
%             xvalues = {'2','3','4','5','6','7','8','9','10'};
%             yvalues = {'1','2','3','4','5','6','7','8','9'};
%             h = heatmap(fliplr(xvalues),fliplr(yvalues),...
%                 squeeze(prob(ind_d,ind_b1,ind_b2,:,:)),...
%                 'MissingDataLabel','','MissingDataColor','white',...
%                 'CellLabelFormat','%0.2f','CellLabelColor','black',...
%                 'FontSize',fs,'Colormap',autumn,'GridVisible','off',...
%                 'ColorLimits',[0 1]);
%             h.Title = strcat(['\beta_1 = D + ', num2str(b1-d), ', \beta_2 = D + ', num2str(b2-d)]);
%             h.XLabel = 'High Threshold (H)';
%             h.YLabel = 'Low Threshold (L)';
%         end
%     end
%     print(figure(4),'Fig4.png','-dpng','-r300');
% end
% 
% for ind_d = len_d:len_d
%     d = arr_d(ind_d);
%     ind_cases = 0;
%     figure(5);
%     set(gcf, 'WindowState', 'maximized');
%     sgtitle(['Odds of system recovery (r: 1-r) with ' cell_d{ind_d} ' dispersal (D = ' num2str(d) ')'],'FontSize',18,'FontWeight','bold');
%     for ind_b1 = 1:len_b1
%         b1 = d + arr_b1(ind_b1);
%         for ind_b2 = 1:len_b2
%             b2 = d + arr_b2(ind_b2);
%             ind_cases = ind_cases + 1;
%             subplot(len_b1,len_b2,ind_cases);
%             xvalues = {'2','3','4','5','6','7','8','9','10'};
%             yvalues = {'1','2','3','4','5','6','7','8','9'};
%             h = heatmap(fliplr(xvalues),fliplr(yvalues),...
%                 squeeze(prob(ind_d,ind_b1,ind_b2,:,:))./(1-squeeze(prob(ind_d,ind_b1,ind_b2,:,:))),...
%                 'CellLabelColor','none','Colormap',hot,'GridVisible','off',...
%                 'MissingDataLabel','','MissingDataColor','white',...
%                 'FontSize',fs,'ColorBarVisible','on','ColorScaling','Log',...
%                 'ColorLimits',[log(5e-1) log(5e4)]);
%             h.Title = strcat(['\beta_1 = D + ', num2str(b1-d), ', \beta_2 = D + ', num2str(b2-d)]);
%             h.XLabel = 'High Threshold (H)';
%             h.YLabel = 'Low Threshold (L)';
%             ax = struct(gca); % ignore warning that this should be avoided
%             cb = ax.Colorbar;
%             cb.TickLabels = {'1','100','10000'};
%         end
%     end
%     print(figure(5),'Fig5.png','-dpng','-r300');
% end