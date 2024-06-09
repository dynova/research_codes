load main.mat;
n = 9;
len_d = 3;
len_b1 = 3;
len_b2 = 3;
arr_d = linspace(0.01,0.99,len_d);
arr_b1 = linspace(0.01,0.99,len_b1);
arr_b2 = linspace(0.01,0.99,len_b2);
cell_d = {'low', 'medium', 'high'};
fs = 14;

figure(7);
set(gcf, 'WindowState', 'maximized');
    
for ind_d = 1
    d = arr_d(ind_d);
    ind_cases = 0;
%     sgtitle(['Probability of system recovery (r) with ' cell_d{ind_d} ' dispersal (D = ' num2str(d) ') and a high threshold of H = N'],'FontSize',18,'FontWeight','bold');
    for ind_b1 = 1:len_b1
        b1 = d + arr_b1(ind_b1);
        for ind_b2 = 1:len_b2
            b2 = d + arr_b2(ind_b2);
            ind_cases = ind_cases + 1;
            subplot(len_b1,len_b2,ind_cases);
            xvalues = {'1','2','3','4','5','6','7'};
            yvalues = {'2','3','4','5','6','7','8'};
            h = heatmap(xvalues,yvalues,...
                squeeze(prob(ind_d,ind_b1,ind_b2,2:n-1,n,1:n-2)),...
                'MissingDataLabel','','MissingDataColor','white',...
                'CellLabelFormat','%0.3f','CellLabelColor','black',...
                'FontSize',fs,'Colormap',autumn,'GridVisible','off',...
                'ColorLimits',[0.744 0.758],'ColorbarVisible','on');
            h.Title = strcat(['\beta_1 = D + ', num2str(b1-d), ', \beta_2 = D + ', num2str(b2-d)]);
            h.XLabel = 'Low Threshold (L)';
            h.YLabel = 'Allee Threshold (A)';
            axs = struct(gca);
            cb = axs.Colorbar;
            cb.TickLabels = {'0.745','0.750','0.755'};
        end
    end
end
export_fig p1f7 -pdf -painters -transparent