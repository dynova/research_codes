load main.mat;
n = 9;
len_d = 3;
len_b1 = 3;
len_b2 = 3;
arr_d = linspace(0.01,0.99,len_d);
arr_b1 = linspace(0.01,0.99,len_b1);
arr_b2 = linspace(0.01,0.99,len_b2);
arr_thresh = [0.75 0.76 0.77];
len_thresh = length(arr_thresh);
proportion = zeros(len_d,len_thresh,len_b1,len_b2);
fs = 16;

for ind_thresh = 1:len_thresh
    for ind_d = 1:len_d
        for ind_b1 = 1:len_b1
            for ind_b2 = 1:len_b2
                thresh = arr_thresh(ind_thresh);
                idx = squeeze(prob(ind_d,ind_b1,ind_b2,:,:,:)) > thresh;
                proportion(ind_d,ind_thresh,ind_b1,ind_b2) = sum(idx(:))/nchoosek(n,3);
            end
        end
    end
end

ind_cases = 0;
figure(3);
set(gcf, 'WindowState', 'maximized');

for ind_d = 1:len_d
    d = arr_d(ind_d);    
%     sgtitle('Proportion (\nu) of parameter space such that r > \eta','FontSize',18,'FontWeight','bold');
    for ind_thresh = 1:len_thresh
        thresh = arr_thresh(ind_thresh);
        ind_cases = ind_cases + 1;
        subplot(len_d,len_thresh,ind_cases);
        xvalues = {'D + 0.01', 'D + 0.5', 'D + 0.99'};
        yvalues = {'D + 0.01', 'D + 0.5', 'D + 0.99'};
        h = heatmap(xvalues,yvalues,squeeze(proportion(ind_d,ind_thresh,:,:)),...
            'MissingDataLabel','','MissingDataColor','white',...
            'CellLabelFormat','%0.2f','CellLabelColor','black',...
            'FontSize',fs,'Colormap',autumn,'GridVisible','off',...
            'ColorLimits',[0 1],'ColorbarVisible','on');
        h.Title = strcat(['D = ', num2str(d), ', \eta = ', num2str(thresh)]);
        h.XLabel = '\beta_1';
        h.YLabel = '\beta_2';
    end
end
export_fig p1f3 -pdf -painters -transparent