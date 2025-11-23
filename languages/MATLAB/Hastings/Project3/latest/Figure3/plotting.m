load NLineC3.mat;
v1 = val;
load NTriangleC3.mat;
v2 = val;
len_d = 3;
len_b = 3;
len_o = 100;
arr_d = linspace(0.01,0.99,len_d);
arr_b = linspace(0.01,0.99,len_b);
arr_o = linspace(0.01,1.00,len_o);
fs = 16;
lw = 3;
fig = figure(3);
set(gcf, 'WindowState', 'fullscreen');
ind_sub = 0;
str = {'Very strong Allee effect; low dispersal','Very strong Allee effect; moderate dispersal','Very strong Allee effect; high dispersal','Moderately strong Allee effect; low dispersal','Moderately strong Allee effect; moderate dispersal','Moderately strong Allee effect; high dispersal','Weakly strong Allee effect; low dispersal','Weakly strong Allee effect; moderate dispersal','Weakly strong Allee effect; high dispersal'};
for ind_d = 1:len_d
    for ind_b = 1:len_b
        ind_sub = ind_sub + 1;
        subplot(len_b,len_d,ind_sub);
        hold all;
        h1 = plot(arr_o,squeeze(v1(ind_d,ind_b,:)),'k*','Linewidth',lw);
        h2 = plot(arr_o,squeeze(v2(ind_d,ind_b,:)),'k-','Linewidth',lw);
        ax = gca;
        set(ax,'FontWeight','bold','LineWidth',lw,'FontSize',fs);
        ax.XAxis.TickLength = [0.025 0.025];
        ax.YAxis.TickLength = [0.025 0.025];
        xlim([0 1]);
        ylim([min([v1(:);v2(:)]),max([v1(:);v2(:)])]);
        if ind_sub >= 7
            xlabel('Probability of accurate detection','FontSize',fs);
        end
        if mod(ind_sub,3) == 1
            ylabel('Reward value','FontSize',fs);
        end
        if ind_sub == 6
            s = {' Line', ' Triangle'};
            f = legend([h1,h2],s);
            set(f,'Location','North','box','off','FontSize',fs);
        end
        title(str{ind_sub},'FontWeight','bold','LineWidth',lw,'FontSize',fs-1);
        hold off;
    end
%     sgtitle('Model outcomes when all locations are observable','FontSize',18,'FontWeight','bold');
end
export_fig p3f3 -pdf -painters -transparent