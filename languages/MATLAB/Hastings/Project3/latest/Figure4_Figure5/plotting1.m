s1 = {'Line: one location (at end)','Line: two locations (at ends)','Line: one location (middle)','Line: two locations (end and middle)','Triangle: one location','Triangle: two locations'};
s2 = {'NLineC1e','NLineC2e','NLineC1m','NLineC2m','NTriangleC1','NTriangleC2'};
fig = figure(4);
set(gcf,'WindowState','fullscreen');
t = tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
for ii = 1:length(s2)
    load(strcat(s2{ii},'.mat'));
    v = val;
    nexttile;
    h = heatmap(v','Colormap',gray,'GridVisible','off');
    ax = gca;
    ax.FontSize = 28;
    h.CellLabelFormat = '%.3f';
    ax.XDisplayLabels = nan(size(ax.XDisplayData));
    ax.YDisplayLabels = nan(size(ax.YDisplayData));
    h.Title = s1{ii};
    if ii == 1
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.TickLabels = {'510.0', '511.0', '512.0'};
    end
    if ii == 2
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.TickLabels = {'745.0', '746.0', '747.0'};
    end
    if ii == 4
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.TickLabels = {'745.0', '746.0', '747.0'};
    end    
    if ii == 3
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = [471.1 471.8 472.5];
        cb.TickLabels = {'471.1', '471.8', '472.5'};
    end
    if ii == 5
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = [471.1 472.0 472.9];
        cb.TickLabels = {'471.1', '472.0', '472.9'};
    end
    if ii == 6
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = [745 746.5 748];
        cb.TickLabels = {'745.0', '746.5', '748.0'};
    end    
end
export_fig p3f4 -pdf -painters -transparent