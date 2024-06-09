s1 = {'Line: three locations','Triangle: three locations'};
s2 = {'NLineC3','NTriangleC3'};
fig = figure(4);
set(gcf,'WindowState','fullscreen');
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
for ii = 1:length(s2)
    load(strcat(s2{ii},'.mat'));
    v = val;
    nexttile;
    h = heatmap(v','Colormap',gray,'GridVisible','off');
    ax = gca;
    ax.FontSize = 30;
    h.CellLabelFormat = '%.3f';
    ax.XDisplayLabels = nan(size(ax.XDisplayData));
    ax.YDisplayLabels = nan(size(ax.YDisplayData));
    h.Title = s1{ii};
    if ii == 1
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = [744.5 745.9 747.3];
        cb.TickLabels = {'744.5', '745.9', '747.3'};
    elseif ii == 2
        axs = struct(gca);
        cb = axs.Colorbar;
        cb.Ticks = [744.5 746.5 748.5];
        cb.TickLabels = {'744.5', '746.5', '748.5'};
    end
end
export_fig p3f5 -pdf -painters -transparent