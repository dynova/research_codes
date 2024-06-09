opts = delimitedTextImportOptions("NumVariables", 10);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["usa_or_hhsregion", "week_ending", "variant", "share", "share_hi", "share_lo", "nchs_or_count_flag", "modeltype", "time_interval", "published_date"];
opts.VariableTypes = ["double", "datetime", "categorical", "double", "double", "double", "categorical", "categorical", "categorical", "datetime"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["variant", "nchs_or_count_flag", "modeltype", "time_interval"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "week_ending", "InputFormat", "MM/dd/yyyy hh:mm:ss aa");
opts = setvaropts(opts, "published_date", "InputFormat", "MM/dd/yyyy hh:mm:ss aa");
tbl = readtable("/Users/amallela/Documents/bmab/SARS-CoV-2_Variant_Proportions.csv", opts);
clear opts;

x = find(tbl.usa_or_hhsregion == 2 & tbl.variant == 'B.1.1.7');
y = sortrows(tbl(x,:),2);
temp1 = unique(y.week_ending);

for i = 3:length(temp1)
    temp2 = find(y.week_ending == temp1(i));
    alpha_prev_nyc(i-2) = mean(y.share(temp2));
end

x = find(tbl.usa_or_hhsregion == 6 & tbl.variant == 'B.1.1.7');
y = sortrows(tbl(x,:),2);
temp1 = unique(y.week_ending);

for i = 3:length(temp1)
    temp2 = find(y.week_ending == temp1(i));
    alpha_prev_dallas(i-2) = mean(y.share(temp2));
    alpha_prev_houston(i-2) = mean(y.share(temp2));
end

x = find(tbl.usa_or_hhsregion == 9 & tbl.variant == 'B.1.1.7');
y = sortrows(tbl(x,:),2);
temp1 = unique(y.week_ending);

for i = 3:length(temp1)
    temp2 = find(y.week_ending == temp1(i));
    alpha_prev_phoenix(i-2) = mean(y.share(temp2));
end

x = find(tbl.usa_or_hhsregion == 2 & tbl.variant == 'B.1.617.2');
y = sortrows(tbl(x,:),2);
temp1 = unique(y.week_ending);

for i = 3:42
    temp2 = find(y.week_ending == temp1(i));
    delta_prev_nyc(i-2) = mean(y.share(temp2));
end

x = find(tbl.usa_or_hhsregion == 6 & tbl.variant == 'B.1.617.2');
y = sortrows(tbl(x,:),2);
temp1 = unique(y.week_ending);

for i = 3:42
    temp2 = find(y.week_ending == temp1(i));
    delta_prev_dallas(i-2) = mean(y.share(temp2));
    delta_prev_houston(i-2) = mean(y.share(temp2));
end

x = find(tbl.usa_or_hhsregion == 9 & tbl.variant == 'B.1.617.2');
y = sortrows(tbl(x,:),2);
temp1 = unique(y.week_ending);

for i = 3:42
    temp2 = find(y.week_ending == temp1(i));
    delta_prev_phoenix(i-2) = mean(y.share(temp2));
end

save prevalence.mat alpha_prev_phoenix alpha_prev_houston alpha_prev_dallas alpha_prev_nyc delta_prev_phoenix delta_prev_houston delta_prev_dallas delta_prev_nyc;