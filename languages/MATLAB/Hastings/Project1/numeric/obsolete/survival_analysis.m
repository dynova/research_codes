load survival.mat
load mfpts.mat
lw = 3;
str = {'(HH \rightarrow LH)', '(HL \rightarrow LL)', ...
    '(HH \rightarrow HL)', '(LH \rightarrow LL)'};
array_t = logspace(-2,6,1001);
for ii = 1:len_HL
    plt(ii) = figure('WindowState','fullscreen');
    for jj = 1:len_d
        subplot(1,len_d,jj)
        for hh = 1:len_b1
            for ll = 1:len_b2
                hold all
                ax = gca;
                set(ax,'FontSize',16,'FontWeight','bold','LineWidth',lw,...
                    'XMinorTick','off','YMinorTick','off',...
                    'XScale','log','YScale','linear');
                xlabel('time','FontWeight','bold','FontSize',20);
                ylabel('survival');
                plot(array_t,arrayfun(matlabFunction(surv{jj,hh,ll,ii}),array_t),'LineWidth',lw);
                s={' \beta_1 = 0.01, \beta_2 = 0.01', ' \beta_1 = 0.01, \beta_2 = 0.99',...
                    ' \beta_1 = 0.99, \beta_2 = 0.01',' \beta_1 = 0.99, \beta_2 = 0.99'};
                legend(s,'Location','Best','FontSize',16,'box','off');
                box on
                hold off
            end
        end
    end
    print(plt(ii), '-djpeg', ['Fig' num2str(ii) '_survival']);
end

for ii = 1:len_HL
    plt(ii) = figure('WindowState','fullscreen');
    for jj = 1:len_d
        subplot(1,len_d,jj)
        for hh = 1:len_b1
            for ll = 1:len_b2
                hold all
                ax = gca;
                set(ax,'FontSize',16,'FontWeight','bold','LineWidth',lw,...
                    'XMinorTick','off','YMinorTick','off',...
                    'XScale','log','YScale','linear');
                xlabel('time','FontWeight','bold','FontSize',20);
                ylabel('fptd');
                plot(array_t,arrayfun(matlabFunction(fptd{jj,hh,ll,ii}),array_t),'LineWidth',lw);
                s={' \beta_1 = 0.01, \beta_2 = 0.01', ' \beta_1 = 0.01, \beta_2 = 0.99',...
                    ' \beta_1 = 0.99, \beta_2 = 0.01',' \beta_1 = 0.99, \beta_2 = 0.99'};
                legend(s,'Location','Best','FontSize',16,'box','off');
                box on
                hold off
            end
        end
    end
    print(plt(ii), '-djpeg', ['Fig' num2str(ii) '_fptd']);
end