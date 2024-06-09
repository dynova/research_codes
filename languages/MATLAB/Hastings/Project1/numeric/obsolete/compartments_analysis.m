% state_prob = zeros(len_d,len_b1,len_b2,len_t,len_init,len_HL/2);
load compartments.mat
counter = 1;
names = {'Fig1.jpeg','Fig2.jpeg','Fig3.jpeg','Fig4.jpeg'};
len_b1 = 2;
len_b2 = 2;
len_d = 2;
lw = 3;
p = [1 5 2 6 3 7 4 8];

for ind_b1 = 1:len_b1
    for ind_b2 = 1:len_b2
        plt(counter) = figure('WindowState','fullscreen');
        counter = counter + 1;
        ii = 0;
        for ind_d = 1:len_d
            for ind_init = 1:len_init
                ii = ii + 1;
                subplot(len_init,len_d,p(ii))
                hold all
                ax = gca;
                set(ax,'FontSize',16,'FontWeight','bold','LineWidth',lw,...
                    'XMinorTick','off','YMinorTick','off',...
                    'XScale','log','YScale','linear');
                xlabel('time','FontWeight','bold','FontSize',20);
                ylabel('probability');
                plot(array_t,squeeze(state_prob(ind_d,ind_b1,ind_b2,:,...
                    ind_init,:)),'linewidth',3);
                s={' HH',' HL',' LH', ' LL'};
                legend(s,'Location','Best','FontSize',16,'box','off');
                box on
                hold off
            end
        end
        print(plt(counter-1), '-djpeg', names{counter-1});
    end
end