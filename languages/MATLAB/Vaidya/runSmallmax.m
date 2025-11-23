%% Figs. 4a,4b,4c
figure(1)
Smallmax0(3.17,0,2.33,730500)
hold on
Smallmax1(3.17,1,2.33,730500)
legend('show',...
      {'u_1^* (A_3=0)', 'u_2^* (A_3=0)', 'u_1^* (A_3=1)', 'u_2^* (A_3=1)'},...
      'Location','Best','FontSize',14);
set(gca,'XTick',0:100:300,'FontSize',14,'FontWeight','Bold');
set(gca,'YTick',0:0.01:0.06,'FontSize',14,'FontWeight','Bold');

figure(2)
Smallmax0(3.17,0,2.33,146100)
hold on
Smallmax1(3.17,1,2.33,146100)
legend('show',...
      {'u_1^* (A_3=0)', 'u_2^* (A_3=0)', 'u_1^* (A_3=1)', 'u_2^* (A_3=1)'},...
      'Location','Best','FontSize',14);
set(gca,'XTick',0:100:300,'FontSize',14,'FontWeight','Bold');
set(gca,'YTick',0:0.01:0.06,'FontSize',14,'FontWeight','Bold');


figure(3)
Smallmax0(3.17,0,2.33,29220)
hold on
Smallmax1(3.17,1,2.33,29220)
legend('show',...
      {'u_1^* (A_3=0)', 'u_2^* (A_3=0)', 'u_1^* (A_3=1)', 'u_2^* (A_3=1)'},...
      'Location','Best','FontSize',14);
set(gca,'XTick',0:100:300,'FontSize',14,'FontWeight','Bold');
set(gca,'YTick',0:0.01:0.06,'FontSize',14,'FontWeight','Bold');


%% Figs. 5a,5b,5c
figure(4)
Smallmax0(0,1,0,730500)
hold on
Smallmax1(3.17,1,2.33,730500)
legend('show',...
      {'u_1^* (A_1=0, A_2=0)', 'u_2^* (A_1=0, A_2=0)', 'u_1^* (A_1=3.17, A_2=2.33)', 'u_2^* (A_1=3.17, A_2=2.33)'},...
      'Location','Best','FontSize',14);
set(gca,'XTick',0:100:300,'FontSize',14,'FontWeight','Bold');
set(gca,'YTick',0:0.01:0.06,'FontSize',14,'FontWeight','Bold');


figure(5)
Smallmax0(0,1,0,146100)
hold on
Smallmax1(3.17,1,2.33,146100)
legend('show',...
      {'u_1^* (A_1=0, A_2=0)', 'u_2^* (A_1=0, A_2=0)', 'u_1^* (A_1=3.17, A_2=2.33)', 'u_2^* (A_1=3.17, A_2=2.33)'},...
      'Location','Best','FontSize',14);
set(gca,'XTick',0:100:300,'FontSize',14,'FontWeight','Bold');
set(gca,'YTick',0:0.01:0.06,'FontSize',14,'FontWeight','Bold');

figure(6)
Smallmax0(0,1,0,29220)
hold on
Smallmax1(3.17,1,2.33,29220)
legend('show',...
      {'u_1^* (A_1=0, A_2=0)', 'u_2^* (A_1=0, A_2=0)', 'u_1^* (A_1=3.17, A_2=2.33)', 'u_2^* (A_1=3.17, A_2=2.33)'},...
      'Location','Best','FontSize',14);
set(gca,'XTick',0:100:300,'FontSize',14,'FontWeight','Bold');
set(gca,'YTick',0:0.01:0.06,'FontSize',14,'FontWeight','Bold');