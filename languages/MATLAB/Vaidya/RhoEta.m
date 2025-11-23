global Lambda mu mu_A mu_T beta tau eta_1 eta_2 lambda_T rho_1 rho_2 gamma

eta_22 = 0:0.005:0.05;
nn = length(eta_22);
FinalAA = zeros(nn,269);
FinalBB = zeros(nn,269);
FinalCC = zeros(nn,269);
FinalBurden = zeros(nn,269);
%% Pre-defined variables
Lambda = 531062;
mu = (1/70)/365;
mu_A = 0.25/365;
mu_T = 0.2/365;
beta = 0.187/365;
tau = 4/365;
eta_1 = 0;
lambda_T = 0.1;
gamma = 1e-3;

for kk = 1:nn
    eta_2 = eta_22(kk);

    t = 269;
    TIME = 365;
    AA = zeros(1,t);
    BB = zeros(1,t);
    CC = zeros(1,t);
    mseries = 1:t;
    
    for m = 1:t
        rho_1 = 1/(m+1);
        rho_2 = (rho_1)./(271.*rho_1-1);
        
        %%%== Initial conditions ==
        S0 = 191564208;
        T0 = 131533276;
        H0 = 2405659;
        C0 = 1805024;
        C10 = 1000000;
        C20 = 1000000;
        CM10 = 500000;
        CM20 = 500000;
        y0 = [S0, T0, H0, C0, C10, C20, CM10, CM20];
        [tf,y] = ode45('Equations',0:TIME,y0);
        S = y(:,1);
        T = y(:,2);
        H = y(:,3);
        C = y(:,4);
        C1 = y(:,5);
        C2 = y(:,6);
        CM1 = y(:,7);
        CM2 = y(:,8);
        N=S+T+H+C+C1+C2+CM1+CM2;
        
        HIVinf1=0:TIME;
        HIVinf2=beta.*(S+T).*(C1+C2)./N;
        AA(m)=trapz(HIVinf1,HIVinf2);
        
        IRIScases1=0:TIME;
        IRIScases2=gamma.*(CM1+(rho_1/(rho_1+rho_2)).*CM2);
        BB(m)=trapz(IRIScases1,IRIScases2);
        
        AIDSdeaths1=0:TIME;
        AIDSdeaths2=mu_A.*(C1+C2);
        CC(m)=trapz(AIDSdeaths1,AIDSdeaths2);
        
        FinalAA(kk,m)=AA(m);
        FinalBB(kk,m)=BB(m);
        FinalCC(kk,m)=CC(m);
        FinalBurden(kk,m)=AA(m)+BB(m)+CC(m);
    end
end

figure(1);
set(gcf,'renderer','zbuffer');
[X1,Y1]=meshgrid(mseries-1,eta_22);
surf(X1,Y1,FinalBurden./1e6,'Edgecolor','none','facecolor','interp');
view(2);
axis([0 max(mseries-1) 0 max(eta_22)])
xlabel('$1/\rho_1$','Interpreter','LaTex','FontSize',15);
ylabel('$\eta_2$','Interpreter','LaTex','FontSize',15);
colorbar