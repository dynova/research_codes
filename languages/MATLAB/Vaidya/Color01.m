global Lambda mu mu_A mu_T beta tau eta_1 eta_2 lambda_T rho_1 rho_2 gamma
%% Pre-defined variables
alpha = 51;
TIME = 365;
Lambda = 531062;
mu = (1/70)/365;
mu_A = 0.25/365;
mu_T = 0.2/365;
beta = 0.187/365;
tau = 4/365;
lambda_T = 0.1;
rho_1 = 1/180;
rho_2 = (rho_1)./(270.*rho_1-1);
gamma = 1e-3;
%% Pre-allocation of output variables
eta_11 = zeros(1,alpha);
eta_22 = zeros(1,alpha);
AA = zeros(1,alpha);
BB = zeros(1,alpha);
CC = zeros(1,alpha);
Burden = zeros(1,alpha);
mseries = 1:alpha;
nseries = 1:alpha;
%% Begin loop code
for n = 1:alpha
    eta_1 = (n-1)./(alpha-1);
    eta_11(n) = (n-1)./(alpha-1);
    
    for m = 1:alpha
        eta_2 = (m-1)./(alpha-1);
        eta_22(m) = (m-1)./(alpha-1);
        
        %%%== Initial conditions ==
        S0 = 191564208;
        T0 = 131533276;
        H0 = 2405659;
        C0 = 1805024;
        C10 = 1000000;
        C20 = 1000000;
        CM10 = 500000;
        CM20 = 500000;
        y0 = [S0,T0,H0,C0,C10,C20,CM10,CM20];
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
        AA(n,m)=trapz(HIVinf1,HIVinf2);
        
        IRIScases1=0:TIME;
        IRIScases2=gamma.*(CM1+(rho_1/(rho_1+rho_2)).*(CM2));
        BB(n,m)=trapz(IRIScases1,IRIScases2);
        
        AIDSdeaths1=0:TIME;
        AIDSdeaths2=mu_A.*(C1+C2);
        CC(n,m)=trapz(AIDSdeaths1,AIDSdeaths2);
        
        Burden(n,m)=AA(n,m)+BB(n,m)+CC(n,m);
    end
end

figure(1);
set(gcf,'renderer','zbuffer');
[X,Y]=meshgrid(eta_11,eta_22);
surf(Y,X,AA./1e6,'Edgecolor','none','facecolor','interp');
view(2);
axis square
xlabel('$\eta_1$','Interpreter','LaTex','FontSize',16);
ylabel('$\eta_2$','Interpreter','LaTex','FontSize',16);
colorbar;

figure(2);
set(gcf,'renderer','zbuffer');
surf(Y,X,BB./1e6,'Edgecolor','none','facecolor','interp');
view(2);
axis square
xlabel('$\eta_1$','Interpreter','LaTex','FontSize',16);
ylabel('$\eta_2$','Interpreter','LaTex','FontSize',16);
colorbar;

figure(3);
set(gcf,'renderer','zbuffer');
surf(Y,X,CC./1e6,'Edgecolor','none','facecolor','interp');
view(2);
axis square
xlabel('$\eta_1$','Interpreter','LaTex','FontSize',16);
ylabel('$\eta_2$','Interpreter','LaTex','FontSize',16);
colorbar;

figure(4);
set(gcf,'renderer','zbuffer');
surf(Y,X,Burden./1e6,'Edgecolor','none','facecolor','interp');
view(2);
axis square
xlabel('$\eta_1$','Interpreter','LaTex','FontSize',16);
ylabel('$\eta_2$','Interpreter','LaTex','FontSize',16);
colorbar;