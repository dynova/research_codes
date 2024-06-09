function Smallmax1(A1,A2,A3,A4)
%% Prelims
KK = 30;
MAX = 0.05;
delta = 1e-3;
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
Time = 365;
test = -1;
M = Time;
tvec = linspace(0,Time,M+1)';
k = 0;
x=zeros(M+1,8);
lambda=zeros(M+1,8);
u1=zeros(M+1,1);
u2=zeros(M+1,1);

%% while loop
while(test < 0 && k < KK)
    k = k + 1;
    disp(k)
    oldx = x;
    oldlambda = lambda;
    oldu1 = u1;
    oldu2 = u2;    
    lambda5=lambda(:,5);
    lambda6=lambda(:,6);
    lambda7=lambda(:,7);
    lambda8=lambda(:,8);
      
    u11 = min(MAX,max(0,(x(:,5).*(lambda5-lambda7)./(2.*A4))));
    u22 = min(MAX,max(0,(x(:,6).*(lambda6-lambda8)./(2.*A4))));
    u1 = 0.5.*(u11 + oldu1);
    u2 = 0.5.*(u22 + oldu2);
    
    solx = ode45(@(t,x) states(t,x,tvec,u1,u2),tvec,[191564208 131533276 2405659 1805024 1000000 1000000 500000 500000]);
    x = deval(solx,tvec)';

    sollamb = ode45(@(t,lambda) adjoints(t,lambda,tvec,x,u1,u2),[Time 0],[0 0 0 0 0 0 0 0]);
    lambda = deval(sollamb,tvec)';
    
    test = min([delta*norm(oldu1,inf)-norm(oldu1-u1,inf) delta*norm(oldu2,inf)-norm(oldu2-u2,inf) delta*norm(oldx,inf)-norm(oldx-x,inf) delta*norm(oldlambda,inf)-norm(oldlambda-lambda,inf)]);
    disp(test)
    
%% Objective functional calculations
% HIVinf1=0:Time;
% HIVinf2=beta.*(x(:,1)+x(:,2)).*(x(:,5)+x(:,6))./(x(:,1)+x(:,2)+x(:,3)+x(:,4)+x(:,5)+x(:,6)+x(:,7)+x(:,8));
% TN=trapz(HIVinf1,HIVinf2);
% 
% HIVdeaths1=0:Time;
% HIVdeaths2=mu_A.*(x(:,5)+x(:,6));
% TD=trapz(HIVdeaths1,HIVdeaths2);
% 
% IRIScases1=0:Time;
% IRIScases2=gamma.*(x(:,7)+(rho_1./(rho_1+rho_2)).*(x(:,8)));
% TI=trapz(IRIScases1,IRIScases2);
% 
% Controls1=0:Time;
% Controls2=(u1).^2+(u2).^2;
% Cost=trapz(Controls1,Controls2);
% 
% Prod1 = A1.*TN;
% Prod2 = A2.*TI;
% Prod3 = A3.*TD;
% Prod4 = A4.*Cost;
% J = Prod1 + Prod2 + Prod3 + Prod4;
end
% save SmallMAX A1 TN A2 TI A3 TD A4 Cost Prod1 Prod2 Prod3 Prod4 J

%% Plot calls
plot(tvec,u1,'k','LineSmoothing','on');
axis([0 Time 0 1.2*MAX]);
xlabel('Time (Days)','FontSize',16,'FontWeight','Bold');
ylabel('Optimal controls','FontSize',16,'FontWeight','Bold');
hold on;
plot(tvec,u2,'b','LineSmoothing','on');

%% ===RECURSIVE CALLS===
%% STATES

function dx = states(t,x,tvec,u1,u2)
u1=pchip(tvec,u1,t);
u2=pchip(tvec,u2,t);
dx=zeros(8,1);
N = x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8);
E = x(3)+x(4)+x(5)+x(6);
dx(1) = Lambda-mu.*x(1)-beta.*E.*(x(1)./N)-tau.*(x(2)+x(4)).*(x(1)./N);
dx(2) = tau.*(x(2)+x(4)).*(x(1)./N)-beta.*E.*(x(2)./N)-(mu+mu_T).*x(2);
dx(3) = beta.*E.*(x(1)./N)-tau.*(x(2)+x(4)).*(x(3)./N)-(mu+mu_A).*x(3);
dx(4) = beta.*E.*(x(2)./N)+tau.*(x(2)+x(4)).*(x(3)./N)-(mu+mu_A+mu_T+lambda_T).*x(4);    
dx(5) = lambda_T.*x(4)-(mu+mu_A+rho_1+u1).*x(5);
dx(6) = rho_1.*x(5)-(mu+mu_A+rho_2+u2).*x(6);
dx(7) = u1.*x(5)-(mu+rho_1+gamma).*x(7);
dx(8) = u2.*x(6)-(rho_2+mu).*x(8)+(rho_1).*x(7)-gamma.*(rho_1/(rho_1+rho_2)).*x(8);
end

%% ADJOINTS
function dlambda = adjoints(t,lambda,tvec,x,u1,u2)
x=interp1(tvec,x,t);
u1=pchip(tvec,u1,t);
u2=pchip(tvec,u2,t);
dlambda=zeros(8,1);
N = x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8);
E = x(3)+x(4)+x(5)+x(6);

dlambda(1) = lambda(1)*(mu + (beta*(E))/(N) + (tau*(x(4) + x(2)))/(N) - (x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta*(E))/(N)^2) - lambda(2)*((tau*(x(4) + x(2)))/(N) - (x(1)*tau*(x(4) + x(2)))/(N)^2 + (x(2)*beta*(E))/(N)^2) + lambda(4)*((x(3)*tau*(x(4) + x(2)))/(N)^2 + (x(2)*beta*(E))/(N)^2) - lambda(3)*((beta*(E))/(N) + (x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta*(E))/(N)^2) - (A1*beta*(x(5) + x(6)))/(N) + (A1*beta*(x(1) + x(2))*(x(5) + x(6)))/(N)^2;
dlambda(2) = lambda(3)*((x(3)*tau)/(N) - (x(3)*tau*(x(4) + x(2)))/(N)^2 + (x(1)*beta*(E))/(N)^2) - lambda(4)*((x(3)*tau)/(N) + (beta*(E))/(N) - (x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) - lambda(1)*((x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*tau)/(N) + (x(1)*beta*(E))/(N)^2) + lambda(2)*(mu + mu_T - (x(1)*tau)/(N) + (beta*(E))/(N) + (x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) - (A1*beta*(x(5) + x(6)))/(N) + (A1*beta*(x(1) + x(2))*(x(5) + x(6)))/(N)^2;
dlambda(3) = lambda(3)*(mu + mu_A - (x(1)*beta)/(N) + (tau*(x(4) + x(2)))/(N) - (x(3)*tau*(x(4) + x(2)))/(N)^2 + (x(1)*beta*(E))/(N)^2) - lambda(4)*((x(2)*beta)/(N) + (tau*(x(4) + x(2)))/(N) - (x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) - lambda(1)*((x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta)/(N) + (x(1)*beta*(E))/(N)^2) + lambda(2)*((x(2)*beta)/(N) + (x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) + (A1*beta*(x(1) + x(2))*(x(5) + x(6)))/(N)^2;
dlambda(4) = lambda(2)*((x(2)*beta)/(N) - (x(1)*tau)/(N) + (x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) - lambda(5)*lambda_T + lambda(4)*(lambda_T + mu + mu_A + mu_T - (x(2)*beta)/(N) - (x(3)*tau)/(N) + (x(3)*tau*(x(4) + x(2)))/(N)^2 + (x(2)*beta*(E))/(N)^2) - lambda(3)*((x(1)*beta)/(N) - (x(3)*tau)/(N) + (x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta*(E))/(N)^2) + lambda(1)*((x(1)*beta)/(N) + (x(1)*tau)/(N) - (x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta*(E))/(N)^2) + (A1*beta*(x(1) + x(2))*(x(5) + x(6)))/(N)^2;
dlambda(5) = lambda(4)*((x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta)/(N) + (x(2)*beta*(E))/(N)^2) - lambda(6)*rho_1 - lambda(7)*u1 - lambda(3)*((x(1)*beta)/(N) + (x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta*(E))/(N)^2) - A3*mu_A - lambda(1)*((x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta)/(N) + (x(1)*beta*(E))/(N)^2) + lambda(2)*((x(2)*beta)/(N) + (x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) + lambda(5)*(mu + mu_A + rho_1 + u1) - (A1*beta*(x(1) + x(2)))/(N) + (A1*beta*(x(1) + x(2))*(x(5) + x(6)))/(N)^2;
dlambda(6) = lambda(4)*((x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta)/(N) + (x(2)*beta*(E))/(N)^2) - lambda(8)*u2 - lambda(3)*((x(1)*beta)/(N) + (x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta*(E))/(N)^2) - A3*mu_A - lambda(1)*((x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta)/(N) + (x(1)*beta*(E))/(N)^2) + lambda(2)*((x(2)*beta)/(N) + (x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) + lambda(6)*(mu + mu_A + rho_2 + u2) - (A1*beta*(x(1) + x(2)))/(N) + (A1*beta*(x(1) + x(2))*(x(5) + x(6)))/(N)^2;
dlambda(7) = lambda(7)*(gamma + mu + rho_1) - lambda(8)*rho_1 - A2*gamma - lambda(3)*((x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta*(E))/(N)^2) + lambda(4)*((x(3)*tau*(x(4) + x(2)))/(N)^2 + (x(2)*beta*(E))/(N)^2) - lambda(1)*((x(1)*tau*(x(4) + x(2)))/(N)^2 + (x(1)*beta*(E))/(N)^2) + lambda(2)*((x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) + (A1*beta*(x(1) + x(2))*(x(5) + x(6)))/(N)^2;
dlambda(8) = lambda(4)*((x(3)*tau*(x(4) + x(2)))/(N)^2 + (x(2)*beta*(E))/(N)^2) - lambda(3)*((x(3)*tau*(x(4) + x(2)))/(N)^2 - (x(1)*beta*(E))/(N)^2) - lambda(1)*((x(1)*tau*(x(4) + x(2)))/(N)^2 + (x(1)*beta*(E))/(N)^2) + lambda(2)*((x(1)*tau*(x(4) + x(2)))/(N)^2 - (x(2)*beta*(E))/(N)^2) + lambda(8)*(mu + rho_2 + (gamma*rho_1)/(rho_1 + rho_2)) - (A2*gamma*rho_1)/(rho_1 + rho_2) + (A1*beta*(x(1) + x(2))*(x(5) + x(6)))/(N)^2;
end
end